#include "blind_rotate.h"
#include "common/mod_arith.h"
#include "common/ntt.h"
#include "common/sampling.h"
#include <stdexcept>

namespace hehub {
namespace fhew {

namespace {
// coeff form 全零多项式。
RnsPolynomial zero_coeff(const tfhe::TfheParams &gp) {
    return get_zero_poly(RnsPolyParams{gp.N, 1, {gp.q}}, PolyRepForm::coeff);
}

// value(NTT) form 的单项式 X^exp (exp∈[0,2N), 负循环: X^N=-1)。
RnsPolynomial monomial_value(size_t exp, const tfhe::TfheParams &gp) {
    auto m = zero_coeff(gp);
    const size_t e = exp % (2 * gp.N);
    if (e < gp.N) {
        m[0][e] = 1;
    } else {
        m[0][e - gp.N] = gp.q - 1; // -1
    }
    ntt_negacyclic_inplace_lazy(m);
    return m;
}

// coeff form 单项式 X^exp (负循环), 供 RGSW 明文 μ 使用。
RnsPolynomial monomial_coeff(size_t exp, const tfhe::TfheParams &gp) {
    auto m = zero_coeff(gp);
    const size_t e = exp % (2 * gp.N);
    if (e < gp.N) {
        m[0][e] = 1;
    } else {
        m[0][e - gp.N] = gp.q - 1; // -1
    }
    return m;
}

// (-x) mod 2N, x 为非负整数(可能很大)。
size_t neg_mod_2n(unsigned long long x, size_t two_n) {
    size_t r = (size_t)(x % two_n);
    return (two_n - r) % two_n;
}
} // namespace

FhewParams make_fhew_params(size_t N, u64 q, size_t base_bits, size_t br_bits) {
    FhewParams fp;
    fp.gp = tfhe::make_tfhe_params(N, q, base_bits);
    fp.br_bits = br_bits;
    // d_r = ceil(log_{B_r}(2N)): 满足 B_r^{d_r} ≥ 2N, 使 a_i∈[0,2N) 的分解精确。
    const size_t two_n = 2 * N;
    size_t digits = 0;
    unsigned long long cap = 1;
    while (cap < two_n) {
        cap <<= br_bits;
        digits++;
    }
    fp.br_digits = digits;
    return fp;
}

FhewBootstrapKey gen_fhew_bootstrap_key(const LweSk &lwe_sk,
                                        const RlweSk &rlwe_sk,
                                        const FhewParams &fp) {
    const size_t n = lwe_sk.s.size();
    const size_t two_n = 2 * fp.gp.N;
    const size_t B_r = fp.base_r();

    FhewBootstrapKey bk;
    bk.params = fp;
    bk.n = n;
    bk.ek.assign(n, {});

    for (size_t i = 0; i < n; i++) {
        const int si = lwe_sk.s[i]; // 三元 {-1,0,1}(AP 法不限二元)
        bk.ek[i].assign(fp.br_digits, {});
        // B_r^j mod 2N (指数在负循环群 Z_{2N} 中周期化)。
        unsigned long long base_pow = 1; // B_r^0
        for (size_t j = 0; j < fp.br_digits; j++) {
            bk.ek[i][j].resize(B_r); // v=0 项留空(外积单位元, 盲旋转跳过)
            for (size_t v = 1; v < B_r; v++) {
                // 指数 = -(v·B_r^j·s_i) mod 2N。s_i 有符号, 用有符号折算。
                long long signed_exp =
                    (long long)((unsigned long long)v * base_pow %
                                (unsigned long long)two_n) *
                    (long long)si;
                size_t e;
                if (signed_exp >= 0) {
                    e = neg_mod_2n((unsigned long long)signed_exp, two_n);
                } else {
                    // -(负数) = 正数
                    e = (size_t)((unsigned long long)(-signed_exp) % two_n);
                }
                auto mu = monomial_coeff(e, fp.gp);
                bk.ek[i][j][v] = tfhe::tfhe_rgsw_encrypt(mu, rlwe_sk, fp.gp);
            }
            base_pow = base_pow * B_r % (unsigned long long)two_n;
        }
    }
    return bk;
}

RlweCt fhew_blind_rotate(const LweCt &ct_mod2N,
                         const RnsPolynomial &testvec_coeff,
                         const FhewBootstrapKey &bk) {
    const auto &fp = bk.params;
    const auto &gp = fp.gp;
    const size_t n = ct_mod2N.a.size();
    if (bk.n != n || bk.ek.size() != n) {
        throw std::invalid_argument("FHEW 自举密钥维度与 LWE 不符。");
    }
    const size_t two_n = 2 * gp.N;
    const size_t B_r = fp.base_r();
    const size_t mask = B_r - 1;

    // ACC = 平凡 RLWE(testvec) = (ntt(testvec), 0), 再乘 X^{-b}。
    RnsPolynomial acc0 = testvec_coeff;
    if (acc0.rep_form == PolyRepForm::coeff) {
        ntt_negacyclic_inplace_lazy(acc0);
    }
    RlweCt acc{std::move(acc0),
               get_zero_poly(RnsPolyParams{gp.N, 1, {gp.q}}, PolyRepForm::value)};
    {
        auto mono = monomial_value((two_n - ct_mod2N.b % two_n) % two_n, gp);
        acc = RlweCt{acc[0] * mono, acc[1] * mono};
    }

    // 逐维、逐数字外积: ACC ← ACC ⊡ ek[i][j][a_{i,j}]。
    // 累积效果 ACC ← ACC·Π X^{-(a_{i,j}·B_r^j·s_i)} = ACC·X^{-Σa_i s_i}。
    for (size_t i = 0; i < n; i++) {
        size_t ai = ct_mod2N.a[i] % two_n; // ∈[0,2N), 分解精确
        for (size_t j = 0; j < fp.br_digits; j++) {
            const size_t digit = ai & mask;
            ai >>= fp.br_bits;
            if (digit == 0) {
                continue; // ek[i][j][0] = RGSW(1), 单位元, 跳过
            }
            acc = tfhe::tfhe_external_product(acc, bk.ek[i][j][digit], gp);
        }
    }
    return acc;
}

LweCt fhew_functional_bootstrap(const LweCt &ct, const RnsPolynomial &lut_poly,
                                const FhewBootstrapKey &bk) {
    const size_t N = bk.params.gp.N;
    const LweCt ct_2n = tfhe::lwe_mod_switch(ct, 2 * N);
    const RlweCt acc = fhew_blind_rotate(ct_2n, lut_poly, bk);
    return sample_extract(acc, 0);
}

} // namespace fhew
} // namespace hehub
