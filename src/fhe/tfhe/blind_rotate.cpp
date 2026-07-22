#include "blind_rotate.h"
#include "common/mod_arith.h"
#include "common/ntt.h"
#include "common/sampling.h"
#include <stdexcept>

namespace hehub {
namespace tfhe {

namespace {
using u128 = unsigned __int128;

// value(NTT) form 的全零多项式。
RnsPolynomial zero_value(const TfheParams &gp) {
    return get_zero_poly(RnsPolyParams{gp.N, 1, {gp.q}}, PolyRepForm::value);
}

// coeff form 的全零多项式。
RnsPolynomial zero_coeff(const TfheParams &gp) {
    return get_zero_poly(RnsPolyParams{gp.N, 1, {gp.q}}, PolyRepForm::coeff);
}

// value(NTT) form 的单项式 X^exp (exp∈[0,2N), 负循环: X^N=-1)。
RnsPolynomial monomial_value(size_t exp, const TfheParams &gp) {
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

// 有符号 base-B gadget 分解: 输入 coeff form 单分量多项式(coeffs∈[0,q)),
// 逐系数分解为 ℓ 个有符号数字(∈[-B/2,B/2)), 各层构成一个 coeff form 多项式(mod q)。
std::vector<RnsPolynomial> gadget_decompose(const RnsPolynomial &poly_coeff,
                                            const TfheParams &gp) {
    const u64 q = gp.q;
    const size_t N = poly_coeff.dimension();
    const u64 B = (u64)1 << gp.base_bits;
    const u64 mask = B - 1;
    const int64_t half = (int64_t)(B >> 1);
    std::vector<RnsPolynomial> digits(gp.levels);
    for (auto &d : digits) {
        d = zero_coeff(gp);
    }
    for (size_t j = 0; j < N; j++) {
        u64 c = poly_coeff[0][j] % q;
        for (size_t l = 0; l < gp.levels; l++) {
            int64_t dig = (int64_t)(c & mask);
            c >>= gp.base_bits;
            if (dig >= half) {
                dig -= (int64_t)B;
                c += 1; // 进位
            }
            u64 v = (dig >= 0) ? ((u64)dig % q) : (q - ((u64)(-dig) % q));
            digits[l][0][j] = v;
        }
    }
    return digits;
}
} // namespace

TfheParams make_tfhe_params(size_t N, u64 q, size_t base_bits) {
    TfheParams gp;
    gp.N = N;
    gp.q = q;
    gp.base_bits = base_bits;
    // levels 满足 B^ℓ ≥ 2q, 即 ℓ ≥ (bits(q)+1)/base_bits。
    size_t qbits = 0;
    for (u64 t = q; t; t >>= 1) {
        qbits++;
    }
    gp.levels = (qbits + 1 + base_bits - 1) / base_bits;
    return gp;
}

RgswCt tfhe_rgsw_encrypt(const RnsPolynomial &mu_coeff, const RlweSk &sk,
                         const TfheParams &gp) {
    if (mu_coeff.rep_form != PolyRepForm::coeff) {
        throw std::invalid_argument("RGSW 明文 μ 需 coeff form。");
    }
    // ntt(μ) 与 ntt(μ·s): sk 已是 value form, 逐点乘即 μ·s 的 value form。
    RnsPolynomial mu_val = mu_coeff;
    ntt_negacyclic_inplace_lazy(mu_val);
    RnsPolynomial mus_val = mu_val * sk;

    RgswCt rgsw(2 * gp.levels);
    const u64 q = gp.q;
    const u64 B = (u64)1 << gp.base_bits;
    u64 Bl = 1 % q; // B^l mod q
    for (size_t l = 0; l < gp.levels; l++) {
        // RLWE'(μ) 第 l 行: 相位 e, 加 μ·B^l 到 c0。
        RlweCt row0 = get_rlwe_sample(sk);
        RnsPolynomial add0 = mu_val;
        add0 *= Bl; // 标量乘, value form
        row0[0] += add0;
        rgsw[l] = row0;

        // RLWE'(μ·s) 第 l 行: 相位 e, 加 μ·s·B^l 到 c0。
        RlweCt row1 = get_rlwe_sample(sk);
        RnsPolynomial add1 = mus_val;
        add1 *= Bl;
        row1[0] += add1;
        rgsw[gp.levels + l] = row1;

        Bl = (u64)((u128)Bl * B % q);
    }
    return rgsw;
}

RlweCt tfhe_external_product(const RlweCt &ct, const RgswCt &rgsw,
                            const TfheParams &gp) {
    if (rgsw.size() != 2 * gp.levels) {
        throw std::invalid_argument("RGSW 行数与 gadget levels 不符。");
    }
    // 转 coeff form 后 gadget 分解 c0, c1。
    RnsPolynomial c0 = ct[0], c1 = ct[1];
    intt_negacyclic_inplace_lazy(c0);
    reduce_strict(c0);
    intt_negacyclic_inplace_lazy(c1);
    reduce_strict(c1);
    auto D0 = gadget_decompose(c0, gp);
    auto D1 = gadget_decompose(c1, gp);
    for (auto &d : D0) {
        ntt_negacyclic_inplace_lazy(d);
    }
    for (auto &d : D1) {
        ntt_negacyclic_inplace_lazy(d);
    }

    auto out0 = zero_value(gp);
    auto out1 = zero_value(gp);
    for (size_t l = 0; l < gp.levels; l++) {
        // RLWE'(μ) 行 (配 c0 分解): 贡献 μ·c0。
        out0 += D0[l] * rgsw[l][0];
        out1 += D0[l] * rgsw[l][1];
        // RLWE'(μ·s) 行 (配 c1 分解): 贡献 μ·s·c1。
        out0 += D1[l] * rgsw[gp.levels + l][0];
        out1 += D1[l] * rgsw[gp.levels + l][1];
    }
    return RlweCt{std::move(out0), std::move(out1)};
}

RlweCt cmux(const RgswCt &bit_rgsw, const RlweCt &ct0, const RlweCt &ct1,
            const TfheParams &gp) {
    RlweCt diff{ct1[0] - ct0[0], ct1[1] - ct0[1]};
    RlweCt prod = tfhe_external_product(diff, bit_rgsw, gp);
    return RlweCt{ct0[0] + prod[0], ct0[1] + prod[1]};
}

std::vector<RgswCt> gen_bootstrap_keys(const LweSk &lwe_sk,
                                       const RlweSk &rlwe_sk,
                                       const TfheParams &gp) {
    std::vector<RgswCt> bk(lwe_sk.s.size());
    for (size_t i = 0; i < lwe_sk.s.size(); i++) {
        const int si = lwe_sk.s[i];
        if (si != 0 && si != 1) {
            throw std::invalid_argument("盲旋转要求二元 LWE 私钥 s_i∈{0,1}。");
        }
        auto mu = zero_coeff(gp);
        mu[0][0] = (u64)si; // 常数多项式 μ = s_i
        bk[i] = tfhe_rgsw_encrypt(mu, rlwe_sk, gp);
    }
    return bk;
}

LweCt lwe_mod_switch(const LweCt &ct, u64 new_mod) {
    LweCt out;
    out.modulus = new_mod;
    out.a.resize(ct.a.size());
    const u64 old = ct.modulus;
    auto sw = [&](u64 x) -> u64 {
        return (u64)(((u128)(x % old) * new_mod + old / 2) / old) % new_mod;
    };
    for (size_t i = 0; i < ct.a.size(); i++) {
        out.a[i] = sw(ct.a[i]);
    }
    out.b = sw(ct.b);
    return out;
}

RlweCt blind_rotate(const LweCt &ct_mod2N, const RnsPolynomial &testvec_coeff,
                    const std::vector<RgswCt> &bk, const TfheParams &gp) {
    const size_t n = ct_mod2N.a.size();
    if (bk.size() != n) {
        throw std::invalid_argument("自举密钥数量与 LWE 维度不符。");
    }
    const size_t two_n = 2 * gp.N;

    // ACC = 平凡 RLWE(testvec): (ntt(testvec), 0), 再乘 X^{-b}。
    RnsPolynomial acc0 = testvec_coeff;
    if (acc0.rep_form == PolyRepForm::coeff) {
        ntt_negacyclic_inplace_lazy(acc0);
    }
    RlweCt acc{std::move(acc0), zero_value(gp)};
    {
        auto mono = monomial_value((two_n - ct_mod2N.b % two_n) % two_n, gp);
        acc = RlweCt{acc[0] * mono, acc[1] * mono};
    }

    // 逐维 CMux: ACC ← CMux(bk[i], ACC, X^{-a_i}·ACC) = ACC·X^{-a_i·s_i}。
    for (size_t i = 0; i < n; i++) {
        auto mono = monomial_value((two_n - ct_mod2N.a[i] % two_n) % two_n, gp);
        RlweCt rot{acc[0] * mono, acc[1] * mono};
        acc = cmux(bk[i], acc, rot, gp);
    }
    return acc;
}

RlweCt tfhe_rlwe_encrypt(const RnsPolynomial &pt, u64 t, const RlweSk &sk,
                         const TfheParams &gp) {
    const u64 q = gp.q;
    const u64 delta = (u64)(((u128)q + t / 2) / t); // round(q/t)
    RnsPolynomial scaled = pt;
    if (scaled.rep_form != PolyRepForm::coeff) {
        throw std::invalid_argument("明文需 coeff form。");
    }
    scaled *= delta; // Δ·pt, coeff form
    ntt_negacyclic_inplace_lazy(scaled);
    auto row = get_rlwe_sample(sk);
    row[0] += scaled;
    return row;
}

RnsPolynomial tfhe_rlwe_decrypt(const RlweCt &ct, u64 t, const RlweSk &sk,
                                const TfheParams &gp) {
    RlwePt phase = decrypt_core(ct, sk); // coeff form, ∈[0,q)
    const u64 q = gp.q;
    auto out = zero_coeff(gp);
    for (size_t j = 0; j < gp.N; j++) {
        u64 x = phase[0][j] % q;
        // round(t·x/q) mod t
        u64 r = (u64)(((u128)x * t + q / 2) / q) % t;
        out[0][j] = r;
    }
    return out;
}

} // namespace tfhe
} // namespace hehub
