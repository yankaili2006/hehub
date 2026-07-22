#include "bfv.h"
#include "bfv_internal.h"
#include "common/mod_arith.h"
#include "common/ntt.h"
#include <stdexcept>
#include <vector>

namespace hehub {
namespace bfv {

namespace {
using i128 = __int128;
using u128 = unsigned __int128;

// 把 RnsPolynomial(coeff form) 逐系数 CRT 重建并中心化到 [-Q/2, Q/2) 的整数(i128)。
// 首版约定: Q 及卷积中间量在 i128 范围内(约 Q<2^58, 密文模个数少)。
std::vector<i128> reconstruct_centered(const RnsPolynomial &poly, const UBInt &Q,
                                       const UBInt &half_Q,
                                       const std::vector<UBInt> &crt_coeff) {
    const auto k = poly.component_count();
    const auto n = poly.dimension();
    std::vector<i128> out(n);
    for (size_t j = 0; j < n; j++) {
        UBInt x;
        for (size_t i = 0; i < k; i++) {
            x += UBInt(poly[i][j]) * crt_coeff[i];
        }
        x %= Q;
        if (x < half_Q) {
            out[j] = (i128)to_u64(x); // x < Q < 2^64
        } else {
            out[j] = -(i128)to_u64(Q - x);
        }
    }
    return out;
}

// 负循环卷积(mod x^n+1): r[k] = Σ_{i+j=k} a·b - Σ_{i+j=k+n} a·b。
std::vector<i128> negacyclic_conv(const std::vector<i128> &a,
                                  const std::vector<i128> &b) {
    const size_t n = a.size();
    std::vector<i128> r(n, 0);
    for (size_t i = 0; i < n; i++) {
        if (a[i] == 0)
            continue;
        for (size_t j = 0; j < n; j++) {
            i128 prod = a[i] * b[j];
            size_t k = i + j;
            if (k < n) {
                r[k] += prod;
            } else {
                r[k - n] -= prod; // x^n = -1
            }
        }
    }
    return r;
}

// 把中心化整数多项式按 t/Q 缩放取整, 再 reduce 到 RNS(coeff form, value=false)。
RnsPolynomial scale_by_t_over_Q(const std::vector<i128> &coeffs,
                                const std::vector<u64> &moduli, const UBInt &Q,
                                const UBInt &half_Q, u64 t) {
    const size_t n = coeffs.size();
    RnsPolynomial out(RnsPolyParams{n, moduli.size(), moduli});
    out.rep_form = PolyRepForm::coeff;
    const UBInt t_big(t);
    for (size_t j = 0; j < n; j++) {
        i128 x = coeffs[j];
        bool neg = x < 0;
        u128 mag = neg ? (u128)(-x) : (u128)x;
        // round(t·mag/Q) = (t·mag + Q/2) / Q
        UBInt r = (t_big * ubint_from_u128(mag) + half_Q) / Q;
        for (size_t i = 0; i < moduli.size(); i++) {
            const u64 q_i = moduli[i];
            u64 v = to_u64(r % UBInt(q_i));
            if (neg && v != 0) {
                v = q_i - v;
            }
            out[i][j] = v;
        }
    }
    return out;
}
} // namespace

BfvQuadraticCt mult_low_level(const BfvCt &ct1, const BfvCt &ct2) {
    if (ct1.plain_modulus != ct2.plain_modulus) {
        throw std::invalid_argument("Plain moduli mismatch.");
    }
    const auto moduli = ct1[0].modulus_vec();
    const u64 t = ct1.plain_modulus;
    const UBInt Q = product_of(moduli);
    const UBInt half_Q = Q / UBInt(2);
    std::vector<UBInt> crt_coeff(moduli.size());
    for (size_t i = 0; i < moduli.size(); i++) {
        const u64 q_i = moduli[i];
        const UBInt M_i = Q / UBInt(q_i);
        crt_coeff[i] =
            M_i * UBInt(inverse_mod_prime(to_u64(M_i % UBInt(q_i)), q_i));
    }

    // 转 coeff form 后中心化重建。
    auto to_int = [&](const RnsPolynomial &p_value) {
        RnsPolynomial p = p_value;
        intt_negacyclic_inplace_lazy(p);
        reduce_strict(p);
        return reconstruct_centered(p, Q, half_Q, crt_coeff);
    };
    auto c0 = to_int(ct1[0]), c1 = to_int(ct1[1]);
    auto d0 = to_int(ct2[0]), d1 = to_int(ct2[1]);

    // tensor: (b0,b1,b2) = (c0d0, c0d1+c1d0, c1d1)。
    auto P00 = negacyclic_conv(c0, d0);
    auto P01 = negacyclic_conv(c0, d1);
    auto P10 = negacyclic_conv(c1, d0);
    auto P11 = negacyclic_conv(c1, d1);
    const size_t n = c0.size();
    std::vector<i128> b1_int(n);
    for (size_t j = 0; j < n; j++) {
        b1_int[j] = P01[j] + P10[j];
    }

    // 各分量按 t/Q 缩放 + 转 value form。
    auto b0 = scale_by_t_over_Q(P00, moduli, Q, half_Q, t);
    auto b1 = scale_by_t_over_Q(b1_int, moduli, Q, half_Q, t);
    auto b2 = scale_by_t_over_Q(P11, moduli, Q, half_Q, t);
    ntt_negacyclic_inplace_lazy(b0);
    ntt_negacyclic_inplace_lazy(b1);
    ntt_negacyclic_inplace_lazy(b2);

    BfvQuadraticCt out;
    out[0] = std::move(b0);
    out[1] = std::move(b1);
    out[2] = std::move(b2);
    out.plain_modulus = t;
    return out;
}

BfvPt decrypt(const BfvQuadraticCt &ct, const RlweSk &rlwe_sk) {
    // phase = b0 + b1·s + b2·s²  (coeff form)。
    RlweCt lin{ct[0], ct[1]};
    auto phase = decrypt_core(lin, rlwe_sk); // b0 + b1·s, coeff form

    // sk 已是 value(NTT) form; s² = sk·sk (value form), b2·s² 再转 coeff 累加。
    auto s2 = rlwe_sk * rlwe_sk;
    auto b2s2 = ct[2] * s2;
    intt_negacyclic_inplace_lazy(b2s2);
    reduce_strict(b2s2);

    phase += b2s2;
    reduce_strict(phase);
    return scale_and_round(phase, ct.plain_modulus);
}

} // namespace bfv
} // namespace hehub
