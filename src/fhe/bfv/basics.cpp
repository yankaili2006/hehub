#include "bfv.h"
#include "bgv/bgv.h"
#include "common/bigint.h"
#include "common/mod_arith.h"
#include "common/ntt.h"
#include "common/rns_transform.h"
#include <algorithm>
#include <stdexcept>

namespace hehub {
namespace bfv {

// 打包编码/解码与 BGV 完全一致(纯明文 CRT 批处理, 与方案无关), 直接复用。
BfvPt simd_encode(const std::vector<u64> &data, const u64 modulus,
                  size_t slot_count) {
    return bgv::simd_encode(data, modulus, slot_count);
}

std::vector<u64> simd_decode(const BfvPt &pt, size_t data_size) {
    return bgv::simd_decode(pt, data_size);
}

namespace {
// 分量原地标量乘: data[j] = data[j] * scalar mod modulus。
inline void scale_component_inplace(u64 *data, size_t n, u64 scalar,
                                    u64 modulus) {
    for (size_t j = 0; j < n; j++) {
        data[j] = (u64)(((unsigned __int128)data[j] * scalar) % modulus);
    }
}

// Q = ∏ ct_moduli (大整数)。
UBInt product_of(const std::vector<u64> &moduli) {
    UBInt q(1);
    for (auto m : moduli) {
        q *= UBInt(m);
    }
    return q;
}
} // namespace

BfvCt encrypt(const BfvPt &pt, const RlweSk &rlwe_sk,
              std::vector<u64> ct_moduli) {
    const auto t = pt.modulus_at(0);
    if (ct_moduli.empty()) {
        ct_moduli = rlwe_sk.modulus_vec();
    }
    if (std::find(ct_moduli.begin(), ct_moduli.end(), t) != ct_moduli.end()) {
        throw std::logic_error(
            "Plaintext modulus t must be coprime with ciphertext moduli.");
    }

    // 零加密样本 {c0, c1}, 满足 c0 + c1·s = e (NTT value form)。
    auto [c0, c1] = get_rlwe_sample(rlwe_sk, ct_moduli.size());

    // 把明文(系数 <t)迁移到密文模 RNS 下(coeff form)。
    auto pt_under_ct = rns_base_transform(pt, ct_moduli);

    // Δ = ⌊Q/t⌋; 逐分量乘 (Δ mod q_i), 得 Δ·pt。
    const UBInt Q = product_of(ct_moduli);
    const UBInt Delta = Q / UBInt(t);
    for (size_t i = 0; i < ct_moduli.size(); i++) {
        const u64 q_i = ct_moduli[i];
        const u64 delta_i = to_u64(Delta % UBInt(q_i));
        scale_component_inplace(pt_under_ct[i].data(), pt_under_ct.dimension(),
                                delta_i, q_i);
    }

    // 转到 NTT value form 与 c0 相加。
    ntt_negacyclic_inplace_lazy(pt_under_ct);
    c0 += pt_under_ct;

    BfvCt ct = RlweCt{std::move(c0), std::move(c1)};
    ct.plain_modulus = t;
    return ct;
}

BfvPt decrypt(const BfvCt &ct, const RlweSk &rlwe_sk) {
    // phase = (c0 + c1·s) mod Q = Δ·pt + e, coeff form, RNS。
    auto phase = decrypt_core(ct, rlwe_sk);
    const auto moduli = phase.modulus_vec();
    const auto k = moduli.size();
    const auto n = phase.dimension();
    const u64 t = ct.plain_modulus;

    // CRT 重建系数 + scale-and-round: pt_j = round(t·x_j/Q) mod t。
    const UBInt Q = product_of(moduli);
    const UBInt half_Q = Q / UBInt(2);
    const UBInt t_big(t);
    std::vector<UBInt> crt_coeff(k); // M_i · (M_i^{-1} mod q_i)
    for (size_t i = 0; i < k; i++) {
        const u64 q_i = moduli[i];
        const UBInt M_i = Q / UBInt(q_i);
        const u64 m_i = to_u64(M_i % UBInt(q_i));
        const u64 g_i = inverse_mod_prime(m_i, q_i);
        crt_coeff[i] = M_i * UBInt(g_i);
    }

    BfvPt pt(RnsPolyParams{n, 1, std::vector<u64>{t}});
    pt.rep_form = PolyRepForm::coeff;
    auto out = pt[0].data();
    for (size_t j = 0; j < n; j++) {
        UBInt x;
        for (size_t i = 0; i < k; i++) {
            x += UBInt(phase[i][j]) * crt_coeff[i];
        }
        x %= Q;
        // round(t·x/Q) = (t·x + Q/2) / Q
        UBInt r = (t_big * x + half_Q) / Q;
        r %= t_big;
        out[j] = to_u64(r);
    }
    return pt;
}

} // namespace bfv
} // namespace hehub
