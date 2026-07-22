#include "bfv.h"
#include "common/bigint.h"
#include "common/ntt.h"
#include "common/rns_transform.h"

namespace hehub {
namespace bfv {

namespace {
// 把明文按 Δ=⌊Q/t⌋ 缩放到密文模 RNS 下并转 NTT value form(供与 c0 加/减)。
RnsPolynomial delta_scaled_pt(const BfvPt &pt, const std::vector<u64> &ct_moduli) {
    auto pt_under_ct = rns_base_transform(pt, ct_moduli);
    UBInt Q(1);
    for (auto m : ct_moduli) {
        Q *= UBInt(m);
    }
    const UBInt Delta = Q / UBInt(pt.modulus_at(0));
    for (size_t i = 0; i < ct_moduli.size(); i++) {
        const u64 q_i = ct_moduli[i];
        const u64 delta_i = to_u64(Delta % UBInt(q_i));
        auto data = pt_under_ct[i].data();
        for (size_t j = 0; j < pt_under_ct.dimension(); j++) {
            data[j] = (u64)(((unsigned __int128)data[j] * delta_i) % q_i);
        }
    }
    ntt_negacyclic_inplace_lazy(pt_under_ct);
    return pt_under_ct;
}
} // namespace

BfvCt add(const BfvCt &ct1, const BfvCt &ct2) {
    BfvCt out = RlweCt{ct1[0] + ct2[0], ct1[1] + ct2[1]};
    out.plain_modulus = ct1.plain_modulus;
    return out;
}

BfvCt sub(const BfvCt &ct1, const BfvCt &ct2) {
    BfvCt out = RlweCt{ct1[0] - ct2[0], ct1[1] - ct2[1]};
    out.plain_modulus = ct1.plain_modulus;
    return out;
}

BfvCt add_plain(const BfvCt &ct, const BfvPt &pt) {
    auto scaled = delta_scaled_pt(pt, ct[0].modulus_vec());
    BfvCt out = ct;
    out[0] += scaled;
    return out;
}

BfvCt sub_plain(const BfvCt &ct, const BfvPt &pt) {
    auto scaled = delta_scaled_pt(pt, ct[0].modulus_vec());
    BfvCt out = ct;
    out[0] -= scaled;
    return out;
}

} // namespace bfv
} // namespace hehub
