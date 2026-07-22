/**
 * @file bfv.h
 * @brief BFV scheme (Brakerski/Fan-Vercauteren) — 整数 SIMD 全同态(leveled)。
 *
 * 与 BGV 同为 RLWE 上的整数方案, 复用 common(RNS/NTT) 与 primitives(RLWE)。
 * 区别: 明文按 Δ=⌊Q/t⌋ 缩放置于高位(MSB), 噪声在低位; 解密 scale-and-round。
 * 首版覆盖: 编码/解码(复用 BGV 打包) + 加密/解密 + 同态加/明文加/减。
 * 乘法(需 tensor + t/Q 缩放 + relin)留待后续增量。
 */

#pragma once
#include "primitives/keys.h"
#include "primitives/rlwe.h"

namespace hehub {
namespace bfv {

/// @brief BFV 明文 = RLWE 明文(系数在 Z_t 上, coeff form)。
using BfvPt = RlwePt;

/// @brief BFV 密文。plain_modulus = t。
struct BfvCt : public RlweCt {
    using RlweCt::RlweCt;
    BfvCt(RlweCt &&other) : RlweCt(std::move(other)) {}
    u64 plain_modulus = 1;
};

/// @brief SIMD 打包编码(与 BGV 一致的 CRT 批处理)。
BfvPt simd_encode(const std::vector<u64> &data, const u64 modulus,
                  size_t slot_count = 0);

/// @brief SIMD 解码。
std::vector<u64> simd_decode(const BfvPt &pt, size_t data_size = 0);

/// @brief 加密: c0 = -a·s + e + Δ·pt, c1 = a (Δ=⌊Q/t⌋)。
BfvCt encrypt(const BfvPt &pt, const RlweSk &rlwe_sk,
              std::vector<u64> ct_moduli = std::vector<u64>{});

/// @brief 解密: phase=(c0+c1·s) mod Q = Δ·pt+e, 再 round(t·phase/Q) mod t。
BfvPt decrypt(const BfvCt &ct, const RlweSk &rlwe_sk);

/// @brief 同态加(密+密)。
BfvCt add(const BfvCt &ct1, const BfvCt &ct2);

/// @brief 明文加(密+明)。
BfvCt add_plain(const BfvCt &ct, const BfvPt &pt);

/// @brief 同态减(密-密)。
BfvCt sub(const BfvCt &ct1, const BfvCt &ct2);

/// @brief 明文减(密-明)。
BfvCt sub_plain(const BfvCt &ct, const BfvPt &pt);

} // namespace bfv

// export types
using BfvPt = bfv::BfvPt;
using BfvCt = bfv::BfvCt;

} // namespace hehub
