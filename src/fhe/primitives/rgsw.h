/**
 * @file rgsw.h
 * @brief RNS-gadget RGSW 密文与 RLWE⊡RGSW 外积 (服务 BGV/BFV/CKKS 的重线性化/密钥切换)。
 *
 * 这里的 RGSW 用 **RNS(CRT 基) 作为 gadget 分解基**: 一个 RgswCt 由若干 RLWE 样本构成,
 * 第 k 个样本在 c0 上叠加 pt·(gadget 分量 g_k), 其中 g_k 取自 decomp_basis(通常是
 * ∏_{i≠k} q_i · [ (∏_{i≠k} q_i)^{-1} mod q_k ] 这类 RNS/CRT gadget 向量)。外积时对输入
 * RLWE 做 RNS 分解逐分量与 RGSW 内积, 从而在不放大噪声的前提下实现"明文乘"。
 *
 * ⚠️ 与 tfhe/fhew 的完整 RGSW 区分: 那套用 **单模数 + power-of-B gadget**, 专供门自举
 * 盲旋转; 本套用 **多模数 RNS gadget**, 单模数下会退化不能当 gadget。见 tfhe/blind_rotate.h。
 */

#pragma once

#include "common/rns.h"
#include "rlwe.h"
#include <array>

namespace hehub {

/**
 * @brief RGSW 规格占位(预留: 维度/模数/gadget 基长等参数聚合)。当前 RGSW 的形状由
 *        入参 decomp_basis 与私钥的 RnsPolyParams 隐式决定, 故此结构暂为空。
 */
struct RgswSpec {};

/// @brief RGSW 密文: 一组 RLWE 样本(每个对应一个 gadget 分量行)。
using RgswCt = std::vector<RlweCt>;

/**
 * @brief RNS-gadget RGSW 加密。生成 decomp_basis.size() 个 RLWE 随机样本, 并在第 k 个
 *        样本的 c0 上叠加 pt_ntt·decomp_basis[k](逐 RNS 分量标量乘), 得到 gadget 阶梯。
 * @param pt_ntt 明文, **必须为 NTT(value) form**(否则抛异常)。
 * @param sk RLWE 私钥(NTT form)。
 * @param decomp_basis gadget 分解基, 每行是一组按 RNS 分量给出的标量(∈各模数)。
 * @return RgswCt 大小等于 decomp_basis.size() 的 RLWE 样本组。
 */
RgswCt rgsw_encrypt(const RlwePt &pt_ntt, const RlweSk &sk,
                    const std::vector<std::vector<u64>> &decomp_basis);

/**
 * @brief 同 rgsw_encrypt, 但把每个多项式额外乘以 Montgomery 常数 (2^64 mod q_i),
 *        使其可直接喂给 ext_prod_montgomery(外积内部用 128 位 Montgomery 归约)。
 * @param pt_ntt 明文(NTT form)。
 * @param sk RLWE 私钥。
 * @param decomp_basis gadget 分解基。
 * @return RgswCt Montgomery 域下的 RGSW 密文。
 */
RgswCt rgsw_encrypt_montgomery(const RlwePt &pt_ntt, const RlweSk &sk,
                               const std::vector<std::vector<u64>> &decomp_basis);

/**
 * @brief RLWE⊡RGSW 外积 (Montgomery 变体)。对输入 pt 做 RNS 分解构成分量矩阵, 与
 *        Montgomery 域 RGSW 内积并用 batched_montgomery_128 归约, 输出扩展模数下的
 *        RLWE 密文(相位 ≈ μ·phase(pt), μ 为 RGSW 承载的明文)。用于 relin/密钥切换。
 * @param pt 输入多项式(单/多模数, value form)。
 * @param rgsw Montgomery 域 RGSW 密文(须与 pt 的模数/维度相容, 且多一个扩展模数)。
 * @return RlweCt 扩展模数空间下的 RLWE 结果(value form)。
 */
RlweCt ext_prod_montgomery(const RlwePt &pt, const RgswCt &rgsw);

} // namespace hehub
