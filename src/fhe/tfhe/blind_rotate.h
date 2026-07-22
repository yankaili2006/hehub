/**
 * @file blind_rotate.h
 * @brief FHEW/TFHE 门自举地基: 完整形式 RGSW、RLWE⊡RGSW 外积、CMux、盲旋转。
 *
 * 与 primitives/rgsw.h 的区别: 那套是 **RNS-gadget、单多项式** 外积(服务 BFV/CKKS
 * relin 密钥切换), 在单模数场景下 RNS 分解退化不能当 gadget。此处是 **单模数 +
 * power-of-B gadget** 的完整 RGSW 与 RLWE⊡RGSW(二者皆密文)外积, 专供门自举。
 *
 * 相位约定沿用 rlwe.cpp: phase(c0,c1) = c0 + c1·s。据此
 *   RGSW(μ) = ( RLWE'(μ), RLWE'(μ·s) ),  外积 ct⊡RGSW(μ) = RLWE(μ·phase(ct))。
 * RGSW 用 RgswCt(=vector<RlweCt>) 承载 2ℓ 行: 前 ℓ 行 = RLWE'(μ) 各加密 μ·B^j 到 c0,
 * 后 ℓ 行 = RLWE'(μ·s) 各加密 μ·s·B^j。全部单分量、value(NTT) form。
 */
#pragma once

#include "common/rns.h"
#include "primitives/lwe.h"
#include "primitives/rgsw.h"
#include "primitives/rlwe.h"
#include <vector>

namespace hehub {
namespace tfhe {

/// @brief TFHE 参数: 环维度 N、单模数 q、gadget 基位宽与层数。
struct TfheParams {
    size_t N = 0;         // 环维度 (2 的幂)
    u64 q = 0;            // 单个 NTT 友好素数 (≡1 mod 2N)
    size_t base_bits = 0; // gadget 基 B = 2^base_bits
    size_t levels = 0;    // gadget 层数 ℓ, 满足 B^ℓ ≥ 2q
};

/// @brief 由 (N, q) 与固定基位宽推导完整 TfheParams (加/解密两侧必须一致)。
TfheParams make_tfhe_params(size_t N, u64 q, size_t base_bits = 8);

/// @brief 完整 RGSW 加密 μ(coeff form, 单分量, mod q): 返回 2ℓ 行 RlweCt。
///        前 ℓ 行加密 μ·B^j, 后 ℓ 行加密 μ·s·B^j。
RgswCt tfhe_rgsw_encrypt(const RnsPolynomial &mu_coeff, const RlweSk &sk,
                         const TfheParams &gp);

/// @brief RLWE⊡RGSW 外积: 输入密文 ct(value form) 与 RGSW(μ), 输出 RLWE(μ·phase(ct))。
RlweCt tfhe_external_product(const RlweCt &ct, const RgswCt &rgsw,
                            const TfheParams &gp);

/// @brief CMux(选择子 RGSW(bit), ct0, ct1) = bit ? ct1 : ct0。
RlweCt cmux(const RgswCt &bit_rgsw, const RlweCt &ct0, const RlweCt &ct1,
            const TfheParams &gp);

/// @brief 生成自举密钥: 二元 LWE 私钥各分量 s_i∈{0,1} 在 RLWE 私钥下的 RGSW(s_i)。
std::vector<RgswCt> gen_bootstrap_keys(const LweSk &lwe_sk,
                                       const RlweSk &rlwe_sk,
                                       const TfheParams &gp);

/// @brief LWE 模切换 q→new_mod (相位按比例缩放取整, 明文/密钥不变)。
LweCt lwe_mod_switch(const LweCt &ct, u64 new_mod);

/// @brief 盲旋转: 输入模 2N 的 LWE(a,b)(二元私钥), 测试多项式 testvec(coeff form),
///        自举密钥 bk; 输出累加器 RLWE, 其相位 = testvec·X^{-(b+Σa_i s_i)}。
RlweCt blind_rotate(const LweCt &ct_mod2N, const RnsPolynomial &testvec_coeff,
                    const std::vector<RgswCt> &bk, const TfheParams &gp);

// ---- 测试/上层便捷封装 (缩放明文, 容忍噪声) ----

/// @brief 缩放 RLWE 加密: phase = Δ·pt + e, Δ=round(q/t)。pt coeff form mod t。
RlweCt tfhe_rlwe_encrypt(const RnsPolynomial &pt, u64 t, const RlweSk &sk,
                         const TfheParams &gp);

/// @brief 缩放 RLWE 解密: round(t·phase/q) mod t (coeff form)。
RnsPolynomial tfhe_rlwe_decrypt(const RlweCt &ct, u64 t, const RlweSk &sk,
                                const TfheParams &gp);

} // namespace tfhe
} // namespace hehub
