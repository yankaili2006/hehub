/**
 * @file lwe.h
 * @brief LWE 密文与基础运算 (FHEW/TFHE 门自举的地基)。
 *
 * LWE 密文 ct=(a, b), a∈Z_q^n, b∈Z_q, 满足 b + <a,s> = Δ·m + e (mod q),
 * Δ=round(q/t)。加/减为线性(不刷新噪声, 需自举); 自举/门运算见 tfhe。
 */
#pragma once

#include "common/type_defs.h"
#include "primitives/rlwe.h"
#include <cstddef>
#include <vector>

namespace hehub {

/// @brief LWE 参数: 维度 n 与模数 q。
struct LweParams {
    size_t dimension = 0; // n
    u64 modulus = 0;      // q
};

/// @brief LWE 密文 (a, b)。
struct LweCt {
    std::vector<u64> a; // 长度 n, 各分量 ∈ [0,q)
    u64 b = 0;
    u64 modulus = 0;
    size_t dimension() const { return a.size(); }
};

/// @brief LWE 三元私钥 s∈{-1,0,1}^n。
struct LweSk {
    std::vector<int> s; // ternary, 长度 n
    u64 modulus = 0;
    LweSk() = default;
    explicit LweSk(const LweParams &params);
};

/// @brief 加密 m∈Z_t: b + <a,s> = Δ·m + e (mod q), Δ=round(q/t)。
LweCt lwe_encrypt(u64 message, u64 plain_modulus, const LweSk &sk);

/// @brief 解密: round(t·(b+<a,s>)/q) mod t。
u64 lwe_decrypt(const LweCt &ct, u64 plain_modulus, const LweSk &sk);

/// @brief 同态加/减/取负 (线性, 逐分量 mod q)。
LweCt lwe_add(const LweCt &x, const LweCt &y);
LweCt lwe_sub(const LweCt &x, const LweCt &y);
LweCt lwe_negate(const LweCt &x);

/// @brief 由 RLWE 私钥(单模数)导出对应的 LWE 私钥(系数即 s_i, 供 sample extraction 后解密)。
LweSk lwe_sk_from_rlwe(const RlweSk &rlwe_sk);

/// @brief 从单模数 RLWE 密文提取第 coeff_index 个系数对应的 LWE 密文。
///        RLWE 相位 = c0 + c1·s (coeff form), 提取后 b+<a,s'> = 该系数的相位(=Δ·m_k+e),
///        s' = RLWE 私钥系数(见 lwe_sk_from_rlwe)。是门自举输出/中间步骤的关键。
LweCt sample_extract(const RlweCt &ct, size_t coeff_index = 0);

} // namespace hehub
