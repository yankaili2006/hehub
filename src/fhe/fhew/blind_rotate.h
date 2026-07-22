/**
 * @file blind_rotate.h
 * @brief FHEW (AP 法) 盲旋转与函数自举。
 *
 * 与 tfhe/blind_rotate.h(GINX 法)的区别在**自举密钥结构**与**私钥支持范围**:
 *
 *   - GINX: 自举密钥 = 每个私钥分量 RGSW(s_i), 盲旋转逐维 CMux。**要求二元私钥**
 *           s_i∈{0,1}(单次 CMux 即得 X^{a_i·s_i})。密钥小(n 个 RGSW), 但私钥受限。
 *   - AP:   自举密钥 = 一张**按数字展开的表** ek[i][j][v] = RGSW(X^{-(v·B_r^j·s_i)}),
 *           盲旋转把 a_i 按刷新基 B_r 分解为数字 a_{i,j}, 逐数字做一次 RLWE⊡RGSW
 *           外积: ACC ← ACC ⊡ ek[i][j][a_{i,j}]。**支持任意(含三元 {-1,0,1})私钥**,
 *           因 s_i 在密钥生成期已烘焙进指数。代价是密钥表大(n·d_r·B_r 个 RGSW)、
 *           外积次数 n·d_r 略多。
 *
 * 两法共用 tfhe 的完整 RGSW(tfhe_rgsw_encrypt)与 RLWE⊡RGSW 外积(tfhe_external_product):
 * FHEW 的密钥项只是 RGSW 加密的单项式 X^e, 外积语义 ACC⊡RGSW(X^e)=RLWE(X^e·phase(ACC))
 * 恰好实现"给累加器的明文乘上单项式"。相位约定同 rlwe.cpp: phase(c0,c1)=c0+c1·s。
 *
 * 盲旋转最终 ACC 的相位 = testvec·X^{-(b+Σa_i s_i)}, 与 GINX 一致(仅内部实现不同),
 * 故 tfhe 侧的 sample_extract / functional_bootstrap 语义可直接复用。
 */
#pragma once

#include "common/rns.h"
#include "primitives/lwe.h"
#include "primitives/rgsw.h"
#include "primitives/rlwe.h"
#include "tfhe/blind_rotate.h"
#include <vector>

namespace hehub {
namespace fhew {

/// @brief FHEW 参数: 复用 tfhe 的 RGSW gadget 参数(N,q,base_bits,levels), 另加刷新基。
struct FhewParams {
    tfhe::TfheParams gp;    // RGSW/外积用的 gadget 参数 (N, q, base_bits, levels)
    size_t br_bits = 0;     // 刷新基 B_r = 2^br_bits (a_i 的数字分解基)
    size_t br_digits = 0;   // d_r = ceil(log_{B_r}(2N)), 满足 B_r^{d_r} ≥ 2N
    size_t base_r() const { return (size_t)1 << br_bits; }
};

/// @brief 由 (N, q) 推导完整 FhewParams。base_bits 为 RGSW gadget 基位宽,
///        br_bits 为刷新基位宽(越大→密钥表越大但外积次数越少→噪声越低)。
FhewParams make_fhew_params(size_t N, u64 q, size_t base_bits = 8,
                            size_t br_bits = 6);

/// @brief FHEW 自举密钥: ek[i][j][v] = RGSW(X^{-(v·B_r^j·s_i) mod 2N})。
///        i∈[n] 私钥分量, j∈[d_r] 数字位, v∈[0,B_r) 数字值。
///        v=0 项恒为 RGSW(1)(外积单位元), 盲旋转时直接跳过, 不实际存储。
struct FhewBootstrapKey {
    // [n][d_r][B_r]; ek[i][j][0] 留空(未使用), 其余为有效 RGSW。
    std::vector<std::vector<std::vector<RgswCt>>> ek;
    FhewParams params;
    size_t n = 0; // LWE 维度
};

/// @brief 生成 FHEW 自举密钥(AP 表)。lwe_sk 可为三元 {-1,0,1}(AP 法的核心优势)。
FhewBootstrapKey gen_fhew_bootstrap_key(const LweSk &lwe_sk,
                                        const RlweSk &rlwe_sk,
                                        const FhewParams &fp);

/// @brief FHEW 盲旋转(AP 法)。输入模 2N 的 LWE(a,b), 测试多项式 testvec(coeff form),
///        自举密钥 bk; 输出累加器 RLWE, 其相位 = testvec·X^{-(b+Σa_i s_i)}。
RlweCt fhew_blind_rotate(const LweCt &ct_mod2N,
                         const RnsPolynomial &testvec_coeff,
                         const FhewBootstrapKey &bk);

/// @brief FHEW 函数自举: 模切换 q→2N → 盲旋转 → 提取常数项 LWE。
///        与 tfhe::functional_bootstrap 等价, 但走 AP 法(支持三元私钥)。
LweCt fhew_functional_bootstrap(const LweCt &ct, const RnsPolynomial &lut_poly,
                                const FhewBootstrapKey &bk);

} // namespace fhew
} // namespace hehub
