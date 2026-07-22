/**
 * @file bfv_internal.h
 * @brief BFV 内部共享工具(不对外导出)。
 */
#pragma once
#include "bfv.h"
#include "common/bigint.h"
#include "common/rns.h"

namespace hehub {
namespace bfv {

/// @brief ∏ moduli 的大整数。
UBInt product_of(const std::vector<u64> &moduli);

/// @brief 从 128 位无符号值构造 UBInt。
UBInt ubint_from_u128(unsigned __int128 v);

/// @brief 相位多项式(coeff form, RNS)做 BFV scale-and-round:
///        逐系数 CRT 重建 x∈[0,Q), 取 round(t·x/Q) mod t, 得明文(coeff form, mod t)。
BfvPt scale_and_round(const RnsPolynomial &phase, u64 t);

} // namespace bfv
} // namespace hehub
