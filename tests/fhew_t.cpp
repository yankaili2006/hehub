#include "catch2/catch.hpp"
#include "common/rns.h"
#include "common/sampling.h"
#include "fhew/blind_rotate.h"
#include "primitives/lwe.h"
#include "primitives/rlwe.h"
#include "tfhe/blind_rotate.h"
#include <vector>

using namespace hehub;

namespace {
// 单模数 q≈2^30, ≡1 mod 2N (N=1024 → 2N=2048)。
constexpr u64 kQ = 1073643521;
constexpr size_t kN = 1024;

RnsPolynomial make_coeff_poly(const std::vector<u64> &vals) {
    auto p = get_zero_poly(RnsPolyParams{kN, 1, {kQ}}, PolyRepForm::coeff);
    for (size_t j = 0; j < vals.size() && j < kN; j++) {
        p[0][j] = vals[j] % kQ;
    }
    return p;
}
} // namespace

TEST_CASE("fhew blind rotation (AP) with ternary LWE key") {
    // AP 法核心优势: 支持**三元**私钥 s_i∈{-1,0,1}(GINX 的 CMux 只能二元)。
    // 负循环 naive FBS: 取 t=8, m∈{1,2,3} → 模切换后相位 φ≈256·m ∈{256,512,768},
    // 距 0/N 边界余量 256 ≫ 噪声。m=0 与近 t/2 边界属负循环固有限制。
    auto fp = fhew::make_fhew_params(kN, kQ);
    RlweSk rlwe_sk(RnsPolyParams{kN, 1, {kQ}});

    const size_t n = 8;
    const u64 t = 8;
    // 三元私钥(LweSk 默认构造即三元采样)。
    LweSk lwe_sk(LweParams{n, kQ});
    auto bk = fhew::gen_fhew_bootstrap_key(lwe_sk, rlwe_sk, fp);

    // 阶梯 LUT(恒等函数): testvec[j] = Δ_out·round(j·t/2N)。
    const u64 delta_out = (u64)(((unsigned __int128)kQ + t / 2) / t);
    std::vector<u64> lut(kN);
    for (size_t j = 0; j < kN; j++) {
        u64 slot = (u64)(((unsigned __int128)j * t + kN) / (2 * kN));
        lut[j] = (u64)(((unsigned __int128)delta_out * slot) % kQ);
    }
    auto lut_poly = make_coeff_poly(lut);

    // 提取用 LWE 私钥 = RLWE 私钥系数。
    LweSk extract_sk = lwe_sk_from_rlwe(rlwe_sk);

    for (u64 m : {(u64)1, (u64)2, (u64)3}) {
        auto ct = lwe_encrypt(m, t, lwe_sk);
        auto boot = fhew::fhew_functional_bootstrap(ct, lut_poly, bk);
        u64 out = lwe_decrypt(boot, t, extract_sk);
        REQUIRE(out == m);
    }
}

TEST_CASE("fhew params: digit count covers 2N") {
    auto fp = fhew::make_fhew_params(kN, kQ, 8, 6);
    // B_r^{d_r} ≥ 2N: 保证 a_i∈[0,2N) 的数字分解精确无溢出。
    unsigned long long cap = 1;
    for (size_t j = 0; j < fp.br_digits; j++) {
        cap *= fp.base_r();
    }
    REQUIRE(cap >= 2 * kN);
    // 且恰好足够(去掉一位不够), 即数字位数最小。
    unsigned long long cap_minus1 = 1;
    for (size_t j = 0; j + 1 < fp.br_digits; j++) {
        cap_minus1 *= fp.base_r();
    }
    REQUIRE(cap_minus1 < 2 * kN);
}
