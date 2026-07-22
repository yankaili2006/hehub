#include "catch2/catch.hpp"
#include "common/rns.h"
#include "common/sampling.h"
#include "primitives/lwe.h"
#include "primitives/rlwe.h"
#include "tfhe/blind_rotate.h"
#include "tfhe/func_boot.h"
#include <vector>

using namespace hehub;
using namespace hehub::tfhe;

namespace {
// 单模数 q≈2^30, ≡1 mod 2N (N=1024 → 2N=2048)。
constexpr u64 kQ = 1073643521;
constexpr size_t kN = 1024;

// 构造 coeff form 单分量多项式(值取自 vals, 其余为 0)。
RnsPolynomial make_coeff_poly(const std::vector<u64> &vals) {
    auto p = get_zero_poly(RnsPolyParams{kN, 1, {kQ}}, PolyRepForm::coeff);
    for (size_t j = 0; j < vals.size() && j < kN; j++) {
        p[0][j] = vals[j] % kQ;
    }
    return p;
}
} // namespace

TEST_CASE("tfhe RGSW external product + CMux") {
    auto gp = make_tfhe_params(kN, kQ);
    RlweSk sk(RnsPolyParams{kN, 1, {kQ}});
    const u64 t = 16;

    // 两个已知小明文多项式。
    std::vector<u64> v0(kN), v1(kN);
    u64 s = 5;
    for (size_t j = 0; j < kN; j++) {
        v0[j] = (s = s * 1103515245 + 12345) % t;
        v1[j] = (s = s * 1103515245 + 12345) % t;
    }
    auto p0 = make_coeff_poly(v0);
    auto p1 = make_coeff_poly(v1);

    auto d0 = tfhe_rlwe_encrypt(p0, t, sk, gp);
    auto d1 = tfhe_rlwe_encrypt(p1, t, sk, gp);

    // RGSW(0) 与 RGSW(1) (常数明文)。
    auto rgsw0 = tfhe_rgsw_encrypt(make_coeff_poly({0}), sk, gp);
    auto rgsw1 = tfhe_rgsw_encrypt(make_coeff_poly({1}), sk, gp);

    SECTION("外积 μ=1 复原相位 (RLWE⊡RGSW(1) ≈ 原密文)") {
        auto prod = tfhe_external_product(d0, rgsw1, gp);
        auto dec = tfhe_rlwe_decrypt(prod, t, sk, gp);
        for (size_t j = 0; j < kN; j++) {
            REQUIRE(dec[0][j] == v0[j]);
        }
    }
    SECTION("CMux(bit=0) 选择 ct0") {
        auto out = cmux(rgsw0, d0, d1, gp);
        auto dec = tfhe_rlwe_decrypt(out, t, sk, gp);
        for (size_t j = 0; j < kN; j++) {
            REQUIRE(dec[0][j] == v0[j]);
        }
    }
    SECTION("CMux(bit=1) 选择 ct1") {
        auto out = cmux(rgsw1, d0, d1, gp);
        auto dec = tfhe_rlwe_decrypt(out, t, sk, gp);
        for (size_t j = 0; j < kN; j++) {
            REQUIRE(dec[0][j] == v1[j]);
        }
    }
}

TEST_CASE("tfhe blind rotation (FBS) with staircase LUT") {
    // 负循环 naive FBS: 要求相位 φ=b+Σa_i s_i ∈ (0,N) 且远离 0/N 边界(否则噪声
    // 使 φ 绕过 2N 触发 X^N=-1 的取负, 命中 LUT 反射区)。取 t=8, m∈{1,2,3} →
    // 模切换后 φ≈2N/t·m = 256·m ∈ {256,512,768}, 余量 256 ≫ 噪声。
    // m=0 与近 t/2 边界属负循环固有限制, 由 FFBS(redundant-MSB) 后续消除。
    auto gp = make_tfhe_params(kN, kQ);
    RlweSk rlwe_sk(RnsPolyParams{kN, 1, {kQ}});

    const size_t n = 16;
    const u64 t = 8;
    LweSk lwe_sk = sample_binary_lwe_sk(LweParams{n, kQ});
    auto bk = gen_bootstrap_keys(lwe_sk, rlwe_sk, gp);

    // 阶梯 LUT(恒等函数): testvec[j] = Δ_out·round(j·t/2N), Δ_out=round(q/t)。
    const u64 delta_out = (u64)(((unsigned __int128)kQ + t / 2) / t);
    std::vector<u64> lut(kN);
    for (size_t j = 0; j < kN; j++) {
        u64 slot = (u64)(((unsigned __int128)j * t + kN) / (2 * kN)); // round(j·t/2N)
        lut[j] = (u64)(((unsigned __int128)delta_out * slot) % kQ);
    }
    auto lut_poly = make_coeff_poly(lut);

    // 提取用的 LWE 私钥 = RLWE 私钥系数。
    LweSk extract_sk = lwe_sk_from_rlwe(rlwe_sk);

    for (u64 m : {(u64)1, (u64)2, (u64)3}) {
        auto ct = lwe_encrypt(m, t, lwe_sk);
        auto boot = functional_bootstrap(ct, lut_poly, bk);
        u64 out = lwe_decrypt(boot, t, extract_sk);
        REQUIRE(out == m);
    }
}
