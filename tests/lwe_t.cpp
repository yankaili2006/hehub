#include "bfv/bfv.h"
#include "catch2/catch.hpp"
#include "common/sampling.h"
#include "primitives/lwe.h"

using namespace hehub;

TEST_CASE("lwe sample extraction from RLWE") {
    // 单模数 RLWE(用 BFV 加密), 提取若干系数的 LWE 并用导出的 LWE 私钥解密。
    std::vector<u64> ct_moduli{1073643521}; // 单个 NTT 友好素数 (≡1 mod 2n)
    size_t dimension = 256;
    RnsPolyParams ct_params{dimension, ct_moduli.size(), ct_moduli};
    RlweSk sk(ct_params);

    u64 t = 17;
    RnsPolyParams pt_params{dimension, 1, std::vector{t}};
    BfvPt pt = get_rand_uniform_poly(pt_params); // coeff form, mod t

    auto ct = bfv::encrypt(pt, sk, ct_moduli);
    auto lwe_sk = lwe_sk_from_rlwe(sk);

    for (size_t k : {(size_t)0, (size_t)1, (size_t)5, dimension - 1}) {
        auto lwe = sample_extract(ct, k);
        u64 m = lwe_decrypt(lwe, t, lwe_sk);
        REQUIRE(m == pt[0][k]);
    }
}

TEST_CASE("lwe key switch") {
    // sk_from(维度 N=256) → sk_to(维度 n=64), 同模数, 明文保持不变。
    u64 q = (u64)1 << 32;
    u64 t = 16;
    LweSk sk_from(LweParams{256, q});
    LweSk sk_to(LweParams{64, q});
    auto ksk = gen_lwe_ksk(sk_from, sk_to, 8);

    for (u64 m = 0; m < t; m++) {
        auto ct = lwe_encrypt(m, t, sk_from);
        auto ct2 = lwe_key_switch(ct, ksk);
        REQUIRE(ct2.a.size() == 64);
        REQUIRE(lwe_decrypt(ct2, t, sk_to) == m);
    }
}

TEST_CASE("lwe basics") {
    LweParams params{512, (u64)1 << 32}; // n=512, q=2^32
    u64 t = 16;
    LweSk sk(params);

    SECTION("encrypt/decrypt round-trip") {
        for (u64 m = 0; m < t; m++) {
            auto ct = lwe_encrypt(m, t, sk);
            REQUIRE(lwe_decrypt(ct, t, sk) == m);
        }
    }
    SECTION("homomorphic add") {
        for (int iter = 0; iter < 20; iter++) {
            u64 m1 = iter % t, m2 = (iter * 3 + 1) % t;
            auto c1 = lwe_encrypt(m1, t, sk);
            auto c2 = lwe_encrypt(m2, t, sk);
            auto cs = lwe_add(c1, c2);
            REQUIRE(lwe_decrypt(cs, t, sk) == (m1 + m2) % t);
        }
    }
    SECTION("sub and negate") {
        u64 m1 = 9, m2 = 5;
        auto c1 = lwe_encrypt(m1, t, sk);
        auto c2 = lwe_encrypt(m2, t, sk);
        REQUIRE(lwe_decrypt(lwe_sub(c1, c2), t, sk) == (m1 + t - m2) % t);
        REQUIRE(lwe_decrypt(lwe_negate(c1), t, sk) == (t - m1) % t);
    }
}
