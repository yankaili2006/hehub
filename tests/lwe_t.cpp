#include "catch2/catch.hpp"
#include "primitives/lwe.h"

using namespace hehub;

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
