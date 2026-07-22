#include "bfv/bfv.h"
#include "catch2/catch.hpp"
#include "common/sampling.h"
#include <vector>

using namespace hehub;

TEST_CASE("bfv encoding") {
    u64 t = 65537;
    size_t n = 1024;
    std::vector<u64> data(n);
    u64 seed = 1;
    for (auto &d : data) {
        d = ((seed++) * 888 + 123) % t;
    }
    auto pt = bfv::simd_encode(data, t);
    auto decoded = bfv::simd_decode(pt);
    REQUIRE(decoded == data);
}

TEST_CASE("bfv encryption") {
    std::vector<u64> ct_moduli{131530753, 130809857};
    size_t dimension = 128;
    RnsPolyParams ct_params{dimension, ct_moduli.size(), ct_moduli};
    RlweSk sk(ct_params);

    SECTION("encrypt/decrypt round-trip") {
        u64 t = 65537;
        RnsPolyParams pt_params{dimension, 1, std::vector{t}};
        BfvPt pt = get_rand_uniform_poly(pt_params);

        auto ct = bfv::encrypt(pt, sk);
        auto pt_rec = bfv::decrypt(ct, sk);
        REQUIRE(pt == pt_rec);
    }

    SECTION("homomorphic add") {
        u64 t = 65537;
        RnsPolyParams pt_params{dimension, 1, std::vector{t}};
        BfvPt pt1 = get_rand_uniform_poly(pt_params);
        BfvPt pt2 = get_rand_uniform_poly(pt_params);

        auto ct1 = bfv::encrypt(pt1, sk);
        auto ct2 = bfv::encrypt(pt2, sk);
        auto ct_sum = bfv::add(ct1, ct2);
        auto pt_sum = bfv::decrypt(ct_sum, sk);

        // expected = (pt1 + pt2) mod t, 逐系数
        for (size_t j = 0; j < dimension; j++) {
            u64 expect = (pt1[0][j] + pt2[0][j]) % t;
            REQUIRE(pt_sum[0][j] == expect);
        }
    }

    SECTION("add_plain") {
        u64 t = 65537;
        RnsPolyParams pt_params{dimension, 1, std::vector{t}};
        BfvPt pt1 = get_rand_uniform_poly(pt_params);
        BfvPt pt2 = get_rand_uniform_poly(pt_params);

        auto ct1 = bfv::encrypt(pt1, sk);
        auto ct_sum = bfv::add_plain(ct1, pt2);
        auto pt_sum = bfv::decrypt(ct_sum, sk);

        for (size_t j = 0; j < dimension; j++) {
            u64 expect = (pt1[0][j] + pt2[0][j]) % t;
            REQUIRE(pt_sum[0][j] == expect);
        }
    }

    SECTION("SIMD encrypt/decrypt + mult_plain (slot-wise)") {
        u64 t = 65537;
        size_t slots = dimension;
        std::vector<u64> a(slots), b(slots);
        u64 s = 7;
        for (size_t i = 0; i < slots; i++) {
            a[i] = (s = s * 1103515245 + 12345) % t;
            b[i] = (s = s * 1103515245 + 12345) % t;
        }
        auto pt_a = bfv::simd_encode(a, t, slots);
        auto pt_b = bfv::simd_encode(b, t, slots);

        auto ct = bfv::encrypt(pt_a, sk);
        // 密文-明文乘: 槽位乘 a[i]*b[i]
        auto ct_prod = bfv::mult_plain(ct, pt_b);
        auto dec = bfv::simd_decode(bfv::decrypt(ct_prod, sk), slots);

        for (size_t i = 0; i < slots; i++) {
            u64 expect = (u64)(((unsigned __int128)a[i] * b[i]) % t);
            REQUIRE(dec[i] == expect);
        }
    }

    SECTION("coprime condition") {
        u64 t = 131530753; // 与某个密文模相同 → 非互素, 应抛异常
        RnsPolyParams pt_params{dimension, 1, std::vector{t}};
        BfvPt pt(pt_params);
        REQUIRE_THROWS(bfv::encrypt(pt, sk));
    }
}
