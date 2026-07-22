#include "lwe.h"
#include <cmath>
#include <random>
#include <stdexcept>

namespace hehub {

namespace {
// 模块级 RNG (LWE 无需 NTT, 自包含随机源)。
std::mt19937_64 &lwe_rng() {
    static std::mt19937_64 eng(0xC0FFEEULL);
    return eng;
}

int sample_ternary() {
    static std::uniform_int_distribution<int> d(-1, 1);
    return d(lwe_rng());
}

// 高斯噪声(取整), std≈3.2。
i64 sample_gaussian(double std_dev = 3.2) {
    static std::normal_distribution<double> d(0.0, 1.0);
    return (i64)std::llround(d(lwe_rng()) * std_dev);
}

u64 rand_mod(u64 q) {
    std::uniform_int_distribution<u64> d(0, q - 1);
    return d(lwe_rng());
}

// (b + <a,s>) mod q, s_i∈{-1,0,1}。
u64 phase(const std::vector<u64> &a, u64 b, const std::vector<int> &s, u64 q) {
    unsigned __int128 acc = b % q;
    for (size_t i = 0; i < a.size(); i++) {
        if (s[i] == 1) {
            acc += a[i];
        } else if (s[i] == -1) {
            acc += (q - a[i] % q) % q;
        }
        acc %= q;
    }
    return (u64)acc;
}
} // namespace

LweSk::LweSk(const LweParams &params)
    : s(params.dimension), modulus(params.modulus) {
    for (auto &si : s) {
        si = sample_ternary();
    }
}

LweCt lwe_encrypt(u64 message, u64 plain_modulus, const LweSk &sk) {
    const u64 q = sk.modulus;
    const size_t n = sk.s.size();
    const u64 t = plain_modulus;
    // Δ = round(q/t)
    const u64 Delta = (u64)std::llround((double)q / (double)t);

    LweCt ct;
    ct.modulus = q;
    ct.a.resize(n);
    for (size_t i = 0; i < n; i++) {
        ct.a[i] = rand_mod(q);
    }
    // b = Δ·m + e - <a,s> (mod q)
    i64 e = sample_gaussian();
    unsigned __int128 body = (unsigned __int128)(Delta % q) * (message % t) % q;
    // + e (处理负)
    i64 em = e % (i64)q;
    if (em < 0)
        em += (i64)q;
    body = (body + (u64)em) % q;
    // - <a,s>
    u64 as = phase(ct.a, 0, sk.s, q);
    ct.b = (u64)((body + (q - as)) % q);
    return ct;
}

u64 lwe_decrypt(const LweCt &ct, u64 plain_modulus, const LweSk &sk) {
    const u64 q = ct.modulus;
    const u64 t = plain_modulus;
    u64 ph = phase(ct.a, ct.b, sk.s, q); // = Δ·m + e
    // round(t·ph/q) mod t
    unsigned __int128 num = (unsigned __int128)t * ph + q / 2;
    u64 m = (u64)(num / q) % t;
    return m;
}

LweCt lwe_add(const LweCt &x, const LweCt &y) {
    if (x.modulus != y.modulus || x.dimension() != y.dimension()) {
        throw std::invalid_argument("LWE add: params mismatch.");
    }
    const u64 q = x.modulus;
    LweCt out;
    out.modulus = q;
    out.a.resize(x.dimension());
    for (size_t i = 0; i < x.dimension(); i++) {
        out.a[i] = (u64)(((unsigned __int128)x.a[i] + y.a[i]) % q);
    }
    out.b = (u64)(((unsigned __int128)x.b + y.b) % q);
    return out;
}

LweCt lwe_sub(const LweCt &x, const LweCt &y) {
    if (x.modulus != y.modulus || x.dimension() != y.dimension()) {
        throw std::invalid_argument("LWE sub: params mismatch.");
    }
    const u64 q = x.modulus;
    LweCt out;
    out.modulus = q;
    out.a.resize(x.dimension());
    for (size_t i = 0; i < x.dimension(); i++) {
        out.a[i] = (u64)(((unsigned __int128)x.a[i] + (q - y.a[i] % q)) % q);
    }
    out.b = (u64)(((unsigned __int128)x.b + (q - y.b % q)) % q);
    return out;
}

LweCt lwe_negate(const LweCt &x) {
    const u64 q = x.modulus;
    LweCt out;
    out.modulus = q;
    out.a.resize(x.dimension());
    for (size_t i = 0; i < x.dimension(); i++) {
        out.a[i] = (x.a[i] == 0) ? 0 : (q - x.a[i]);
    }
    out.b = (x.b == 0) ? 0 : (q - x.b);
    return out;
}

} // namespace hehub
