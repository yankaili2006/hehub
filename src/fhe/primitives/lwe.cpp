#include "lwe.h"
#include "common/mod_arith.h"
#include "common/ntt.h"
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

LweSk sample_binary_lwe_sk(const LweParams &params) {
    LweSk sk;
    sk.modulus = params.modulus;
    sk.s.resize(params.dimension);
    std::uniform_int_distribution<int> coin(0, 1);
    for (auto &si : sk.s) {
        si = coin(lwe_rng());
    }
    return sk;
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

LweSk lwe_sk_from_rlwe(const RlweSk &rlwe_sk) {
    if (rlwe_sk.component_count() != 1) {
        throw std::invalid_argument(
            "sample extraction 需单模数 RLWE 私钥(component_count==1)。");
    }
    // RlweSk 存为 NTT value form, 转 coeff 得系数, 再中心化到 {-1,0,1}。
    RnsPolynomial s = rlwe_sk;
    intt_negacyclic_inplace_lazy(s);
    reduce_strict(s);
    const u64 q = s.modulus_at(0);
    const size_t n = s.dimension();
    LweSk out;
    out.modulus = q;
    out.s.resize(n);
    for (size_t i = 0; i < n; i++) {
        u64 c = s[0][i];
        out.s[i] = (c == 0) ? 0 : (c == 1 ? 1 : -1); // 三元: 0, 1, 或 q-1(=-1)
    }
    return out;
}

namespace {
// 无缩放加密: b + <a,s> = value + e (mod q)。供 key-switch 钥。
LweCt lwe_encrypt_raw(u64 value, const LweSk &sk) {
    const u64 q = sk.modulus;
    const size_t n = sk.s.size();
    LweCt ct;
    ct.modulus = q;
    ct.a.resize(n);
    for (size_t i = 0; i < n; i++) {
        ct.a[i] = rand_mod(q);
    }
    i64 e = sample_gaussian();
    i64 em = e % (i64)q;
    if (em < 0)
        em += (i64)q;
    u64 body = (u64)(((unsigned __int128)(value % q) + (u64)em) % q);
    u64 as = phase(ct.a, 0, sk.s, q);
    ct.b = (u64)(((unsigned __int128)body + (q - as)) % q);
    return ct;
}

// out += scalar · ct (逐分量 mod q)。
void lwe_axpy(LweCt &out, u64 scalar, const LweCt &ct) {
    const u64 q = out.modulus;
    for (size_t i = 0; i < out.a.size(); i++) {
        out.a[i] = (u64)(((unsigned __int128)out.a[i] +
                          (unsigned __int128)scalar * ct.a[i]) %
                         q);
    }
    out.b = (u64)(((unsigned __int128)out.b + (unsigned __int128)scalar * ct.b) %
                  q);
}
} // namespace

LweKsk gen_lwe_ksk(const LweSk &sk_from, const LweSk &sk_to, size_t base_bits) {
    if (sk_from.modulus != sk_to.modulus) {
        throw std::invalid_argument("KSK: 模数需一致。");
    }
    const u64 q = sk_from.modulus;
    const size_t N = sk_from.s.size();
    const u64 base = (u64)1 << base_bits;
    size_t digits = 0;
    for (u64 tmp = q; tmp > 0; tmp >>= base_bits) {
        digits++;
    }

    LweKsk ksk;
    ksk.base_bits = base_bits;
    ksk.modulus = q;
    ksk.data.resize(N);
    for (size_t i = 0; i < N; i++) {
        ksk.data[i].resize(digits);
        u64 bj = 1; // base^j
        for (size_t j = 0; j < digits; j++) {
            // value = s_from_i · base^j (mod q), s_from_i ∈ {-1,0,1}
            u64 val;
            u64 sb = (u64)(((unsigned __int128)bj) % q);
            if (sk_from.s[i] == 1) {
                val = sb;
            } else if (sk_from.s[i] == -1) {
                val = (sb == 0) ? 0 : (q - sb);
            } else {
                val = 0;
            }
            ksk.data[i][j] = lwe_encrypt_raw(val, sk_to);
            bj = (u64)(((unsigned __int128)bj << base_bits) % q);
        }
    }
    return ksk;
}

LweCt lwe_key_switch(const LweCt &ct, const LweKsk &ksk) {
    const u64 q = ct.modulus;
    const size_t N = ct.a.size();
    const size_t n = ksk.data.empty() ? 0 : ksk.data[0][0].a.size();
    const u64 mask = ((u64)1 << ksk.base_bits) - 1;

    LweCt out;
    out.modulus = q;
    out.a.assign(n, 0);
    out.b = ct.b; // (0, b)
    for (size_t i = 0; i < N; i++) {
        u64 ai = ct.a[i];
        for (size_t j = 0; j < ksk.data[i].size(); j++) {
            u64 digit = (ai >> (j * ksk.base_bits)) & mask;
            if (digit != 0) {
                lwe_axpy(out, digit, ksk.data[i][j]);
            }
        }
    }
    return out;
}

LweCt sample_extract(const RlweCt &ct, size_t coeff_index) {
    if (ct[0].component_count() != 1 || ct[1].component_count() != 1) {
        throw std::invalid_argument("sample extraction 需单模数 RLWE 密文。");
    }
    // 转 coeff form。
    RnsPolynomial c0 = ct[0], c1 = ct[1];
    intt_negacyclic_inplace_lazy(c0);
    reduce_strict(c0);
    intt_negacyclic_inplace_lazy(c1);
    reduce_strict(c1);

    const u64 q = c0.modulus_at(0);
    const size_t n = c0.dimension();
    const size_t k = coeff_index;

    LweCt out;
    out.modulus = q;
    out.b = c0[0][k];
    out.a.resize(n);
    // (c1·s)_k = Σ_i a[i]·s_i, a[i] = c1_{k-i} (i≤k) 或 -c1_{k-i+n} (i>k)。
    for (size_t i = 0; i < n; i++) {
        if (i <= k) {
            out.a[i] = c1[0][k - i];
        } else {
            u64 v = c1[0][k - i + n];
            out.a[i] = (v == 0) ? 0 : (q - v); // 负循环取负
        }
    }
    return out;
}

} // namespace hehub
