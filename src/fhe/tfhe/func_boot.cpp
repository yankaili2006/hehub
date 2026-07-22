#include "func_boot.h"
#include "blind_rotate.h"

namespace hehub {
namespace tfhe {

LweCt functional_bootstrap(const LweCt &ct, const RnsPolynomial &lut_poly,
                           const std::vector<RgswCt> &bootstrap_keys) {
    // 参数由 LUT 多项式的维度与模数推导 (与自举密钥生成侧一致的固定基位宽)。
    const size_t N = lut_poly.dimension();
    const u64 q = lut_poly.modulus_at(0);
    const TfheParams gp = make_tfhe_params(N, q);

    // 模切换 q→2N, 盲旋转, 再从累加器提取常数项对应的 LWE。
    const LweCt ct_2n = lwe_mod_switch(ct, 2 * N);
    const RlweCt acc = blind_rotate(ct_2n, lut_poly, bootstrap_keys, gp);
    return sample_extract(acc, 0);
}

} // namespace tfhe
} // namespace hehub
