#ifndef ATSAMPLEMETHODS_H
#define ATSAMPLEMETHODS_H

#include <memory>

namespace RandomSample {
/**
 * @brief. Methods of random sampling.
 *
 * All methods implemented that can be constructed by the factory method CreateSampler(SampleMethod).
 * @ingroup AtHitSampling
 */
enum class SampleMethod { kUniform = 0, kChargeWeighted = 1, kGaussian = 2, kWeightedGaussian = 3 };
class AtSample;

template <typename... Ts>
std::unique_ptr<AtSample> CreateSampler(SampleMethod method, Ts &&...params);

} // namespace RandomSample

#endif // #ifndef ATSAMPLEMETHODS_H
