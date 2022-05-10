#include "AtSampleMethods.h"

#include "AtChargeWeighted.h"
#include "AtGaussian.h"
#include "AtSample.h"
#include "AtUniform.h"
#include "AtWeightedGaussian.h"
using namespace RandomSample;

/**
 * @brief Static factory method for creating instances of AtSample.
 *
 * Factory for creating sampling methods. Will pass parameters to the constructors
 *Only returns a nullptr for sub-types of AtSample that do not require additional information.
 * @ingroup AtHitSampling
 */
template <typename... Ts>
std::unique_ptr<AtSample> CreateSampler(SampleMethod method, Ts &&...params)
{
   switch (method) {
   case SampleMethod::kUniform: return std::make_unique<AtUniform>();
   case SampleMethod::kChargeWeighted: return std::make_unique<AtChargeWeighted>();
   case SampleMethod::kGaussian: return std::make_unique<AtGaussian>(std::forward<Ts>(params)...);
   case SampleMethod::kWeightedGaussian: return std::make_unique<AtGaussian>(std::forward<Ts>(params)...);
   default: return nullptr;
   }
}
