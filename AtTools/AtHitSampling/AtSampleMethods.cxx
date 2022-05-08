#include "AtSampleMethods.h"

#include "AtChargeWeighted.h"
#include "AtGaussian.h"
#include "AtSample.h"
#include "AtUniform.h"
#include "AtWeightedGaussian.h"

/**
 * Static factory method for creating instances of AtSample
 *
 * Only returns a nullptr for sub-types of AtSample that do not require additional information.
 * @ingroup AtHitSampling
 */
std::unique_ptr<AtTools::AtSample> AtTools::AtSample::CreateSampler(SampleMethod method)
{
   switch (method) {
   case SampleMethod::kUniform: return std::make_unique<AtUniform>();
   case SampleMethod::kChargeWeighted: return std::make_unique<AtChargeWeighted>();
   default: return nullptr;
   }
}
