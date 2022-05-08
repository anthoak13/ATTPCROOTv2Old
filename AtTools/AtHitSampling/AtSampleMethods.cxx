#include "AtSampleMethods.h"

std::unique_ptr<AtTools::AtSample> AtTools::AtSample::CreateSampler(SampleMethod method)
{
   switch (method) {
   case SampleMethod::kUniform: return std::make_unique<AtUniform>();
   case SampleMethod::kChargeWeighted: return std::make_unique<AtChargeWeighted>();
   default: return nullptr;
   }
}
