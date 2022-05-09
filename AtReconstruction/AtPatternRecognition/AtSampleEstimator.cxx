#include "AtSampleEstimator.h"

#include "AtEstimatorMethods.h"
#include "AtPattern.h"

#include <cmath>

int AtSampleEstimator::EvaluateModel(AtPatterns::AtPattern *model, const std::vector<AtHit> &hits, double distThresh,
                                     Estimators estimator = Estimators::kRANSAC)
{
   switch (estimator) {
   case (Estimators::kRANSAC): return AtEstimatorRansac::EvaluateModel(model, hits, distThresh);
   case (Estimators::kLMedS): return AtEstimatorLmeds::EvaluateModel(model, hits, distThresh);
   case (Estimators::kMLESAC): return AtEstimatorMlesac::EvaluateModel(model, hits, distThresh);
   default: return 0;
   }
}
