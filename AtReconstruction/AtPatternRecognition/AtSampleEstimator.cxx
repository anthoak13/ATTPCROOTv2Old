#include "AtSampleEstimator.h"

#include "AtEstimatorMethods.h"
#include "AtPattern.h"

#include <cmath>

/**
 * @brief Evaluate how well model describes hits
 *
 * Checks all hits for consistancy with modeled pattern, and sets the model Chi2 based on
 * the estimator.
 *
 * @ingroup AtSampleEstimator
 *
 * @param[in] model Model to evaluate
 * @param[in] hits Hits to compare to model
 * @param[in] distThresh How close a point must be to be consistent with the model
 * @return Number of points consistent with model in hits
 */

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
