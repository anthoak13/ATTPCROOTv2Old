#ifndef ATESTIMATORMETHODS_H
#define ATESTIMATORMETHODS_H

#include "AtSampleEstimator.h"

/**
 * @brief Implemented estimators for sample consensus
 */
enum class AtSampleEstimator::Estimators { kRANSAC, kLMedS, kMLESAC };

/**
 * Implementation of RANSAC estimator
 * @ingroup AtSampleEstimator
 */
class AtEstimatorRansac {
public:
   static int EvaluateModel(AtPatterns::AtPattern *model, const std::vector<AtHit> &hitArray, double distanceThreshold);
};

/**
 * Implementation of MLESAC estimator
 * @ingroup AtSampleEstimator
 */

class AtEstimatorMlesac {
public:
   static int EvaluateModel(AtPatterns::AtPattern *model, const std::vector<AtHit> &hitArray, double distanceThreshold);
};

/**
 * Implementation of LMedS estimator
 * @ingroup AtSampleEstimator
 */
class AtEstimatorLmeds {
public:
   static int EvaluateModel(AtPatterns::AtPattern *model, const std::vector<AtHit> &hitArray, double distanceThreshold);
   static double GetMedian(std::vector<double> &vec);
};

#endif //#ifndef ATESTIMATORMETHODS_H
