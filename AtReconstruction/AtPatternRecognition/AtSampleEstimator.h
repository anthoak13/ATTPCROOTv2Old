#ifndef ATSAMPLEESTIMATOR_H
#define ATSAMPLEESTIMATOR_H

#include <vector>
namespace AtPatterns {
class AtPattern;
}
class AtHit;

/**
 * @defgroup AtSampleEstimator Estimators
 *
 * Group of classes and functions for sample consensus estimators. Should be used by including the file
 * AtEstimatorMethods.h
 */

/**
 * Static class for calling the correct estimator based on the enum
 * Enum definition and implementation is in AtEstimatorMethods.h
 * @ingroup AtSampleEstimator
 */
class AtSampleEstimator {

public:
   enum class Estimators;
   static int
   EvaluateModel(AtPatterns::AtPattern *model, const std::vector<AtHit> &hits, double distThresh, Estimators estimator);
};

#endif //#ifndef ATSAMPLEESTIMATOR_H
