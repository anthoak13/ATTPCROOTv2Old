#ifndef ATSAMPLEESTIMATOR_H
#define ATSAMPLEESTIMATOR_H

#include <vector>

class AtPattern;
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
 */
class AtSampleEstimator {

public:
   enum class Estimators;
   /**
    * @brief Evaluate how well model describes hits
    *
    * Checks all hits for consistancy with modeled pattern, and sets the model Chi2 based on
    * the estimator.
    *
    * @param[in] model Model to evaluate
    * @param[in] hits Hits to compare to model
    * @param[in] distThresh How close a point must be to be consistent with the model
    * @return Number of points consistent with model in hits
    */
   static int EvaluateModel(AtPattern *model, const std::vector<AtHit> &hits, double distThresh, Estimators estimator);
};

#endif //#ifndef ATSAMPLEESTIMATOR_H
