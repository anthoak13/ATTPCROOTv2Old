#ifndef ATSAMPLEESTIMATOR_H
#define ATSAMPLEESTIMATOR_H

#include "AtHit.h"

class AtTrackModel;

class AtSampleEstimator {

public:
   enum class Estimator { kRANSAC, kLMedS, kMLESAC };
   /**
    * @brief Evaluate how well model describes hits
    *
    * Checks all hits for consistancy with model, and sets the model Chi2 based on
    * the estimator.
    *
    * @param[in] model Model to evaluate
    * @param[in] hits Hits to compare to model
    * @param[in] distThresh How close a point must be to be consistent with the model
    * @return Number of points consistent with model in hits
    */
   static int EvaluateModel(AtTrackModel *model, const std::vector<AtHit> &hits, double distThresh,
                            Estimator estimator = Estimator::kRANSAC);
};

#endif //#ifndef ATSAMPLEESTIMATOR_H
