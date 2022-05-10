/*******************************************************************
// Basic Sample consensus class
// Based on https://github.com/jczamorac/Tracking_RANSAC (https://doi.org/10.1016/j.nima.2020.164899)
// Author: A.K. Anthony
********************************************************************/

#ifndef ATSAMPLECONSENSUS_H
#define ATSAMPLECONSENSUS_H

#include "AtPattern.h"
#include "AtPatternTypes.h"
#include "AtSample.h"
#include "AtSampleEstimator.h"
#include "AtTrack.h" // for AtTrack

#include <stdio.h> // for size_t

#include <algorithm> // for max
#include <set>
#include <utility> // for pair
#include <vector>  // for vector
class AtEvent;
class AtPatternEvent;

/**
 * @defgroup SampleConsensus Sample Consensus
 *
 * @brief Group of classes and functions for sample consensus estimators.

 * AtSampleConsensus will perform a sample consensus calculation using an estimator (implemented estimaros are listed in
 *  Estimators) and an AtPattern.
 *
 */
namespace SampleConsensus {

/**
 * @brief Perform a sample consensus on a cloud of AtHits
 * @ingroup SampleConsensus
 *
 * Construct a sample consensus using an estimator, pattern type, and method for randomly sampling AtHit cloud.
 */
class AtSampleConsensus final {
private:
   using Estimators = SampleConsensus::Estimators;
   using PatternType = AtPatterns::PatternType;
   using SampleMethod = RandomSample::SampleMethod;
   using AtPattern = AtPatterns::AtPattern;
   using PatternPtr = std::unique_ptr<AtPattern>;
   using AtSamplePtr = std::unique_ptr<RandomSample::AtSample>;

   PatternType fPatternType; //< Type of pattern to find
   Estimators fEstimator;    //< Estimator to evaluate pattern
   AtSamplePtr fRandSampler; //< Sampling Method (defaults to uniform)

   float fIterations{500};       //< Number of interations of sample consensus
   float fMinPatternPoints{30};  //< Required number of points to form a pattern
   float fDistanceThreshold{15}; //< Distance a point must be from pattern to be an inlier

   /**
    * @brief Min charge for charge weighted fit.
    *
    * Minimum charge to include point in charge weighted fit. If set to -1 it
    * disables charge weighting when fitting the track.
    */
   double fChargeThres{-1};

public:
   AtSampleConsensus();
   AtSampleConsensus(Estimators estimator, PatternType patternType, SampleMethod sampleMethod);

   AtPatternEvent Solve(AtEvent *event);
   AtPatternEvent Solve(const std::vector<AtHit> &hitArray);

   void SetRandomSample(AtSamplePtr mode) { fRandSampler = std::move(mode); };
   void SetPatternType(PatternType type) { fPatternType = type; }
   void SetEstimator(Estimators estimator) { fEstimator = estimator; }

   void SetNumIterations(Int_t niterations) { fIterations = niterations; };
   void SetMinHitsPattern(Int_t nhits) { fMinPatternPoints = nhits; };
   void SetDistanceThreshold(Float_t threshold) { fDistanceThreshold = threshold; };
   void SetChargeThreshold(double value) { fChargeThres = value; };

private:
   PatternPtr GeneratePatternFromHits(const std::vector<AtHit> &hitArray);
   std::vector<AtHit> movePointsInPattern(AtPattern *pattern, std::vector<AtHit> &indexes);
   // void SaveTrack(AtPattern *pattern, std::vector<AtHit> &indexes, AtPatternEvent *event);
   AtTrack CreateTrack(AtPattern *pattern, std::vector<AtHit> &indexes);
};
} // namespace SampleConsensus
#endif
