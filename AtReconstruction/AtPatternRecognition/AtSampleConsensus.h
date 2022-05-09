/*******************************************************************
// Basic Sample consensus class
// Based on https://github.com/jczamorac/Tracking_RANSAC (https://doi.org/10.1016/j.nima.2020.164899)
// Author: A.K. Anthony
********************************************************************/

#ifndef ATSAMPLECONSENSUS_H
#define ATSAMPLECONSENSUS_H

#include "AtPattern.h"
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

class AtSampleConsensus final {
private:
   using Estimators = AtSampleEstimator::Estimators;
   using PatternType = AtPatterns::PatternType;
   using SampleMethod = AtTools::AtSample::SampleMethod;
   using AtPattern = AtPatterns::AtPattern;
   using PatternPtr = std::unique_ptr<AtPattern>;
   using AtSamplePtr = std::unique_ptr<AtTools::AtSample>;

   PatternPtr fPatternProto; //< Prototype of pattern to find
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
   void SetPatternType(PatternPtr type) { fPatternProto = std::move(type); }
   void SetEstimator(Estimators estimator) { fEstimator = estimator; }

   void SetNumIterations(Int_t niterations) { fIterations = niterations; };
   void SetMinHitsPattern(Int_t nhits) { fMinPatternPoints = nhits; };
   void SetDistanceThreshold(Float_t threshold) { fDistanceThreshold = threshold; };
   void SetChargeThreshold(double value) { fChargeThres = value; };

private:
   PatternPtr GeneratePatternFromHits(const std::vector<AtHit> &hitArray);
   std::vector<AtHit> movePointsInPattern(AtPattern *pattern, std::vector<AtHit> &indexes);
   void SaveTrack(AtPattern *pattern, std::vector<AtHit> &indexes, AtPatternEvent *event);
   AtTrack CreateTrack(AtPattern *pattern, std::vector<AtHit> &indexes);
};

#endif
