/*******************************************************************
// Basic Sample consensus class
// Based on https://github.com/jczamorac/Tracking_RANSAC (https://doi.org/10.1016/j.nima.2020.164899)
// Author: A.K. Anthony
********************************************************************/

#ifndef ATSAMPLECONSENSUS_H
#define ATSAMPLECONSENSUS_H

#include "AtPattern.h"
#include "AtRandomSample.h"
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
   using Estimator = AtSampleEstimator::Estimator;
   using SampleMethod = AtRandomSample::SampleMethod;
   using PatternType = AtPattern::Type;
   using PatternPtr = std::unique_ptr<AtPattern>;

   PatternType fPatternType;                            //< Pattern to find
   Estimator fEstimator{Estimator::kRANSAC};            //< Estimator to evaluate pattern
   SampleMethod fRandSamplMode{SampleMethod::kUniform}; //< Sampling method

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

   AtPatternEvent Solve(AtEvent *event);
   AtPatternEvent Solve(const std::vector<AtHit> &hitArray);

   void SetRanSamMode(SampleMethod mode) { fRandSamplMode = mode; };
   void SetPatternType(PatternType type) { fPatternType = type; }
   void SetEstimator(Estimator estimator) { fEstimator = estimator; }

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
