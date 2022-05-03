/*******************************************************************
// Basic RANSAC Class                                              *
// Author: J.C. Zamora, jczamorac@gmail.com                        *
// University of Sao Paulo, 26-08-2020                             *
********************************************************************/

#ifndef AtSAMPLECONSENSUS_H
#define AtSAMPLECONSENSUS_H

#include "AtModelFactory.h"
#include "AtRandomSample.h"
#include "AtSampleEstimator.h"
#include "AtTrack.h" // for AtTrack
#include "AtTrackModel.h"

#include <stdio.h> // for size_t

#include <algorithm> // for max
#include <set>
#include <utility> // for pair
#include <vector>  // for vector
class AtEvent;
class AtPatternEvent;

class AtSampleConsensus final {
protected:
   AtModelType fModelType{AtModelType::kLINE};
   AtSampleEstimator::Estimator fSampleEstimator{AtSampleEstimator::Estimator::kRANSAC};
   AtRandomSample::SampleMethod fRandSamplMode{AtRandomSample::SampleMethod::kUniform};

   float fIterations{500};       //< Number of interations of sample consensus
   float fMinModelPoints{30};    //< Required number of points to form a model
   float fDistanceThreshold{15}; //< Distance a point must be from model to be an inlier

   /**
    * @brief Min charge for charge weighted fit.
    *
    * Minimum charge to include point in charge weighted fit. If set to -1 it
    * disables charge weighting when fitting the track.
    */
   double fChargeThres{-1};

public:
   AtSampleConsensus();
   virtual ~AtSampleConsensus();

   AtPatternEvent Solve(AtEvent *event);
   AtPatternEvent Solve(const std::vector<AtHit> &hitArray);

   void SetRanSamMode(AtRandomSample::SampleMethod mode) { fRandSamplMode = mode; };
   void SetModelType(AtModelType type) { fModelType = type; }
   void SetEstimator(AtSampleEstimator::Estimator estimator) { fSampleEstimator = estimator; }

   void SetNumIterations(Int_t niterations) { fIterations = niterations; };
   void SetMinHitsModel(Int_t nhits) { fMinModelPoints = nhits; };
   void SetDistanceThreshold(Float_t threshold) { fDistanceThreshold = threshold; };
   void SetChargeThreshold(double value) { fChargeThres = value; };

private:
   using ModelPtr = std::unique_ptr<AtTrackModel>;
   std::unique_ptr<AtTrackModel> GenerateModel(const std::vector<AtHit> &hitArray);

   std::vector<AtHit> movePointsInModel(AtTrackModel *model, std::vector<AtHit> &indexes);
   AtTrack CreateTrack(AtTrackModel *model, std::vector<AtHit> &indexes);

   void SaveTrack(AtTrackModel *model, std::vector<AtHit> &indexes, AtPatternEvent *event);
};

#endif
