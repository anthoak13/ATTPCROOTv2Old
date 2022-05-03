/*******************************************************************
// Basic RANSAC Class                                              *
// Author: J.C. Zamora, jczamorac@gmail.com                        *
// University of Sao Paulo, 26-08-2020                             *
********************************************************************/

#ifndef AtRANSACMOD_H
#define AtRANSACMOD_H

#ifdef _OPENMP
#include <omp.h>
#endif

#include "AtModelFactory.h"
#include "AtRandomSample.h"
#include "AtTrack.h" // for AtTrack
#include "AtTrackModel.h"

#include <Rtypes.h> // for Int_t, Double_t, THashConsistencyHolder, ClassDef

#include <stdio.h> // for size_t

#include <algorithm> // for max
#include <set>
#include <utility> // for pair
#include <vector>  // for vector

class AtEvent;
class AtPatternEvent;

class AtRansacMod {
protected:
   AtModelType fModelType{AtModelType::kLINE};
   AtRandomSample::SampleMethod fRandSamplMode;

   float fRANSACMaxIteration{500}; //< Number of interations of ransac
   float fRANSACMinPoints{30};     //< Required number of points to form a track
   float fRANSACThreshold{15};     //< Distance
   // double fChargeThres{0};

public:
   AtRansacMod();
   virtual ~AtRansacMod();

   void Solve(AtEvent *event, AtPatternEvent *patternEvent);
   void Solve(const std::vector<AtHit> &hitArray, AtPatternEvent *event);

   void SetRanSamMode(AtRandomSample::SampleMethod mode) { fRandSamplMode = mode; };
   void SetModelType(AtModelType type) { fModelType = type; }
   void SetDistanceThreshold(Float_t threshold) { fRANSACThreshold = threshold; };
   void SetMinHitsLine(Int_t nhits) { fRANSACMinPoints = nhits; };
   void SetNumItera(Int_t niterations) { fRANSACMaxIteration = niterations; };
   // void SetChargeThres(double value) { fChargeThres = value; };

protected:
   using ModelPtr = std::unique_ptr<AtTrackModel>;
   std::unique_ptr<AtTrackModel> GenerateModel(const std::vector<AtHit> &hitArray);

   virtual int evaluateModel(AtTrackModel *model, const std::vector<AtHit> &hitsToCheck);
   std::vector<AtHit> movePointsInModel(AtTrackModel *model, std::vector<AtHit> &indexes);
   AtTrack CreateTrack(AtTrackModel *model, std::vector<AtHit> &indexes);

   void SaveTrack(AtTrackModel *model, std::vector<AtHit> &indexes, AtPatternEvent *event);
};

#endif
