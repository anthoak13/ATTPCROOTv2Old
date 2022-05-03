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

#include <Rtypes.h>   // for Int_t, Double_t, THashConsistencyHolder, ClassDef
#include <TObject.h>  // for TObject
#include <TVector3.h> // for TVector3

#include <stdio.h> // for size_t

#include <algorithm> // for max
#include <utility>   // for pair
#include <vector>    // for vector

class AtEvent;
class TBuffer;
class TClass;
class TMemberInspector;

class AtRansacMod : public TObject {
public:
   struct Cluster {
      double ClusterStrength;        //< strength
      size_t ClusterSize;            //< size
      double ClusterChi2;            //< Chi2
      std::vector<int> ClusterIndex; //< Indices
      std::vector<double> fitPar;
   };

   using AllClusters = std::vector<Cluster>;
   using PotentialModels = std::vector<std::unique_ptr<AtTrackModel>>;

protected:
   AtModelType fModelType{AtModelType::kLINE};
   AtRandomSample::SampleMethod fRandSamplMode{0}; //!
   const std::vector<AtHit> *fHitArray{nullptr};   //!

   // Set in constructor
   float fRANSACMaxIteration{500};
   float fRANSACMinPoints{30};
   float fRANSACThreshold{15};
   Int_t fLineDistThreshold{40};
   double fChargeThres{0};

   // Data structure
   std::vector<AtTrack> fTrackCand; // Candidate tracks
   AllClusters cluster_vector;

public:
   AtRansacMod();
   virtual ~AtRansacMod();

   // Behavior
   void CalcRANSACMod(AtEvent *event);

   // Getters
   std::vector<AtTrack> GetTrackCand() const { return fTrackCand; };
   inline AllClusters GetClusters() const { return cluster_vector; }

   // Setters
   void SetRanSamMode(AtRandomSample::SampleMethod mode) { fRandSamplMode = mode; };
   void SetModelType(AtModelType type) { fModelType = type; }
   void SetDistanceThreshold(Float_t threshold) { fRANSACThreshold = threshold; };
   void SetMinHitsLine(Int_t nhits) { fRANSACMinPoints = nhits; };
   void SetNumItera(Int_t niterations) { fRANSACMaxIteration = niterations; };
   void SetChargeThres(double value) { fChargeThres = value; };

protected:
   // Virtual behavior functions
   virtual int evaluateModel(AtTrackModel *model, const std::vector<int> &pointsToCheck);

   void Reset();
   void Solve();
   void doIteration(PotentialModels &IdxMod1);
   std::vector<int> getPointsInModel(const std::vector<int> &indexes, AtTrackModel *model);
   void removePoints(std::vector<int> &vectorToModify, const std::vector<int> &pointsToRemove);
   void Init(AtEvent *event);

   std::vector<AtTrack *> Clusters2Tracks(AllClusters NClusters, AtEvent *event);
   void SetCluster(const std::vector<int> samplesIdx, const double cost, const double Chi2, std::vector<double> fitPar);

   // For some reason these also are what actually set the tracks to be saved...
   void FindVertex(std::vector<AtTrack *> tracks);
   void FindVertexOneTrack(std::vector<AtTrack *> tracks);

   ClassDef(AtRansacMod, 2);
};

#endif
