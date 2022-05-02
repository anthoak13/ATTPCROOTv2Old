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
   // pair.first is "goodness"
   // pair.second is indices defining the model
   using PotentialModels = std::vector<std::pair<double, std::vector<int>>>;

   enum class SampleMethod;

protected:
   std::unique_ptr<AtTrackModel> fModel{nullptr}; //!
   AtTrackModel::SampleMethod fRandSamplMode{0};

   // Set in constructor
   float fRANSACMaxIteration{500};
   float fRANSACMinPoints{30};
   float fRANSACThreshold{15};
   Int_t fLineDistThreshold{40};
   double fChargeThres{0};
   // float fRANSACPointThreshold{};
   // float fRANSACChargeThreshold{};

   TVector3 fVertex_1{-10000, -10000, -10000};
   TVector3 fVertex_2{-10000, -10000, -10000};
   TVector3 fVertex_mean;
   Double_t fVertexTime{-10000};
   std::pair<Int_t, Int_t> fVertex_tracks; // ID of the tracks that form the best vertex
   Int_t fVertexMode{};

   // Data structure
   std::vector<AtTrack> fTrackCand; // Candidate tracks
   AllClusters cluster_vector;

public:
   AtRansacMod();
   virtual ~AtRansacMod();

   // Behavior
   void CalcRANSACMod(AtEvent *event);

   // Getters
   // double GetAvCharge() const { return Avcharge; };
   TVector3 GetVertex1() const { return fVertex_1; };
   TVector3 GetVertex2() const { return fVertex_2; };
   Double_t GetVertexTime() const { return fVertexTime; };
   TVector3 GetVertexMean() const { return fVertex_mean; };
   std::vector<AtTrack> GetTrackCand() const { return fTrackCand; };
   inline AllClusters GetClusters() const { return cluster_vector; }

   // Setters
   // void SetAvCharge(double charge) { Avcharge = charge; };
   void SetRanSamMode(AtTrackModel::SampleMethod mode) { fRandSamplMode = mode; };
   void SetDistanceThreshold(Float_t threshold) { fRANSACThreshold = threshold; };
   void SetMinHitsLine(Int_t nhits) { fRANSACMinPoints = nhits; };
   void SetNumItera(Int_t niterations) { fRANSACMaxIteration = niterations; };
   void SetChargeThres(double value) { fChargeThres = value; };
   void SetVertexMod(Int_t mode) { fVertexMode = mode; };

protected:
   // Virtual behavior functions
   virtual std::pair<double, int> evaluateModel(const std::vector<int> &pointsToCheck);

   void Reset();
   void Solve();
   void doIteration(PotentialModels &IdxMod1);
   std::vector<int> getPointsInModel(const std::vector<int> &indexes);
   void removePoints(std::vector<int> &vectorToModify, const std::vector<int> &pointsToRemove);
   void Init(AtEvent *event);

   TVector3 ClosestPoint2Lines(TVector3 d1, TVector3 pt1, TVector3 d2, TVector3 pt2);
   std::vector<AtTrack *> Clusters2Tracks(AllClusters NClusters, AtEvent *event);
   void FindVertex(std::vector<AtTrack *> tracks);
   void FindVertexOneTrack(std::vector<AtTrack *> tracks);
   void SetCluster(const std::vector<int> samplesIdx, const double cost, const double Chi2, std::vector<double> fitPar);

   ClassDef(AtRansacMod, 1);
};

#endif
