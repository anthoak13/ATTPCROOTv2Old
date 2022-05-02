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
      TVector3 ClusterFitP1;         //< point 1 from the fitted line
      TVector3 ClusterFitP2;         //< point 2 from the fitted line
   };

   using AllClusters = std::vector<Cluster>;
   enum class SampleMethod;

protected:
   // Set in constructor
   float fRANSACMaxIteration{500};
   float fRANSACMinPoints{30};
   float fRANSACThreshold{15};
   Int_t fLineDistThreshold{40};
   SampleMethod fRandSamplMode{0};

   TVector3 fVertex_1{-10000, -10000, -10000};
   TVector3 fVertex_2{-10000, -10000, -10000};
   Double_t fVertexTime{-10000};
   Double_t fMinimum{-1};
   double fChargeThres{0};

   // Not set in constructor
   TVector3 fVertex_mean;
   float fRANSACPointThreshold{};
   float fRANSACChargeThreshold{};
   std::vector<AtTrack> fTrackCand;        // Candidate tracks
   std::pair<Int_t, Int_t> fVertex_tracks; // ID of the tracks that form the best vertex

   Int_t fVertexMod{};

   std::vector<double> vX, vY, vZ, vQ;
   std::vector<double> vTrackCharge;
   int fNumberOfTracksMax{};
   int fOriginalCloudSize{};
   double fTotalCharge{};
   int fVerbose{};
   double Avcharge{};

   TVector3 Vs;
   TVector3 Ps;
   AllClusters cluster_vector;

public:
   AtRansacMod();
   virtual ~AtRansacMod();

   // Behavior
   void CalcRANSACMod(AtEvent *event);

   // Getters
   double GetAvCharge() const { return Avcharge; };
   TVector3 GetVertex1() const { return fVertex_1; };
   TVector3 GetVertex2() const { return fVertex_2; };
   Double_t GetVertexTime() const { return fVertexTime; };
   TVector3 GetVertexMean() const { return fVertex_mean; };
   std::vector<AtTrack> GetTrackCand() const { return fTrackCand; };
   inline AllClusters GetClusters() const { return cluster_vector; }

   // Setters
   void SetAvCharge(double charge) { Avcharge = charge; };
   void SetRanSamMode(SampleMethod mode) { fRandSamplMode = mode; };
   void SetDistanceThreshold(Float_t threshold) { fRANSACThreshold = threshold; };
   void SetMinHitsLine(Int_t nhits) { fRANSACMinPoints = nhits; };
   void SetNumItera(Int_t niterations) { fRANSACMaxIteration = niterations; };
   void SetChargeThres(double value) { fChargeThres = value; };
   void SetVertexMod(Int_t mode) { fVertexMod = mode; };

protected:
   // Virtual behavior functions
   virtual std::pair<double, int> evaluateModel(const std::vector<int> &pointsToCheck);

   void Reset();
   void Solve();
   void doIteration(std::vector<std::pair<double, int>> &IdxMod1, std::vector<std::pair<double, int>> &IdxMod2);
   std::vector<int> getPointsInModel(const std::vector<int> &indexes);
   void removePoints(std::vector<int> &vectorToModify, const std::vector<int> &pointsToRemove);
   void Init(AtEvent *event);

   /***** Begining of model-specific methods *****/
   // This pair<int, int> that is being passed around will be a member of the
   // model (and probably a vector not a pair for models that require more than two
   // points to define. Will also have to be an XYZPoint instead of an index)
   std::pair<int, int> sampleModelPoints(std::vector<int> indX, SampleMethod mode);
   std::pair<int, int> sampleUniform(const std::vector<int> &indX);
   std::pair<int, int> sampleGaussian(const std::vector<int> &indX);
   std::pair<int, int> sampleWeighted(const std::vector<int> &indX);
   std::pair<int, int> sampleWeightedGaussian(const std::vector<int> &indX);
   void setModel(const std::pair<int, int> samplesIdx);
   double Fit3D(std::vector<int> inliners, TVector3 &V1, TVector3 &V2);
   double distanceToModel(int i);
   /***** End of model-specific methods ******/

   /***** Begining of model methods (in base class) *****/
   std::vector<double> GetPDF(const std::vector<int> samplesIdx);
   /***** End of model methods (in base class) *****/

   TVector3 ClosestPoint2Lines(TVector3 d1, TVector3 pt1, TVector3 d2, TVector3 pt2);
   std::vector<AtTrack *> Clusters2Tracks(AllClusters NClusters, AtEvent *event);
   void FindVertex(std::vector<AtTrack *> tracks);
   void FindVertexOneTrack(std::vector<AtTrack *> tracks);
   void SetCluster(const std::vector<int> samplesIdx, const double cost, const double Chi2, TVector3 CP1, TVector3 CP2);

   ClassDef(AtRansacMod, 1);
};

#endif
