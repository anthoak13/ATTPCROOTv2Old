#include "AtRansacMod.h"

#include "AtEvent.h" // for AtEvent
#include "AtHit.h"   // for AtHit
#include "AtModelLine.h"

#include <Math/Point3D.h> // for PositionVector3D
#include <TMath.h>        // for Pi
#include <TRandom.h>      // for TRandom, gRandom

#include <cmath>    // for cos, sin, pow, sqrt, fabs, exp, acos, atan
#include <fstream>  // for std
#include <iterator> // for insert_iterator, inserter
#include <memory>   // for allocator_traits<>::value_type
constexpr auto cRED = "\033[1;31m";
constexpr auto cYELLOW = "\033[1;33m";
constexpr auto cNORMAL = "\033[0m";
constexpr auto cGREEN = "\033[1;32m";

using namespace std;

ClassImp(AtRansacMod);

enum class AtRansacMod::SampleMethod { kUniform = 0, kGaussian = 1, kWeighted = 2, kWeightedGaussian = 3 };

AtRansacMod::AtRansacMod() : fModel(std::make_unique<AtModelLine>()) {}
AtRansacMod::~AtRansacMod() = default;

void AtRansacMod::Init(AtEvent *event)
{
   Reset();

   auto hitArray = event->GetHitArray();
   fModel->SetHitArray(&(event->GetHitArray()));
}

void AtRansacMod::Reset()
{
   cluster_vector.clear();
}

/**
 * @brief Check the goodness of the model
 *
 * Evaluates how good the model is in comparison to the points in pointsToCheck.
 *
 * @param[in] pointsToCheck List of indices of points to compare to the model
 *
 * @return A pair where the the smaller the first number, the better the model
 * and the second element is the number of inliers defined by fRANSACThreshold
 */
std::pair<double, int> AtRansacMod::evaluateModel(const std::vector<int> &pointsToCheck)
{
   int nbInliers = 0;
   double weight = 0;

   for (auto index : pointsToCheck) {
      // double error = distanceToModel(index);
      double error = fModel->DistanceToModel(index);
      error = error * error;
      if (error < (fRANSACThreshold * fRANSACThreshold)) {
         nbInliers++;
         weight += error;
      }
   }
   return {weight / nbInliers, nbInliers};
}

void AtRansacMod::doIteration(PotentialModels &IdxModel)
{
   std::vector<int> remainIndex;
   for (size_t i = 0; i < fModel->GetNumHits(); i++)
      remainIndex.push_back(i);

   if (remainIndex.size() < fRANSACMinPoints)
      return;

   // Sample and set the model (rn linear)
   fModel->ConstructRandomModel(fRandSamplMode);

   // auto Rsamples = sampleModelPoints(remainIndex, fRandSamplMode);
   // setModel(Rsamples);

   auto evaluation = evaluateModel(remainIndex);

   // If the model is consistent with enough points, save it
   if (evaluation.second > fRANSACMinPoints) {
      LOG(debug) << "Adding model " << evaluation.first;
      double scale = evaluation.first;
      IdxModel.emplace_back(scale, fModel->GetIndices());
   }
}
void AtRansacMod::Solve()
{
   // Vectors to store the indices used to define the model and their
   // "scale" of the model

   PotentialModels models;

   // Populate vectors of randomly sampled potential models
   for (int i = 0; i < fRANSACMaxIteration; i++)
      doIteration(models);

   // sort clusters by "goodness" of the models
   sort(models.begin(), models.end());

   std::vector<int> remainIndex;
   for (size_t i = 0; i < fModel->GetNumHits(); i++)
      remainIndex.push_back(i);

   // Loop through each model, and extract the points that fit each model
   for (const auto &model : models) {
      if (remainIndex.size() < fRANSACMinPoints)
         break;
      fModel->ConstructModel(model.second);

      // Get a vector of every point that fits this model
      std::vector<int> inlIdxR = getPointsInModel(remainIndex);

      // If there are enough points that fit the model save it
      if (inlIdxR.size() > fRANSACMinPoints) {
         std::vector<double> fitPar;
         double chi2 = fModel->Fit3D(inlIdxR, fitPar);
         SetCluster(inlIdxR, model.first, chi2, fitPar);
      }

      // Remove all the points that fit this model (even if it wasn't saved?)
      removePoints(remainIndex, inlIdxR);
   }
}
void AtRansacMod::removePoints(std::vector<int> &toModify, const std::vector<int> &toRemove)
{
   std::vector<int> tempRemain;
   std::set_difference(toModify.begin(), toModify.end(), toRemove.begin(), toRemove.end(),
                       std::inserter(tempRemain, tempRemain.begin()));
   toModify = std::move(tempRemain);
}

std::vector<int> AtRansacMod::getPointsInModel(const std::vector<int> &indexes)
{
   std::vector<int> retVec;
   for (auto index : indexes) {
      double error = fModel->DistanceToModel(index);
      if ((error * error) < (fRANSACThreshold * fRANSACThreshold))
         retVec.push_back(index);
   }
   return retVec;
}

void AtRansacMod::CalcRANSACMod(AtEvent *event)
{

   // std::cout << "Inicializa" << '\n';
   Init(event);
   // std::cout << "Resuelve" << '\n';
   if (fModel->GetNumHits() > fRANSACMinPoints) {
      Solve();
      AllClusters myClusters = GetClusters();
      // I know this step is stupid, but this part is meant to be in the future a clustering process
      // std::cout << "Escribe tracks" << '\n';
      std::vector<AtTrack *> tracks = Clusters2Tracks(myClusters, event);

      Int_t tracksSize = tracks.size();
      // std::cout<<"RansacMod tracks size : "<<tracksSize<<std::endl;
      for (Int_t ntrack = 0; ntrack < tracks.size(); ntrack++) {
         tracks.at(ntrack)->SetTrackID(ntrack);
      }

      switch (fVertexMode) {
      case 1:
         if (tracksSize > 1)
            FindVertex(tracks);
         break;
      default:
         if (tracksSize > 0)
            FindVertexOneTrack(tracks); // find a vertex for each track
      }
   }
}

void AtRansacMod::SetCluster(const std::vector<int> samplesIdx, const double cost, const double Chi2,
                             std::vector<double> fitPar)
{
   Cluster cstr;
   cstr.ClusterIndex = samplesIdx;
   cstr.ClusterSize = samplesIdx.size();
   cstr.ClusterStrength = cost;
   cstr.ClusterChi2 = Chi2;
   cstr.fitPar = std::move(fitPar);
   cluster_vector.push_back(cstr);
}

std::vector<AtTrack *> AtRansacMod::Clusters2Tracks(AllClusters NClusters, AtEvent *event)
{

   std::vector<AtTrack *> tracks;

   auto hits = event->GetHitArray();

   int numclus = NClusters.size();
   // std::cout << "numero de clusters "<<numclus << '\n';

   for (int i = 0; i < numclus; i++) {
      size_t clustersize = NClusters[i].ClusterSize;
      std::vector<int> indicesCluster = NClusters[i].ClusterIndex;
      // double costo =  NClusters[i].ClusterStrength;
      double Chi2 = NClusters[i].ClusterChi2;
      const auto &fitPar = NClusters[i].fitPar;
      TVector3 punto1 = {fitPar[0], fitPar[1], fitPar[2]};
      TVector3 punto2 = {fitPar[3], fitPar[4], fitPar[5]};
      TVector3 pdiff = punto2 - punto1;

      auto *track = new AtTrack();

      for (int j = 0; j < clustersize; j++)
         track->AddHit(hits.at(indicesCluster[j]));

      std::vector<Double_t> par;
      par.push_back(punto1.X()); // 0
      par.push_back(pdiff.X());  // 1
      par.push_back(punto1.Y()); // 2
      par.push_back(pdiff.Y());  // 3
      par.push_back(punto1.Z()); // 4
      par.push_back(pdiff.Z());  // 5

      track->SetFitPar(par);
      track->SetMinimum(Chi2);
      track->SetNFree(clustersize - 6);

      tracks.push_back(track);

      // delete track;
   }

   return tracks;
}

void AtRansacMod::FindVertex(std::vector<AtTrack *> tracks)
{

   // Assumes the minimum distance between two lines, with respect a given threshold, the first vertex candidate. Then
   // evaluates the distance of each remaining line with respect to the others (vertex) to decide the particles of the
   // reaction.
   // std::cout<<" New find vertex call "<<std::endl;

   Double_t mad = 10; // Minimum approach distance. This is the minimum distance between the lines in 3D. Must be bigger
                      // than fLineDistThreshold
   // ROOT::Math::XYZVector c_1(-1000,-1000,-1000);
   // ROOT::Math::XYZVector c_2(-1000,-1000,-1000);
   TVector3 c_1(-1000, -1000, -1000);
   TVector3 c_2(-1000, -1000, -1000);
   // std::vector<AtTrack*> *TrackCand;

   // Current  parametrization
   // x = p[0] + p[1]*t;
   // y = p[2] + p[3]*t;
   // z = p[4] + p[5]*t;
   //  (x,y,z) = (p[0],p[2],p[4]) + (p[1],p[3],p[5])*t

   // Vector of the beam determined from the experimental data
   TVector3 BeamDir(0., 0., 1.0);
   TVector3 BeamPoint(0., 0., 500.0);

   std::vector<Bool_t> IsFilled;
   for (Int_t i = 0; i < int(tracks.size()); i++)
      IsFilled.push_back(kFALSE);

   // Test each line against the others to find a vertex candidate
   for (Int_t i = 0; i < int(tracks.size()) - 1; i++) {

      AtTrack *track = tracks.at(i);
      std::vector<Double_t> p = track->GetFitPar();

      if (p.size() == 0)
         continue;

      TVector3 p1(p[0], p[2], p[4]); // p1
      TVector3 e1(p[1], p[3], p[5]); // d1

      for (Int_t j = i + 1; j < tracks.size(); j++) {
         AtTrack *track_f = tracks.at(j);
         std::vector<Double_t> p_f = track_f->GetFitPar();
         if (track->GetTrackID() == track_f->GetTrackID())
            continue;

         if (p_f.size() == 0)
            continue;

         TVector3 p2(p_f[0], p_f[2], p_f[4]); // p2
         TVector3 e2(p_f[1], p_f[3], p_f[5]); // d2
         double angle = e1.Angle(e2) * 180. / 3.1415;

         TVector3 n = e1.Cross(e2);
         double sdist = fabs(n.Dot(p1 - p2) / n.Mag());
         TVector3 vertexbuff = ClosestPoint2Lines(e1, p1, e2, p2);

         if (sdist < fLineDistThreshold) {
            TVector3 meanVer = vertexbuff;
            double radius = sqrt(vertexbuff.X() * vertexbuff.X() + vertexbuff.Y() * vertexbuff.Y());

            if (vertexbuff.Z() > 0 && vertexbuff.Z() < 1000 && radius < 25.0) {

               track->SetTrackVertex(XYZPoint(vertexbuff));
               track_f->SetTrackVertex(XYZPoint(vertexbuff));
               fVertex_1.SetXYZ(vertexbuff.X(), vertexbuff.Y(), vertexbuff.Z());
               fVertex_2.SetXYZ(vertexbuff.X(), vertexbuff.Y(), vertexbuff.Z());
               fVertex_mean = vertexbuff;

               if (!IsFilled[track->GetTrackID()]) {
                  // std::cout<<" Add track  "<<track->GetTrackID()<<std::endl;
                  IsFilled[track->GetTrackID()] = kTRUE;
                  fTrackCand.push_back(*track);
               }

               if (!IsFilled[track_f->GetTrackID()]) {
                  // std::cout<<" Add trackf  "<<track_f->GetTrackID()<<std::endl;
                  IsFilled[track_f->GetTrackID()] = kTRUE;
                  fTrackCand.push_back(*track_f);
               }

               // condition for almost parallel tracks and vertex out of the chamber
            } else if (angle < 10 && angle > 170) {

               TVector3 vertexbuff1 = ClosestPoint2Lines(e1, p1, BeamDir, BeamPoint);
               TVector3 vertexbuff2 = ClosestPoint2Lines(e2, p2, BeamDir, BeamPoint);
               TVector3 vertexbuffav = 0.5 * (vertexbuff1 + vertexbuff2);
               TVector3 vertexbuffdif = vertexbuff1 - vertexbuff2;

               if (vertexbuffdif.Mag() > 2 * fLineDistThreshold)
                  continue;

               track->SetTrackVertex(XYZPoint(vertexbuffav));
               track_f->SetTrackVertex(XYZPoint(vertexbuffav));
               fVertex_1.SetXYZ(vertexbuffav.X(), vertexbuffav.Y(), vertexbuffav.Z());
               fVertex_2.SetXYZ(vertexbuffav.X(), vertexbuffav.Y(), vertexbuffav.Z());
               fVertex_mean = vertexbuffav;

               if (!IsFilled[track->GetTrackID()]) {
                  // std::cout<<" Add track  "<<track->GetTrackID()<<std::endl;
                  IsFilled[track->GetTrackID()] = kTRUE;
                  fTrackCand.push_back(*track);
               }

               if (!IsFilled[track_f->GetTrackID()]) {
                  // std::cout<<" Add trackf  "<<track_f->GetTrackID()<<std::endl;
                  IsFilled[track_f->GetTrackID()] = kTRUE;
                  fTrackCand.push_back(*track_f);
               }
            }

         } // if fLineDistThreshold

      } // End of track_f

   } // Loop over the tracks

   if (fTrackCand.size() > 5)
      fTrackCand.resize(5); // Truncate the maximum number of lines to 5
}

void AtRansacMod::FindVertexOneTrack(std::vector<AtTrack *> tracks)
{

   // Assumes the minimum distance between two lines, with respect a given threshold, the first vertex candidate. Then
   // evaluates the distance of each remaining line with respect to the others (vertex) to decide the particles of the
   // reaction.
   // std::cout<<" New find vertex call "<<std::endl;

   Double_t mad = 10; // Minimum approach distance. This is the minimum distance between the lines in 3D. Must be bigger
                      // than fLineDistThreshold
   // ROOT::Math::XYZVector c_1(-1000,-1000,-1000);
   // ROOT::Math::XYZVector c_2(-1000,-1000,-1000);
   TVector3 c_1(-1000, -1000, -1000);
   TVector3 c_2(-1000, -1000, -1000);
   // std::vector<AtTrack*> *TrackCand;

   // Current  parametrization
   // x = p[0] + p[1]*t;
   // y = p[2] + p[3]*t;
   // z = p[4] + p[5]*t;
   //  (x,y,z) = (p[0],p[2],p[4]) + (p[1],p[3],p[5])*t

   // Vector of the beam determined from the experimental data
   TVector3 BeamDir(0., 0., 1.0);
   TVector3 BeamPoint(0., 0., 500.0);

   std::vector<Bool_t> IsFilled;
   for (Int_t i = 0; i < int(tracks.size()); i++)
      IsFilled.push_back(kFALSE);

   // Test each line against the others to find a vertex candidate
   for (auto track : tracks) {

      std::vector<Double_t> p = track->GetFitPar();

      if (p.size() == 0)
         continue;

      TVector3 p1(p[0], p[2], p[4]); // p1
      TVector3 e1(p[1], p[3], p[5]); // d1

      // double angle = e1.Angle(BeamDir) * 180. / 3.1415;

      TVector3 n = e1.Cross(BeamDir);
      double sdist = fabs(n.Dot(p1 - BeamPoint) / n.Mag());
      TVector3 vertexbuff = ClosestPoint2Lines(e1, p1, BeamDir, BeamPoint);

      if (sdist < fLineDistThreshold) {
         TVector3 meanVer = vertexbuff;
         double radius = sqrt(vertexbuff.X() * vertexbuff.X() + vertexbuff.Y() * vertexbuff.Y());

         // if(vertexbuff.Z()>0 && vertexbuff.Z()<1000 && radius<25.0){
         if (radius < 25.0) {
            track->SetTrackVertex(XYZPoint(vertexbuff));
            fVertex_1.SetXYZ(vertexbuff.X(), vertexbuff.Y(), vertexbuff.Z());
            fVertex_mean = vertexbuff;

            if (!IsFilled[track->GetTrackID()]) {
               // std::cout<<" Add track  "<<track->GetTrackID()<<std::endl;
               IsFilled[track->GetTrackID()] = kTRUE;
               fTrackCand.push_back(*track);
            }

            // condition for almost parallel tracks and vertex out of the chamber
         }

      } // if fLineDistThreshold

   } // Loop over the tracks

   if (fTrackCand.size() > 5)
      fTrackCand.resize(5); // Truncate the maximum number of lines to 5
}

TVector3 AtRansacMod::ClosestPoint2Lines(TVector3 d1, TVector3 pt1, TVector3 d2, TVector3 pt2)
{
   TVector3 n1 = d1.Cross(d2.Cross(d1));
   TVector3 n2 = d2.Cross(d1.Cross(d2));
   double t1 = (pt2 - pt1).Dot(n2) / (d1.Dot(n2));
   double t2 = (pt1 - pt2).Dot(n1) / (d2.Dot(n1));
   TVector3 c1 = pt1 + t1 * d1;
   TVector3 c2 = pt2 + t2 * d2;
   TVector3 meanpoint = 0.5 * (c1 + c2);

   return meanpoint;
}
