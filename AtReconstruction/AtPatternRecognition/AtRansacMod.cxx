#include "AtRansacMod.h"

#include "AtEvent.h" // for AtEvent
#include "AtHit.h"   // for AtHit
#include "AtModelFactory.h"

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

AtRansacMod::AtRansacMod() = default;
AtRansacMod::~AtRansacMod() = default;

void AtRansacMod::Init(AtEvent *event)
{
   Reset();
   fHitArray = &(event->GetHitArray());
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
int AtRansacMod::evaluateModel(AtTrackModel *model, const std::vector<int> &pointsToCheck)
{
   int nbInliers = 0;
   double weight = 0;

   for (auto index : pointsToCheck) {
      // double error = distanceToModel(index);
      auto &pos = fHitArray->at(index).GetPosition();
      double error = model->DistanceToModel(pos);
      error = error * error;
      if (error < (fRANSACThreshold * fRANSACThreshold)) {
         nbInliers++;
         weight += error;
      }
   }
   model->SetChi2(weight / nbInliers);
   return nbInliers;
}

void AtRansacMod::doIteration(PotentialModels &IdxModel)
{
   std::vector<int> remainIndex;
   for (size_t i = 0; i < fHitArray->size(); i++)
      remainIndex.push_back(i);

   if (remainIndex.size() < fRANSACMinPoints)
      return;

   auto testModel = AtModelFactory::CreateModel(fModelType);

   auto randPoints = sampleModelPoints(testModel->GetNumPoints(), fRandSamplMode);
   testModel->ConstructModel(randPoints);

   auto nInliers = evaluateModel(testModel.get(), remainIndex);

   // If the model is consistent with enough points, save it
   if (nInliers > fRANSACMinPoints) {
      LOG(debug) << "Adding model with nInliers: " << nInliers;
      IdxModel.push_back(std::move(testModel));
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
   sort(models.begin(), models.end(), [](const auto &a, const auto &b) { return a->GetChi2() < b->GetChi2(); });

   std::vector<int> remainIndex;
   for (size_t i = 0; i < fHitArray->size(); i++)
      remainIndex.push_back(i);

   // Loop through each model, and extract the points that fit each model
   for (const auto &model : models) {
      if (remainIndex.size() < fRANSACMinPoints)
         break;

      // fModel->ConstructModel(model.second);

      // Get a vector of every point that fits this model
      auto inliers = getPointsInModel(remainIndex, model.get());

      // If there are enough points that fit the model save it
      if (inliers.size() > fRANSACMinPoints) {
         std::vector<XYZPoint> pointsIn;
         for (auto index : inliers)
            pointsIn.push_back(fHitArray->at(index).GetPosition());
         auto cost = model->GetChi2();
         double chi2 = model->FitModel(pointsIn);
         SetCluster(inliers, cost, chi2, model->GetModelPar());
      }

      // Remove all the points that fit this model (even if it wasn't saved?)
      removePoints(remainIndex, inliers);
   }
}
void AtRansacMod::removePoints(std::vector<int> &toModify, const std::vector<int> &toRemove)
{
   std::vector<int> tempRemain;
   std::set_difference(toModify.begin(), toModify.end(), toRemove.begin(), toRemove.end(),
                       std::inserter(tempRemain, tempRemain.begin()));
   toModify = std::move(tempRemain);
}

std::vector<int> AtRansacMod::getPointsInModel(const std::vector<int> &indexes, AtTrackModel *model)
{
   std::vector<int> retVec;
   for (auto index : indexes) {
      auto pos = fHitArray->at(index).GetPosition();
      double error = model->DistanceToModel(pos);
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
   if (fHitArray->size() > fRANSACMinPoints) {
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

/**** Beginning of Sampling info ******/

/**
 * @brief Sample points required to define a model.
 *
 * Sample N points, where N is the minimum number of points required to define
 * the model being used (right now, that's a linear model so N = 2)
 *
 * @todo Generalize the functions being called so they return and arbitrary number
 * of points based on the requirements of the model
 *
 * @param[in] mode flag instructing what sampling algorithm to use
 * @return vector of indices of the points describing the model
 */
std::vector<XYZPoint> AtRansacMod::sampleModelPoints(int numPoints, SampleMethod mode = SampleMethod::kUniform)
{
   switch (mode) {
   case (SampleMethod::kUniform): return sampleUniform(numPoints);
   case (SampleMethod::kGaussian): return sampleGaussian(numPoints);
   case (SampleMethod::kWeighted): return sampleWeighted(numPoints);
   case (SampleMethod::kWeightedGaussian): return sampleWeightedGaussian(numPoints);
   default: return {};
   }
}

/**
 * @brief Sample two random points from indX
 */
std::vector<XYZPoint> AtRansacMod::sampleUniform(int numPoints)
{
   //-------Uniform sampling
   int ind1 = gRandom->Uniform(0, fHitArray->size());
   int ind2;
   do {
      ind2 = gRandom->Uniform(0, fHitArray->size());
   } while (ind1 == ind2);

   return {fHitArray->at(ind1).GetPosition(), fHitArray->at(ind2).GetPosition()};
}

/**
 * @brief Sample two points from indX
 *
 * The first point is sampled randomly. The second points is sampled according to
 * a gaussian distribition around the first point with a sigma of 30. If in 20 samples it
 * doesn't find a point close enough, it defaults to a uniform sample.
 *
 * @TODO We should probably set sigma so it scales with the size of the pad plane or make it
 * a tunable parameter.
 */
std::vector<XYZPoint> AtRansacMod::sampleGaussian(int numPoints)
{
   //--------Gaussian sampling
   double sigma = 30.0;
   double y = 0;
   double gauss = 0;
   int counter = 0;
   int p1 = gRandom->Uniform(0, fHitArray->size());
   int p2;
   auto &P1 = fHitArray->at(p1).GetPosition();

   do {
      p2 = gRandom->Uniform(0, fHitArray->size());
      auto &P2 = fHitArray->at(p2).GetPosition();
      auto dist = std::sqrt((P2 - P1).Mag2());

      gauss = 1.0 * exp(-1.0 * pow(dist / sigma, 2.0));
      y = (gRandom->Uniform(0, 1));
      counter++;
      if (counter > 20 && p2 != p1)
         break;
   } while (p2 == p1 || y > gauss);

   return {fHitArray->at(p1).GetPosition(), fHitArray->at(p2).GetPosition()};
}
/**
 * @ brief Sample two points based on charge
 */
std::vector<XYZPoint> AtRansacMod::sampleWeighted(int numPoints)
{
   //-------Weighted sampling
   auto Proba = getPDF();
   bool validSecondPoint = false;
   int counter = 0;
   int p1 = gRandom->Uniform(0, fHitArray->size());
   int p2 = p1;

   do {
      validSecondPoint = false;
      counter++;
      if (counter > 30 && p2 != p1)
         break;

      p2 = gRandom->Uniform(0, fHitArray->size());
      double TwiceAvCharge = 2 * fAvgCharge;

      if (Proba.size() == fHitArray->size())
         validSecondPoint = (Proba[p2] >= gRandom->Uniform(0, TwiceAvCharge));
      else
         validSecondPoint = true;

   } while (p2 == p1 || validSecondPoint == false);

   return {fHitArray->at(p1).GetPosition(), fHitArray->at(p2).GetPosition()};
}

/**
 * @brief Sample two points from indX
 *
 * The first point is sampled randomly by charge. The second points is sampled according to
 * a gaussian distribition around the first point with a sigma of 30. If in 20 samples it
 * doesn't find a point close enough, it defaults to a charge weighted sample.
 *
 * @TODO We should probably set sigma so it scales with the size of the pad plane or make it
 * a tunable parameter.
 */
std::vector<XYZPoint> AtRansacMod::sampleWeightedGaussian(int numPoints)
{
   //-------Weighted sampling + Gauss dist.
   auto Proba = getPDF();
   bool validSecondPoint = false;
   double sigma = 30.0;
   double y = 0;
   double gauss = 0;
   int counter = 0;

   int p1 = gRandom->Uniform(0, fHitArray->size());
   int p2 = p1;
   auto &P1 = fHitArray->at(p1).GetPosition();

   do {
      p2 = gRandom->Uniform(0, fHitArray->size());
      auto &P2 = fHitArray->at(p2).GetPosition();
      auto dist = std::sqrt((P2 - P1).Mag2());

      gauss = 1.0 * exp(-1.0 * pow(dist / sigma, 2));
      y = (gRandom->Uniform(0, 1));

      counter++;
      if (counter > 30 && p2 != p1)
         break;

      validSecondPoint = false;
      double TwiceAvCharge = 2 * fAvgCharge;

      if (Proba.size() == fHitArray->size())
         validSecondPoint = (Proba[p2] >= gRandom->Uniform(0, TwiceAvCharge));
      else
         validSecondPoint = true;

   } while (p2 == p1 || validSecondPoint == false || y > gauss);

   return {fHitArray->at(p1).GetPosition(), fHitArray->at(p2).GetPosition()};
}

std::vector<double> AtRansacMod::getPDF()
{
   double Tcharge = 0;
   for (const auto &hit : *fHitArray)
      Tcharge += hit.GetCharge();

   fAvgCharge = Tcharge / fHitArray->size();
   std::vector<double> w;
   if (Tcharge > 0)
      for (const auto &hit : *fHitArray)
         w.push_back(hit.GetCharge() / Tcharge);

   return w;
}
