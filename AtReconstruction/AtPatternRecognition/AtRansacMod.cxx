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

   auto points = AtRandomSample::SamplePoints(testModel->GetNumPoints(), *fHitArray, fRandSamplMode);
   testModel->ConstructModel(points);

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
      std::cout << "RansacMod tracks size : " << tracksSize << std::endl;
      for (Int_t ntrack = 0; ntrack < tracks.size(); ntrack++) {
         tracks.at(ntrack)->SetTrackID(ntrack);
         fTrackCand.push_back(*tracks.at(ntrack));
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

      LOG(debug) << "Saving track with " << clustersize << " hits";

      auto *track = new AtTrack();

      for (int j = 0; j < clustersize; j++)
         track->AddHit(hits.at(indicesCluster[j]));

      track->SetFitPar(NClusters[i].fitPar);
      track->SetMinimum(NClusters[i].ClusterChi2);
      track->SetNFree(clustersize - 6);

      tracks.push_back(track);

      // delete track;
   }

   return tracks;
}
