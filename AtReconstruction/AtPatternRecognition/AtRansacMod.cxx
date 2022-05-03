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
int AtRansacMod::evaluateModel(AtTrackModel *model, const std::vector<int> &pointsToCheck,
                               const std::vector<AtHit> &hitArray)
{
   int nbInliers = 0;
   double weight = 0;

   for (auto index : pointsToCheck) {
      auto &pos = hitArray.at(index).GetPosition();
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

std::unique_ptr<AtTrackModel> AtRansacMod::GenerateModel(const std::vector<AtHit> &hitArray)
{

   std::vector<int> remainIndex;
   for (size_t i = 0; i < hitArray.size(); i++)
      remainIndex.push_back(i);

   if (remainIndex.size() < fRANSACMinPoints) {
      return nullptr;
   }

   auto model = AtModelFactory::CreateModel(fModelType);

   auto points = AtRandomSample::SamplePoints(model->GetNumPoints(), hitArray, fRandSamplMode);
   model->ConstructModel(points);

   LOG(debug) << "Testing model" << std::endl;
   auto nInliers = evaluateModel(model.get(), remainIndex, hitArray);
   LOG(debug) << "Found " << nInliers << " inliers";

   // If the model is consistent with enough points, save it
   if (nInliers > fRANSACMinPoints) {
      LOG(debug) << "Adding model with nInliers: " << nInliers;
      return model;
   }

   return nullptr;
}

void AtRansacMod::Solve(AtEvent *event)
{
   if (event->IsGood())
      Solve(event->GetHitArray());
}

void AtRansacMod::Solve(const std::vector<AtHit> &hitArray)
{
   if (hitArray.size() < fRANSACMinPoints) {
      LOG(error) << "Not enough points to solve. Requires" << fRANSACMinPoints;
      return;
   }

   auto cmp = [](const ModelPtr &a, const ModelPtr &b) { return a->GetChi2() < b->GetChi2(); };
   auto models = std::set<ModelPtr, decltype(cmp)>(cmp);
   // SortedModels models(decltype(cmp));

   LOG(debug2) << "Generating models";
   for (int i = 0; i < fRANSACMaxIteration; i++) {
      auto model = GenerateModel(hitArray);
      if (model != nullptr)
         models.insert(std::move(model));
   }
   LOG(debug2) << "Created " << models.size() << " valid models.";

   std::vector<int> remainIndex;
   for (size_t i = 0; i < hitArray.size(); i++)
      remainIndex.push_back(i);

   // Loop through each model, and extract the points that fit each model
   for (const auto &model : models) {
      if (remainIndex.size() < fRANSACMinPoints)
         break;

      LOG(debug2) << "Examining model: " << fTrackCand.size();
      // Get indices of hits in this model
      auto inliers = getPointsInModel(remainIndex, model.get(), hitArray);

      // If there are enough points that fit the model save it as a track
      if (inliers.size() > fRANSACMinPoints) {
         LOG(debug) << "Saving model: " << fTrackCand.size();
         fTrackCand.emplace_back();
         AtTrack &track = fTrackCand.back();
         track.SetTrackID(fTrackCand.size() - 1);

         // Add inliers to our ouput track
         std::vector<XYZPoint> pointsIn;
         for (auto index : inliers) {
            track.AddHit(hitArray.at(index));
            pointsIn.push_back(hitArray.at(index).GetPosition());
         }

         model->FitModel(pointsIn);
         track.SetFitPar(model->GetModelPar());
         track.SetMinimum(model->GetChi2());
         track.SetNFree(model->GetNFree());
      }

      // Remove all the points that fit this model (even if it wasn't saved?)
      removePoints(remainIndex, inliers);
   }
}

std::vector<int>
AtRansacMod::getPointsInModel(const std::vector<int> &indexes, AtTrackModel *model, const std::vector<AtHit> &hitArray)
{
   std::vector<int> retVec;
   for (auto index : indexes) {
      auto pos = hitArray.at(index).GetPosition();
      double error = model->DistanceToModel(pos);
      if ((error * error) < (fRANSACThreshold * fRANSACThreshold))
         retVec.push_back(index);
   }
   return retVec;
}
