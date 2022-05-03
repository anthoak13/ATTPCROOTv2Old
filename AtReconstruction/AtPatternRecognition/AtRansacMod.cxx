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
int AtRansacMod::evaluateModel(AtTrackModel *model, const std::vector<AtHit> &hitArray)
{
   int nbInliers = 0;
   double weight = 0;

   for (const auto &hit : hitArray) {
      auto &pos = hit.GetPosition();
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
   auto nInliers = evaluateModel(model.get(), hitArray);
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
   auto remainHits = hitArray;
   for (const auto &model : models) {
      if (remainHits.size() < fRANSACMinPoints)
         break;

      auto inlierHits = movePointsInModel(model.get(), remainHits);
      if (inlierHits.size() > fRANSACMinPoints)
         SaveTrack(model.get(), inlierHits);
   }
}
void AtRansacMod::SaveTrack(AtTrackModel *model, std::vector<AtHit> &inliers)
{
   LOG(debug) << "Saving model: " << fTrackCand.size();

   fTrackCand.emplace_back();
   AtTrack &track = fTrackCand.back();
   track.SetTrackID(fTrackCand.size() - 1);

   // Add inliers to our ouput track
   for (auto &hit : inliers) {
      std::cout << hit.GetHitID() << std::endl;
      track.AddHit(std::move(hit));
   }

   model->FitModel(inliers, false);
   track.SetFitPar(model->GetModelPar());
   track.SetMinimum(model->GetChi2());
   track.SetNFree(model->GetNFree());
}
std::vector<AtHit> AtRansacMod::movePointsInModel(AtTrackModel *model, std::vector<AtHit> &hits)
{

   std::vector<AtHit> retVec;
   auto itStartEqualRange = hits.end();

   for (auto it = hits.begin(); it != hits.end(); ++it) {

      double error = model->DistanceToModel(it->GetPosition());
      auto isInModel = (error * error) < (fRANSACThreshold * fRANSACThreshold);

      // Start of sub-vector with hits in model
      if (isInModel && itStartEqualRange == hits.end()) {
         itStartEqualRange = it;
         continue;
      }

      // End of sub-vector with hits in model.
      // Move hits in this range to retVec then delete the empty entries
      if (itStartEqualRange != hits.end() && !isInModel) {
         retVec.insert(retVec.end(), std::make_move_iterator(itStartEqualRange), std::make_move_iterator(it));
         hits.erase(itStartEqualRange, it);
         it = itStartEqualRange;
         itStartEqualRange = hits.end();
         continue;
      }
   }

   // If the last chunk of the array was in the model, move it and delete empty entries
   if (itStartEqualRange != hits.end()) {
      auto it = hits.end();
      retVec.insert(retVec.end(), std::make_move_iterator(itStartEqualRange), std::make_move_iterator(it));
      hits.erase(itStartEqualRange, it);
   }
   return retVec;
}
