#include "AtSampleConsensus.h"

#include "AtEvent.h" // for AtEvent
#include "AtHit.h"   // for AtHit
#include "AtModelFactory.h"
#include "AtPatternEvent.h"
#include "AtRandomSample.h"

#include <Math/Point3D.h> // for PositionVector3D
#include <TMath.h>        // for Pi
#include <TRandom.h>      // for TRandom, gRandom

#include <cmath>    // for cos, sin, pow, sqrt, fabs, exp, acos, atan
#include <fstream>  // for std
#include <iterator> // for insert_iterator, inserter
#include <memory>   // for allocator_traits<>::value_type

AtSampleConsensus::AtSampleConsensus() = default;
AtSampleConsensus::~AtSampleConsensus() = default;

std::unique_ptr<AtTrackModel> AtSampleConsensus::GenerateModel(const std::vector<AtHit> &hitArray)
{

   if (hitArray.size() < fMinModelPoints) {
      return nullptr;
   }

   auto model = AtModelFactory::CreateModel(fModelType);

   auto points = AtRandomSample::SamplePoints(model->GetNumPoints(), hitArray, fRandSamplMode);
   model->ConstructModel(points);

   LOG(debug) << "Testing model" << std::endl;
   auto nInliers = AtSampleEstimator::EvaluateModel(model.get(), hitArray, fDistanceThreshold, fSampleEstimator);
   LOG(debug) << "Found " << nInliers << " inliers";

   // If the model is consistent with enough points, save it
   if (nInliers > fMinModelPoints) {
      LOG(debug) << "Adding model with nInliers: " << nInliers;
      return model;
   }

   return nullptr;
}

AtPatternEvent AtSampleConsensus::Solve(AtEvent *event)
{
   if (event->IsGood())
      return Solve(event->GetHitArray());
   return AtPatternEvent();
}

AtPatternEvent AtSampleConsensus::Solve(const std::vector<AtHit> &hitArray)
{
   if (hitArray.size() < fMinModelPoints) {
      LOG(error) << "Not enough points to solve. Requires" << fMinModelPoints;
   }

   auto cmp = [](const ModelPtr &a, const ModelPtr &b) { return a->GetChi2() < b->GetChi2(); };
   auto sortedModels = std::set<ModelPtr, decltype(cmp)>(cmp);

   LOG(debug2) << "Generating " << fIterations << " models";
   for (int i = 0; i < fIterations; i++) {
      auto model = GenerateModel(hitArray);
      if (model != nullptr)
         sortedModels.insert(std::move(model));
   }
   LOG(debug2) << "Created " << sortedModels.size() << " valid models.";

   // Loop through each model, and extract the points that fit each model
   auto remainHits = hitArray;
   AtPatternEvent retEvent;
   for (const auto &model : sortedModels) {
      if (remainHits.size() < fMinModelPoints)
         break;

      auto inlierHits = movePointsInModel(model.get(), remainHits);
      if (inlierHits.size() > fMinModelPoints) {
         auto track = CreateTrack(model.get(), inlierHits);
         track.SetTrackID(retEvent.GetTrackCand().size());
         retEvent.AddTrack(track);
      }
   }
   return retEvent;
}

AtTrack AtSampleConsensus::CreateTrack(AtTrackModel *model, std::vector<AtHit> &inliers)
{
   AtTrack track;

   // Add inliers to our ouput track
   for (auto &hit : inliers)
      track.AddHit(std::move(hit));

   model->FitModel(inliers, fChargeThres);

   track.SetFitPar(model->GetModelPar());
   track.SetMinimum(model->GetChi2());
   track.SetNFree(model->GetNFree());
   return track;
}
/**
 * Moves the entries of hits that are consistent with the model to the returned vector
 *
 * @param [in] model
 * @param[in/out] hits Hits returned are removed from this vector
 * @return vector containing the AtHits consistent with the model
 *
 */
std::vector<AtHit> AtSampleConsensus::movePointsInModel(AtTrackModel *model, std::vector<AtHit> &hits)
{

   std::vector<AtHit> retVec;
   auto itStartEqualRange = hits.end();

   for (auto it = hits.begin(); it != hits.end(); ++it) {

      double error = model->DistanceToModel(it->GetPosition());
      auto isInModel = (error * error) < (fDistanceThreshold * fDistanceThreshold);

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
