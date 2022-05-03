#include "AtRansacMod.h"

#include "AtHit.h" // for AtHit

#include <Math/Point3D.h> // for PositionVector3D

/**
 * @brief Check the goodness of the model
 *
 * Evaluates how good the model is in comparison to the points in pointsToCheck.
 *
 * @param[in] hitArray AtHits to compare to the model
 * @param[in/out] model Model to check. Sets Chi2 of the model
 * @return thenumber of inliers defined by fDistanceThreshold
 */
int AtRansacMod::evaluateModel(AtTrackModel *model, const std::vector<AtHit> &hitArray)
{
   int nbInliers = 0;
   double weight = 0;

   for (const auto &hit : hitArray) {
      auto &pos = hit.GetPosition();
      double error = model->DistanceToModel(pos);
      error = error * error;
      if (error < (fDistanceThreshold * fDistanceThreshold)) {
         nbInliers++;
         weight += error;
      }
   }
   model->SetChi2(weight / nbInliers);
   return nbInliers;
}
