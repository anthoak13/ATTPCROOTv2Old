#include "AtLmedsMod.h"

#include "AtEvent.h" // for AtEvent
#include "AtHit.h"   // for AtHit

#include <Math/Point3D.h> // for PositionVector3D
#include <TMath.h>        // for Pi
#include <TRandom.h>      // for TRandom, gRandom

#include <cmath>    // for cos, sin, pow, sqrt, fabs, exp, acos, atan
#include <fstream>  // for std
#include <iterator> // for insert_iterator, inserter
#include <memory>   // for allocator_traits<>::value_type

int AtLmedsMod::evaluateModel(AtTrackModel *model, const std::vector<AtHit> &hitArray)
{
   std::vector<double> errorsVec;
   // Loop through point and if it is an inlier, then add the error**2 to weight
   for (const auto &hit : hitArray) {

      double error = model->DistanceToModel(hit.GetPosition());
      error = error * error;
      if (error < (fDistanceThreshold * fDistanceThreshold))
         errorsVec.push_back(error);
   }
   model->SetChi2(GetMedian(errorsVec) / errorsVec.size());
   return errorsVec.size();
}

/**
 * @brief Sort vec and returnt the median
 */
double AtLmedsMod::GetMedian(std::vector<double> &vec)
{
   auto vsize = vec.size();

   if (vsize == 0) {
      return 0;
   } else {
      sort(vec.begin(), vec.end());
      if (vsize % 2 == 0) {
         return (vec[vsize / 2 - 1] + vec[vsize / 2]) / 2;
      } else {
         return vec[vsize / 2];
      }
   }
}
