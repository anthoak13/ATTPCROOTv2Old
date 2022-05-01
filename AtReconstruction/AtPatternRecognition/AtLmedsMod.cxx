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

constexpr auto cRED = "\033[1;31m";
constexpr auto cYELLOW = "\033[1;33m";
constexpr auto cNORMAL = "\033[0m";
constexpr auto cGREEN = "\033[1;32m";

using namespace std;

ClassImp(AtLmedsMod);

void AtLmedsMod::doIteration(std::vector<std::pair<double, int>> &IdxModel1,
                             std::vector<std::pair<double, int>> &IdxModel2)
{
   std::vector<int> remainIndex;
   for (size_t i = 0; i < vX.size(); i++)
      remainIndex.push_back(i);

   if (remainIndex.size() < fRANSACMinPoints)
      return;

   // Sample 2 points from all availible. Store the indices of the points in a std::pair
   auto Rsamples = sampleModelPoints(remainIndex, fRandSamplMode); // random sampling

   // Set two vectors defining a line between the points (Vs and Ps)
   setModel(Rsamples); // estimate the linear model

   int nbInliers = 0;
   std::vector<double> errorsVec;
   // Loop through point and if it is an inlier, then add the error**2 to weight
   for (auto index : remainIndex) {

      double error = distanceToModel(index);
      error = error * error;
      if (error < (fRANSACThreshold * fRANSACThreshold)) {
         nbInliers++;
         errorsVec.push_back(error);
      }
   }

   // If there are enough inliers for this to be a potential line, save the indices defining the
   // model with a "scale" indicating the "goodness" of the model
   // scale = median(error_i**2)/nPoints
   double med = GetMedian(errorsVec);
   if (nbInliers > fRANSACMinPoints) {
      // getting the best models
      double scale = med / nbInliers;
      IdxModel1.emplace_back(scale, Rsamples.first);
      IdxModel2.emplace_back(scale, Rsamples.second);
   }
}

double AtLmedsMod::GetMedian(std::vector<double> errvec)
{
   size_t vsize = errvec.size();

   if (vsize == 0) {
      return 0;
   } else {
      sort(errvec.begin(), errvec.end());
      if (vsize % 2 == 0) {
         return (errvec[vsize / 2 - 1] + errvec[vsize / 2]) / 2;
      } else {
         return errvec[vsize / 2];
      }
   }
}
