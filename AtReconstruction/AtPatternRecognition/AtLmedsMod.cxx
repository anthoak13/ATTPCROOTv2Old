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

std::pair<double, int> AtLmedsMod::evaluateModel(const std::vector<int> &pointsToCheck)
{
   std::vector<double> errorsVec;
   // Loop through point and if it is an inlier, then add the error**2 to weight
   for (auto index : pointsToCheck) {

      double error = fModel->DistanceToModel(index);
      error = error * error;
      if (error < (fRANSACThreshold * fRANSACThreshold))
         errorsVec.push_back(error);
   }
   return {GetMedian(errorsVec) / errorsVec.size(), errorsVec.size()};
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
