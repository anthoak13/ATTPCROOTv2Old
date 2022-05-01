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

void AtLmedsMod::Reset()
{
   AtRansacMod::Reset();
   errorsVec.clear();
}
void AtLmedsMod::Solve()
{

   // std::cout << "numero de puntos  "<<vX.size()<< '\n';
   std::vector<int> remainIndex;
   for (size_t i = 0; i < vX.size(); i++)
      remainIndex.push_back(i);

   TVector3 V1, V2;
   std::vector<int> inliners;
   inliners.clear();
   std::vector<std::pair<double, int>> IdxMod1;
   std::vector<std::pair<double, int>> IdxMod2;

   for (int i = 0; i < fRANSACMaxIteration; i++) {

      if (remainIndex.size() < fRANSACMinPoints)
         break;

      auto Rsamples = sampleModelPoints(remainIndex, fRandSamplMode); // random sampling
      setModel(Rsamples);                                             // estimate the linear model

      // std::vector<int> inlIdxR;
      int nbInliers = 0;

      for (int &j : remainIndex) {

         double error = distanceToModel(j); // error of each point relative to the model
         error = error * error;

         if (error < (fRANSACThreshold * fRANSACThreshold)) {
            // inlIdxR.push_back(*j);
            nbInliers++;
            errorsVec.push_back(error);
         }
      }

      double med = GetMedian(errorsVec);
      errorsVec.clear();

      if (nbInliers > fRANSACMinPoints) {
         // getting the best models
         double scale = med / nbInliers;
         IdxMod1.emplace_back(scale, Rsamples.first);
         IdxMod2.emplace_back(scale, Rsamples.second);
      } // if a cluster was found

   } // for Lmeds interactions

   // sort clusters
   sort(IdxMod1.begin(), IdxMod1.end());
   sort(IdxMod2.begin(), IdxMod2.end());

   remainIndex.clear(); // track remaining points
   for (size_t i = 0; i < vX.size(); i++)
      remainIndex.push_back(i);

   // extract inliers using the models
   for (int i = 0; i < IdxMod1.size(); ++i) {
      std::pair<int, int> ModInx = {IdxMod1[i].second, IdxMod2[i].second};
      setModel(ModInx);
      std::vector<int> inlIdxR;

      if (remainIndex.size() < fRANSACMinPoints)
         break;

      int counter = 0;

      for (int &j : remainIndex) {
         double error = distanceToModel(j);

         if ((error * error) < (fRANSACThreshold * fRANSACThreshold)) {
            inlIdxR.push_back(j);
            counter++;
         }
      }

      if (counter > fRANSACMinPoints) {
         TVector3 v1, v2;
         double chi2 = Fit3D(inlIdxR, v1, v2);
         SetCluster(inlIdxR, IdxMod1[i].first, chi2, v1, v2);
         v1.Clear();
         v2.Clear();
      }
      std::vector<int> tempRemain;
      std::set_difference(remainIndex.begin(), remainIndex.end(), inlIdxR.begin(), inlIdxR.end(),
                          std::inserter(tempRemain, tempRemain.begin()));
      remainIndex = tempRemain;
      inlIdxR.clear();
      tempRemain.clear();
   }

   IdxMod1.clear();
   IdxMod2.clear();
   remainIndex.clear();
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
