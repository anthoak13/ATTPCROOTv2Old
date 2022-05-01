#include "AtMlesacMod.h"

#include "AtEvent.h" // for AtEvent
#include "AtHit.h"   // for AtHit

#include <Math/Point3D.h> // for PositionVector3D
#include <TMath.h>        // for Pi
#include <TRandom.h>      // for TRandom, gRandom

#include <cmath>    // for cos, sin, pow, sqrt, exp, fabs, acos, atan
#include <fstream>  // for std
#include <iterator> // for insert_iterator, inserter
#include <memory>   // for allocator_traits<>::value_type

constexpr auto cRED = "\033[1;31m";
constexpr auto cYELLOW = "\033[1;33m";
constexpr auto cNORMAL = "\033[0m";
constexpr auto cGREEN = "\033[1;32m";

using namespace std;

ClassImp(AtMlesacMod);

void AtMlesacMod::Solve()
{

   // std::cout << "numero de puntos  "<<vX.size()<< '\n';

   std::vector<int> remainIndex;
   for (size_t i = 0; i < vX.size(); i++)
      remainIndex.push_back(i);

   TVector3 V1, V2;
   std::vector<int> inliners;
   inliners.clear();

   double sigma = fRANSACThreshold / 1.96;
   double dataSigma2 = sigma * sigma;
   std::vector<std::pair<double, int>> IdxMod1;
   std::vector<std::pair<double, int>> IdxMod2;

   for (int i = 0; i < fRANSACMaxIteration; i++) {

      if (remainIndex.size() < fRANSACMinPoints)
         break;

      std::vector<int> Rsamples = RandSam(remainIndex, fRandSamplMode);
      EstimModel(Rsamples);

      // Calculate squared errors
      double minError = 1e5, maxError = -1e5;
      for (int &j : remainIndex) {
         double error = EstimError(j);
         if (error < minError)
            minError = error;
         if (error > maxError)
            maxError = error;
      }

      // Estimate the inlier ratio using EM
      const double nu = maxError - minError;
      double gamma = 0.5;
      for (int iter = 0; iter < 5; iter++) {
         double sumPosteriorProb = 0;
         const double probOutlier = (1 - gamma) / nu;
         const double probInlierCoeff = gamma / sqrt(2 * TMath::Pi() * dataSigma2);

         for (int &j : remainIndex) {
            double error = EstimError(j);
            double probInlier = probInlierCoeff * exp(-0.5 * error * error / dataSigma2);
            sumPosteriorProb += probInlier / (probInlier + probOutlier);
         }
         gamma = sumPosteriorProb / remainIndex.size();
         // gamma = sumPosteriorProb / totalNbSamples;
      }

      double sumLogLikelihood = 0;
      int nbInliers = 0;
      // std::vector<int> inlIdxR;

      // Evaluate the model
      const double probOutlier = (1 - gamma) / nu;
      const double probInlierCoeff = gamma / sqrt(2 * TMath::Pi() * dataSigma2);
      for (int &j : remainIndex) {
         double error = EstimError(j);
         double probInlier = probInlierCoeff * exp(-0.5 * error * error / dataSigma2);
         // if((probInlier + probOutlier)>0) sumLogLikelihood = sumLogLikelihood - log(probInlier + probOutlier);

         if (error * error < dataSigma2) {
            if ((probInlier + probOutlier) > 0)
               sumLogLikelihood = sumLogLikelihood - log(probInlier + probOutlier);
            nbInliers++;
            // inlIdxR.push_back(*j);
         }
      }

      // std::cout <<sumLogLikelihood<< "/* message */" << '\n';

      if (nbInliers > fRANSACMinPoints) {
         double scale = sumLogLikelihood / nbInliers;
         // std::cout <<sumLogLikelihood<< " Likelihood  " << '\n';
         if (sumLogLikelihood < 0 || std::isinf(sumLogLikelihood))
            scale = 0;
         // std::cout <<sumLogLikelihood<< " Likelihood  "<< scale<< << '\n';
         IdxMod1.emplace_back(scale, Rsamples[0]);
         IdxMod2.emplace_back(scale, Rsamples[1]);
      }

   } // for Mlesac interactions

   // sort clusters
   sort(IdxMod1.begin(), IdxMod1.end());
   sort(IdxMod2.begin(), IdxMod2.end());

   remainIndex.clear(); // track remaining points
   for (size_t i = 0; i < vX.size(); i++)
      remainIndex.push_back(i);

   // extract inliers using the models
   for (int i = 0; i < IdxMod1.size(); ++i) {
      std::vector<int> ModInx = {IdxMod1[i].second, IdxMod2[i].second};
      EstimModel(ModInx);
      std::vector<int> inlIdxR;
      ModInx.clear();

      if (remainIndex.size() < fRANSACMinPoints)
         break;

      int counter = 0;

      for (int &j : remainIndex) {
         double error = EstimError(j);

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
