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
int AtMlesacMod::evaluateModel(AtTrackModel *model, const std::vector<int> &pointsToCheck,
                               const std::vector<AtHit> &hitArray)
{
   double sigma = fRANSACThreshold / 1.96;
   double dataSigma2 = sigma * sigma;

   // Calculate min and max errors
   double minError = 1e5, maxError = -1e5;
   for (int j : pointsToCheck) {
      double error = model->DistanceToModel(hitArray.at(j).GetPosition());
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

      for (int j : pointsToCheck) {
         double error = model->DistanceToModel(hitArray.at(j).GetPosition());
         double probInlier = probInlierCoeff * exp(-0.5 * error * error / dataSigma2);
         sumPosteriorProb += probInlier / (probInlier + probOutlier);
      }
      gamma = sumPosteriorProb / pointsToCheck.size();
   }

   double sumLogLikelihood = 0;
   int nbInliers = 0;

   // Evaluate the model
   const double probOutlier = (1 - gamma) / nu;
   const double probInlierCoeff = gamma / sqrt(2 * TMath::Pi() * dataSigma2);
   for (int j : pointsToCheck) {
      double error = model->DistanceToModel(hitArray.at(j).GetPosition());
      double probInlier = probInlierCoeff * exp(-0.5 * error * error / dataSigma2);
      // if((probInlier + probOutlier)>0) sumLogLikelihood = sumLogLikelihood - log(probInlier + probOutlier);

      if (error * error < dataSigma2) {
         if ((probInlier + probOutlier) > 0)
            sumLogLikelihood = sumLogLikelihood - log(probInlier + probOutlier);
         nbInliers++;
      }
   }
   double scale = sumLogLikelihood / nbInliers;
   // std::cout <<sumLogLikelihood<< " Likelihood  " << '\n';
   if (sumLogLikelihood < 0 || std::isinf(sumLogLikelihood))
      scale = 0;

   model->SetChi2(scale);
   return nbInliers;
}
