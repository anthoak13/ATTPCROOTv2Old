#include "AtEstimatorMethods.h"

#include "AtPattern.h"
using namespace SampleConsensus;

int EvaluateRansac(AtPatterns::AtPattern *model, const std::vector<AtHit> &hitArray, double distanceThreshold)
{
   int nbInliers = 0;
   double weight = 0;

   for (const auto &hit : hitArray) {
      auto &pos = hit.GetPosition();
      double error = model->DistanceToPattern(pos);
      error = error * error;
      if (error < (distanceThreshold * distanceThreshold)) {
         nbInliers++;
         weight += error;
      }
   }
   model->SetChi2(weight / nbInliers);
   return nbInliers;
}

int EvaluateMlesac(AtPatterns::AtPattern *model, const std::vector<AtHit> &hitArray, double distanceThreshold)
{
   double sigma = distanceThreshold / 1.96;
   double dataSigma2 = sigma * sigma;

   // Calculate min and max errors
   double minError = 1e5, maxError = -1e5;
   for (const auto &hit : hitArray) {
      double error = model->DistanceToPattern(hit.GetPosition());
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
      const double probInlierCoeff = gamma / sqrt(2 * M_PI * dataSigma2);

      for (const auto &hit : hitArray) {
         double error = model->DistanceToPattern(hit.GetPosition());
         double probInlier = probInlierCoeff * exp(-0.5 * error * error / dataSigma2);
         sumPosteriorProb += probInlier / (probInlier + probOutlier);
      }
      gamma = sumPosteriorProb / hitArray.size();
   }

   double sumLogLikelihood = 0;
   int nbInliers = 0;

   // Evaluate the model
   const double probOutlier = (1 - gamma) / nu;
   const double probInlierCoeff = gamma / sqrt(2 * M_PI * dataSigma2);
   for (const auto &hit : hitArray) {
      double error = model->DistanceToPattern(hit.GetPosition());
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

double GetMedian(std::vector<double> &vec)
{
   if (vec.empty())
      return 0;

   auto n = vec.size() / 2;

   std::nth_element(vec.begin(), vec.begin() + n, vec.end());
   auto med = vec[n];
   if (vec.size() % 2 == 0) { // if even, get average of middle two entries
      med = (*std::max_element(vec.begin(), vec.begin() + n) + med) / 2.0;
   }
   return med;
}

int EvaluateLmeds(AtPatterns::AtPattern *model, const std::vector<AtHit> &hitArray, double distanceThreshold)
{
   std::vector<double> errorsVec;
   // Loop through point and if it is an inlier, then add the error**2 to weight
   for (const auto &hit : hitArray) {

      double error = model->DistanceToPattern(hit.GetPosition());
      error = error * error;
      if (error < (distanceThreshold * distanceThreshold))
         errorsVec.push_back(error);
   }
   model->SetChi2(GetMedian(errorsVec) / errorsVec.size());
   return errorsVec.size();
}
