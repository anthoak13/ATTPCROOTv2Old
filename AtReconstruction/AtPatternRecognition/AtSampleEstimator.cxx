#include "AtSampleEstimator.h"

#include "AtPattern.h"

#include <cmath>

/**
 * Implementation of RANSAC estimator
 */
class AtEstimatorRansac {
public:
   static int EvaluateModel(AtPattern *model, const std::vector<AtHit> &hitArray, double distanceThreshold)
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
};

/**
 * Implementation of MLESAC estimator
 */
class AtEstimatorMlesac {
public:
   static int EvaluateModel(AtPattern *model, const std::vector<AtHit> &hitArray, double distanceThreshold)
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
};

/**
 * Implementation of LMedS estimator
 */
class AtEstimatorLmeds {
public:
   static int EvaluateModel(AtPattern *model, const std::vector<AtHit> &hitArray, double distanceThreshold)
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

   static double GetMedian(std::vector<double> &vec)
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
};

int AtSampleEstimator::EvaluateModel(AtPattern *model, const std::vector<AtHit> &hits, double distThresh,
                                     Estimator estimator)
{
   switch (estimator) {
   case (Estimator::kRANSAC): return AtEstimatorRansac::EvaluateModel(model, hits, distThresh);
   case (Estimator::kLMedS): return AtEstimatorLmeds::EvaluateModel(model, hits, distThresh);
   case (Estimator::kMLESAC): return AtEstimatorMlesac::EvaluateModel(model, hits, distThresh);
   default: return 0;
   }
}
