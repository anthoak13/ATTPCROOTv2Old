#include "AtSample.h"

#include "AtHit.h"

#include <TRandom3.h>

#include <numeric>
using namespace AtTools;

/**
 * Sample spacial locations (XYZPoints) from fHits.
 */
std::vector<ROOT::Math::XYZPoint> AtSample::SamplePoints(int N)
{
   std::vector<ROOT::Math::XYZPoint> ret;
   for (const auto &hit : SampleHits(N))
      ret.push_back(hit.GetPosition());
   return ret;
}

/**
 * Get a list of indices aampled according to the constructed fCDF (assumes FillCDF() has been called already)
 *
 * @param[in] N number of points to sample. If sampling without replacement this should be small compared to the total
 * number of hits in the hit array.
 * @param[in] vetoed Indices to not sample (even if sampling with replacement).
 */
std::vector<int> AtSample::sampleIndicesFromCDF(int N, const std::vector<int> &vetoed)
{
   std::vector<int> sampledInd;
   while (sampledInd.size() < N) {
      auto r = gRandom->Uniform();
      // Get the index i where CDF[i] >= r and CDF[i-1] < r
      int hitInd = std::distance(fCDF.begin(), std::lower_bound(fCDF.begin(), fCDF.end(), r));
      // Save the index if it has not already been sampled
      auto canSample = fWithReplacement || !isInVector(hitInd, sampledInd);
      if (canSample && !isInVector(hitInd, vetoed))
         sampledInd.push_back(hitInd);
   }

   return sampledInd;
}

/**
 * Fill the cumulitive distribution function to sample using the marginal PDFs returned by the
 * function PDF(const AtHit &hit) from every entry in the vector fHits.
 *
 * Each marginal PDF is normalized to construct the final CDF. Assumes the marginal PDFs are independent
 * (i.e. the joint PDF is just the product of all marginal PDFs)
 */
void AtSample::FillCDF()
{
   std::vector<double> normalization;
   fCDF.clear();
   for (const auto &hit : *fHits) {

      // Get the unnormalized marginal and joint PDFs
      auto pdfMarginal = PDF(hit);
      auto pdfJoint = std::accumulate(pdfMarginal.begin(), pdfMarginal.end(), 1.0,
                                      std::multiplies<>()); // Has to be 1.0 not 1 or return type is deduced as int

      // If this is the first hit, setup the normalization vector
      if (normalization.size() == 0)
         normalization.assign(pdfMarginal.size(), 0);

      for (int i = 0; i < pdfMarginal.size(); ++i)
         normalization[i] += pdfMarginal[i];

      if (fCDF.size() == 0)
         fCDF.push_back(pdfJoint);
      else
         fCDF.push_back(pdfJoint + fCDF.back());
   }

   // Get the total normalization from the marginal PDFs, and normalize the CDF
   auto norm = std::accumulate(normalization.begin(), normalization.end(), 1.0, std::multiplies<>());
   for (auto &elem : fCDF) {
      elem /= norm;
   }
}
