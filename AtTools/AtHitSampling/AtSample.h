#ifndef ATHITSAMPLER_H
#define ATHITSAMPLER_H

#include <Math/Point3D.h>

#include <algorithm>
#include <vector>
class AtHit;

namespace AtTools {

/**
 * Interface for randomly sampling a collection of AtHits according
 * to the cumulitive distribution function fCDF.
 *
 * @ingroup AtHitSampling
 */
class AtSample {
protected:
   const std::vector<AtHit> *fHits; //< Hits to sample from
   std::vector<double> fCDF;        //< Cummulative distribution function for hits
   bool fWithReplacement{false};    //< If we should sample with replacement

public:
   enum class SampleMethod;

   virtual std::vector<AtHit> SampleHits(int N);
   std::vector<ROOT::Math::XYZPoint> SamplePoints(int N);

   virtual void SetHitsToSample(const std::vector<AtHit> *hits) = 0;

   void SetSampleWithReplacement(bool val) { fWithReplacement = val; }

   static std::unique_ptr<AtSample> CreateSampler(SampleMethod method);

protected:
   /**
    * Computes the unnormalized marginal PDFs at the hit.
    *
    * For example, a charge-weighted and spacial gaussian would looke like
    * `return {charge, gaussian_pdf(distance, sigma)};`
    *
    * @return vector where each element is a different maginal pdf
    */
   virtual std::vector<double> PDF(const AtHit &hit) = 0;
   void FillCDF();

   std::vector<int> sampleIndicesFromCDF(int N, std::vector<int> vetoed = {});
   int getIndexFromCDF(double r, double rmCFD, std::vector<int> vetoed);
   template <typename T>
   static inline bool isInVector(T val, std::vector<T> vec)
   {
      if (vec.size() == 0)
         return false;
      return std::find(vec.begin(), vec.end(), val) != vec.end();
   }
};
} // namespace AtTools

#endif //#ifndef ATHITSAMPLER_H
