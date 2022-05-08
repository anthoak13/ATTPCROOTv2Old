#include "AtIndependentSample.h"

#include "AtHit.h"

#include <TRandom3.h>
using namespace AtTools;
/**
 * Assumes fCDF is already setup and we are just using it.
 */
std::vector<AtHit> AtIndependentSample::SampleHits(int N)
{
   // Using the sampled indices, return a vector of positions
   std::vector<AtHit> ret;
   for (auto ind : sampleIndicesFromCDF(N))
      ret.push_back(fHits->at(ind));
   return ret;
}

void AtIndependentSample::SetHitsToSample(const std::vector<AtHit> *hits)
{
   fHits = hits;
   FillCDF();
}
