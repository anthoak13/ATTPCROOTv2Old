#include "AtSampleFromReference.h"

#include "AtHit.h"

#include <TRandom3.h>
using namespace AtTools;
std::vector<AtHit> AtSampleFromReference::SampleHits(int N)
{
   SampleReferenceHit();
   FillCDF();

   // Using the sampled indices, return a vector of positions
   std::vector<AtHit> ret;
   for (auto ind : sampleIndicesFromCDF(N))
      ret.push_back(fHits->at(ind));
   return ret;
}

void AtSampleFromReference::SampleReferenceHit()
{
   int refIndex = gRandom->Uniform() * fHits->size();
   SetReferenceHit(fHits->at(refIndex));
}

void AtSampleFromReference::SetHitsToSample(const std::vector<AtHit> *hits)
{
   fHits = hits;
}
