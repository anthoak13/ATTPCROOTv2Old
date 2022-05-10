#include "AtIndependentSample.h"

#include "AtHit.h"

#include <TRandom3.h>
using namespace RandomSample;

void AtIndependentSample::SetHitsToSample(const std::vector<AtHit> *hits)
{
   fHits = hits;
   FillCDF();
}
