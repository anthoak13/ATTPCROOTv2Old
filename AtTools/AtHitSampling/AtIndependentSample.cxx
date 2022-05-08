#include "AtIndependentSample.h"

#include "AtHit.h"

#include <TRandom3.h>
using namespace AtTools;

void AtIndependentSample::SetHitsToSample(const std::vector<AtHit> *hits)
{
   fHits = hits;
   FillCDF();
}
