#include "AtUniform.h"

#include "AtSample.h" // for RandomSample

class AtHit;

using namespace RandomSample;

std::vector<double> AtUniform::PDF(const AtHit &hit)
{
   return {1};
}
