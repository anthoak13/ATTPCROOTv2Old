#include "AtChargeWeighted.h"

#include "AtHit.h"
using namespace RandomSample;
std::vector<double> AtChargeWeighted::PDF(const AtHit &hit)
{
   return {hit.GetCharge()};
}
