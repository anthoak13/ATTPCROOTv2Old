#include "AtChargeWeighted.h"

#include "AtHit.h"
using namespace AtTools;
std::vector<double> AtChargeWeighted::PDF(const AtHit &hit)
{
   return {hit.GetCharge()};
}
