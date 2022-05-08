#ifndef ATUNIFORM_H
#define ATUNIFORM_H

#include "AtIndependentSample.h"

namespace AtTools {

/**
 * Class for uniformly sampling a collection of AtHits
 *
 * @ingroup AtHitSampling
 */
class AtUniform : public AtIndependentSample {

protected:
   virtual std::vector<double> PDF(const AtHit &hit) override;
};
} // namespace AtTools
#endif //#ifndef ATSAMPLEUNIFORM_H
