#ifndef ATUNIFORM_H
#define ATUNIFORM_H

#include "AtIndependentSample.h"

#include <vector> // for vector
class AtHit;

namespace RandomSample {

/**
 * @brief Uniformly sample a collection of AtHits
 *
 * @ingroup AtHitSampling
 */
class AtUniform : public AtIndependentSample {

protected:
   virtual std::vector<double> PDF(const AtHit &hit) override;
};
} // namespace RandomSample
#endif //#ifndef ATSAMPLEUNIFORM_H
