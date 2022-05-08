#ifndef ATINDEPENDENTSAMPLE_H
#define ATINDEPENDENTSAMPLE_H
#include "AtSample.h"
namespace AtTools {

/**
 * Interface for sampling a collection of AtHits where the points sampled are independent of one another.
 * @ingroup AtHitSampling
 */

class AtIndependentSample : public AtSample {
public:
   virtual void SetHitsToSample(const std::vector<AtHit> *hits) override;

protected:
   virtual std::vector<double> PDF(const AtHit &hit) override = 0;
};
} // namespace AtTools

#endif //#ifndef ATINDEPENDENTSAMPLE_H
