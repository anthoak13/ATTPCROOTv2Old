#ifndef ATCHARGEWEIGHTED_H
#define ATCHARGEWEIGHTED_H

#include "AtIndependentSample.h"
namespace AtTools {

/**
 * Class for sampling a collection of AtHits according to charge
 *
 * Follows the PDF: P(q) = q/TotalCharge
 * @ingroup AtHitSampling
 */
class AtChargeWeighted : public AtIndependentSample {

protected:
   virtual std::vector<double> PDF(const AtHit &hit) override;
};
} // namespace AtTools
#endif // ATCHARGEWEIGHTED
