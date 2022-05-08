#ifndef ATGAUSSIAN_H
#define ATGAUSSIAN_H
#include "AtSampleFromReference.h"
namespace AtTools {

/**
 * Class for sampling a collection of AtHits in a spacial gaussian distribution around
 * the reference hit.
 *
 * @ingroup AtHitSampling
 */
class AtGaussian : public AtSampleFromReference {
protected:
   double fSigma; //< Sigma of gaussian around fReferencehit to sample [mm]

public:
   AtGaussian(double sigma) : fSigma(sigma) {}

protected:
   virtual std::vector<double> PDF(const AtHit &hit) override;
};
} // namespace AtTools
#endif // ATSAMPLEGAUSSIAN_H
