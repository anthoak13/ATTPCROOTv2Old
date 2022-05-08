#ifndef ATSAMPLEFROMREFERENCE_H
#define ATSAMPLEFROMREFERENCE_H

#include "AtHit.h"
#include "AtSample.h"
namespace AtTools {

/**
 * Interface for sampling a collection of AtHits where the PDF depends on
 * a reference hit.
 *
 * @ingroup AtHitSampling
 */
class AtSampleFromReference : public AtSample {
protected:
   AtHit fReferenceHit; //< Hit to use to construct the CDF/PDF

public:
   virtual std::vector<AtHit> SampleHits(int N) override;
   virtual void SetHitsToSample(const std::vector<AtHit> *hits) override;
   void SetReferenceHit(AtHit hit) { fReferenceHit = std::move(hit); }
   const AtHit &GetReferenceHit() const { return fReferenceHit; }

protected:
   virtual void SampleReferenceHit();
};
} // namespace AtTools
#endif //#ifndef ATSAMPLEFROMREFERENCE_H
