#include "AtWeightedGaussian.h"

#include "AtHit.h"

#include <Math/PdfFuncMathCore.h>
#include <TRandom3.h>
using namespace AtTools;
void AtWeightedGaussian::SetHitsToSample(const std::vector<AtHit> *hits)
{
   AtSampleFromReference::SetHitsToSample(hits);
   fChargeSample.SetHitsToSample(hits);
}

std::vector<double> AtWeightedGaussian::PDF(const AtHit &hit)
{
   auto dist = (fReferenceHit.GetPosition() - hit.GetPosition()).Mag2();
   dist = std::sqrt(dist);
   return {ROOT::Math::gaussian_pdf(dist, fSigma), hit.GetCharge()};
}

void AtWeightedGaussian::SampleReferenceHit()
{
   AtChargeWeighted charge;
   charge.SetHitsToSample(fHits);
   SetReferenceHit(std::move(charge.SampleHits(1)[0]));
}
