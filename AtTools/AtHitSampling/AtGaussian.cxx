#include "AtGaussian.h"

#include "AtHit.h"

#include <Math/PdfFuncMathCore.h>
#include <TRandom3.h>
using namespace AtTools;
std::vector<double> AtGaussian::PDF(const AtHit &hit)
{
   auto dist = (fReferenceHit.GetPosition() - hit.GetPosition()).Mag2();
   dist = std::sqrt(dist);
   return {ROOT::Math::gaussian_pdf(dist, fSigma)};
}
