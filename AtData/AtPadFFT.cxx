#include "AtPadFFT.h"

#include <TVirtualFFT.h>

#include <cassert>
std::unique_ptr<AtPad> AtPadFFT::Clone()
{
   return std::make_unique<AtPadFFT>(*this);
}

Double_t AtPadFFT::GetPointRe(int i)
{

   if (i < fRe.size())
      return fRe[i];
   else
      return fRe.at(fAdc.size() - i);
}

Double_t AtPadFFT::GetPointIm(int i)
{
   if (i < fIm.size())
      return fIm[i];
   else
      return -fIm.at(fAdc.size() - i);
}
void AtPadFFT::SetPointRe(int i, Double_t val)
{
   if (i < fRe.size())
      fRe[i] = val;
   else
      fRe.at(fAdc.size() - i) = val;
}
void AtPadFFT::SetPointIm(int i, Double_t val)
{
   if (i < fIm.size())
      fIm[i] = val;
   else
      fIm.at(fAdc.size() - i) = -val;
}

void AtPadFFT::SetData(TraceTrans re, TraceTrans im)
{
   fRe = std::move(re);
   fIm = std::move(im);
}

void AtPadFFT::GetDataFromFFT(TVirtualFFT *fft)
{
   assert(fft->GetN()[0] / 2 + 1 == fRe.size());
   for (int i = 0; i < fRe.size(); ++i)
      fft->GetPointComplex(i, fRe.at(i), fIm.at(i));
   for (auto &re : fRe)
      re /= fft->GetN()[0];
   for (auto &im : fIm)
      im /= fft->GetN()[0];
}

void AtPadFFT::SetFFTData(TVirtualFFT *fft)
{
   assert(fft->GetN()[0] / 2 + 1 == fRe.size());
   fft->SetPointsComplex(fRe.data(), fIm.data());
}
ClassImp(AtPadFFT);
