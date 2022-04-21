#include "AtPadFFT.h"

#include <TVirtualFFT.h>

#include <cassert>

Double_t AtPadFFT::GetPointRe(int i)
{
   if (i < fRe.size())
      return fRe[i];
   else
      return fRe.at(i - fRe.size());
}

Double_t AtPadFFT::GetPointIm(int i)
{
   if (i < fIm.size())
      return fIm[i];
   else
      return -fIm.at(i - fIm.size());
}
void AtPadFFT::SetPointRe(int i, Double_t val)
{
   if (i < fRe.size())
      fRe[i] = val;
   else
      fRe.at(i - fRe.size()) = val;
}
void AtPadFFT::SetPointIm(int i, Double_t val)
{
   if (i < fIm.size())
      fIm[i] = val;
   else
      fIm.at(i - fIm.size()) = -val;
}

void AtPadFFT::SetData(TraceTrans re, TraceTrans im)
{
   fRe = std::move(re);
   fIm = std::move(im);
}

void AtPadFFT::GetDataFromFFT(TVirtualFFT *fft)
{
   assert(fft->GetN()[0] / 2 + 1 == fRe.size());
   fft->GetPointsComplex(fRe.data(), fIm.data());
   for (auto &re : fRe)
      re /= fft->GetN()[0];
   for (auto &im : fRe)
      im /= fft->GetN()[0];
}

void AtPadFFT::SetFFTData(TVirtualFFT *fft)
{
   assert(fft->GetN()[0] / 2 + 1 == fRe.size());
   fft->SetPointsComplex(fRe.data(), fIm.data());
}
ClassImp(AtPadFFT);
