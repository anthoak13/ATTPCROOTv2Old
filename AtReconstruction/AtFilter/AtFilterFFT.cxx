#include "AtFilterFFT.h"

#include "AtPad.h"
#include "AtPadFFT.h"
#include "AtRawEvent.h"

#include <FairLogger.h>
#include <FairRootManager.h>

#include <Rtypes.h>
#include <TCanvas.h>
#include <TComplex.h>
#include <TH1.h>
#include <TVirtualFFT.h>

bool AtFilterFFT::AddFreqRange(AtFreqRange range)
{
   auto canAdd = isValidFreqRange(range);
   if (canAdd) {
      fFreqRanges.push_back(std::move(range));
      auto dFact = (range.fEndFact - range.fBeginFact) / (range.fEndFreq - range.fBeginFreq);
      for (int i = range.fBeginFreq; i <= range.fEndFreq; ++i)
         fFactors[i] = range.fBeginFact + dFact * (i - range.fBeginFreq);
   }
   return canAdd;
}

void AtFilterFFT::Init()
{
   std::vector<Int_t> dimSize = {512};

   // Create a FFT object that we own ("K"), that will optimize the transform ("M"),
   // and is a forward transform from real data to complex ("R2C")
   fFFT = std::unique_ptr<TVirtualFFT>(TVirtualFFT::FFT(1, dimSize.data(), "R2C M K"));

   // Create a FFT object that we own ("K"), that will optimize the transform ("M"),
   // and is a backwards transform from complex to Reak ("C2R")
   fFFTbackward = std::unique_ptr<TVirtualFFT>(TVirtualFFT::FFT(1, dimSize.data(), "C2R M K"));

   if(fSaveTransform) // Hide this in an if so it can run outside of a FairRoot instance
      FairRootManager::Instance()->Register("AtEventFFT", "AtTPC", &fTransformArray, true);
}

AtRawEvent *AtFilterFFT::ConstructOutputEvent(TClonesArray *outputEventArray, AtRawEvent *inputEvent)
{
   // Just do a copy if we're not intrested in saving the transformed output
   if (!fSaveCutTransform)
      return new ((*outputEventArray)[0]) AtRawEvent(*inputEvent);

   // Make a copy of the event changing the pad type to AtPadFFT (for trasnformation information)
   // and return that
   fFilteredEvent = dynamic_cast<AtRawEvent *>(new ((*outputEventArray)[0]) AtRawEvent());
   fFilteredEvent->CopyAllButData(inputEvent);
   for (auto &pad : inputEvent->GetPads())
      fFilteredEvent->AddPad(std::unique_ptr<AtPad>(new AtPadFFT(*pad)));
   return fFilteredEvent;
}

void AtFilterFFT::InitEvent(AtRawEvent *event)
{
   fTransformArray.Clear();
   fTransformedEvent = dynamic_cast<AtRawEvent *>(fTransformArray.ConstructedAt(0));
   fTransformedEvent->CopyAllButData(event);
}

void AtFilterFFT::Filter(AtPad *pad)
{
   // Get data and transform
   fFFT->SetPoints(pad->GetADC().data());
   fFFT->Transform();
   applyFrequencyCutsAndSetInverseFFT();
   fFFTbackward->Transform();

   if (fSaveTransform) {
      // Copy the pad to the unfiltered event before updating pad's adc values
      auto transformedPad =
         dynamic_cast<AtPadFFT *>(fTransformedEvent->AddPad(std::unique_ptr<AtPad>(new AtPadFFT(*pad))));
      // Add the result of the transformation
      transformedPad->GetDataFromFFT(fFFT.get());
   }

   if (fSaveCutTransform)
      dynamic_cast<AtPadFFT *>(pad)->GetDataFromFFT(fFFTbackward.get());

   for (int i = 0; i < pad->GetADC().size(); ++i)
      pad->SetADC(i, fFFTbackward->GetPointReal(i));
}

bool AtFilterFFT::isValidFreqRange(const AtFreqRange &range)
{
   bool isValid = true;
   isValid &= range.fBeginFreq <= range.fEndFreq;
   isValid &= range.fBeginFact >= 0 && range.fBeginFact <= 1;
   isValid &= range.fEndFact >= 0 && range.fEndFact <= 1;
   isValid &= !doesFreqRangeOverlap(range);
   return isValid;
}

bool AtFilterFFT::doesFreqRangeOverlap(const AtFreqRange &newRange)
{
   for (auto &range : fFreqRanges) {
      // Check to make sure the minimum frequency does not exist within the bounds
      // of any added frequency, except for the max frequency when the factors are
      // the same
      if (newRange.fBeginFreq >= range.fBeginFreq && newRange.fBeginFreq < range.fEndFreq)
         return true;
      if (newRange.fBeginFreq == range.fEndFreq && newRange.fBeginFact != range.fEndFact)
         return true;

      // Check to make sure the maximum frequency does not exist within the bounds of any added
      // frequency, except for when it matches the minimum frequency with the same factor
      if (newRange.fEndFreq > range.fBeginFreq && newRange.fEndFreq <= range.fEndFreq)
         return true;
      if (newRange.fEndFreq == range.fBeginFreq && newRange.fEndFact != range.fBeginFreq)
         return true;
   }

   return false;
}

void AtFilterFFT::applyFrequencyCutsAndSetInverseFFT()
{
   for (int i = 0; i < fFFT->GetN()[0]; ++i) {
      Double_t re, im;
      fFFT->GetPointComplex(i, re, im);
      if (fFactors.find(i) != fFactors.end()) {
         re *= fFactors[i];
         im *= fFactors[i];
      }
      fFFTbackward->SetPoint(i, re, im);
   }
}

bool operator<(const AtFilterFFT::AtFreqRange &lhs, const AtFilterFFT::AtFreqRange &rhs)
{
   return lhs.fBeginFreq < rhs.fBeginFreq;
}

void AtFilterFFT::DumpFactors()
{
   for (const auto &pair : fFactors)
      std::cout << pair.first << " " << pair.second << std::endl;
}
