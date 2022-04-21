#include "AtFFTFilter.h"

#include "AtPad.h"
#include "AtPadFFT.h"
#include "AtRawEvent.h"

#include <FairLogger.h>
#include <FairRootManager.h>

#include <Rtypes.h>
#include <TH1.h>
#include <TVirtualFFT.h>

bool AtFFTFilter::AddFreqRange(AtFreqRange range)
{
   auto canAdd = isValidFreqRange(range);
   if (canAdd) {
      fFreqRanges.push_back(std::move(range));
      auto dFact = (range.fEndFact - range.fBeginFact) / (range.fEndFreq - range.fBeginFreq - 1);
      for (int i = range.fBeginFreq; i < range.fEndFreq; ++i)
         fFactors[i] = range.fBeginFact + dFact * (i - range.fBeginFreq);
   }
   return canAdd;
}

void AtFFTFilter::Init()
{
   std::vector<Int_t> dimSize = {512};

   // Create a FFT object that we own ("K"), that will optimize the transform ("M"),
   // and is a forward transform from real data to complex ("R2C")
   fFFT = std::unique_ptr<TVirtualFFT>(TVirtualFFT::FFT(1, dimSize.data(), "R2C M K"));

   // Create a FFT object that we own ("K"), that will optimize the transform ("M"),
   // and is a backwards transform from complex to Reak ("C2R")
   fFFTbackward = std::unique_ptr<TVirtualFFT>(TVirtualFFT::FFT(1, dimSize.data(), "C2R M K"));

   FairRootManager::Instance()->Register("AtEventFFT", "AtTPC", &fTransformArray, fSaveTransform);
}

AtRawEvent *AtFFTFilter::ConstructOutputEvent(TClonesArray *outputEventArray, AtRawEvent *inputEvent)
{
   // Just do a copy if we're not intrested in saving the transformed output
   if (!fSaveCutTransform)
      return new ((*outputEventArray)[0]) AtRawEvent(*inputEvent);

   // Make a copy of the event changing the pad type to AtPadFFT (for trasnformation information)
   auto coppiedEvt = dynamic_cast<AtRawEvent *>(new ((*outputEventArray)[0]) AtRawEvent());
   coppiedEvt->CopyAllButData(inputEvent);
   for (auto &pad : inputEvent->GetPads())
      coppiedEvt->AddPad(std::unique_ptr<AtPad>(new AtPadFFT(*pad)));
   std::cout << "Copied raw event to " << coppiedEvt << std::endl;

   return coppiedEvt;
}

void AtFFTFilter::InitEvent(AtRawEvent *event)
{
   fTransformArray.Clear();
   auto transformedEvent = dynamic_cast<AtRawEvent *>(fTransformArray.ConstructedAt(0));
   transformedEvent->CopyAllButData(event);
}

void AtFFTFilter::Filter(AtPad *pad)
{
   // Get the event we will use to store the transformation as we go
   auto transformedEvent = dynamic_cast<AtRawEvent *>(fTransformArray.ConstructedAt(0));

   // Add the pad to the transformed event, and then fill it with the transformed data
   auto transfPad = dynamic_cast<AtPadFFT *>(transformedEvent->AddPad(std::unique_ptr<AtPad>(new AtPadFFT(*pad))));
   fFFT->SetPoints(transfPad->GetADC().data());
   fFFT->Transform();
   transfPad->GetDataFromFFT(fFFT.get());

   auto filteredPad = dynamic_cast<AtPadFFT *>(pad);

   // Apply the frequency cuts and setup the FFT with the filtered data
   filteredPad->GetDataFromFFT(fFFT.get());
   applyFrequencyCuts(filteredPad);
   filteredPad->SetFFTData(fFFTbackward.get());

   // Transform back to time-space and update the trace
   fFFTbackward->Transform();
   for (int i = 0; i < filteredPad->GetADC().size(); ++i)
      filteredPad->SetADC(i, fFFTbackward->GetPointReal(i));
}

bool AtFFTFilter::isValidFreqRange(const AtFreqRange &range)
{
   bool isValid = true;
   isValid &= range.fBeginFreq <= range.fEndFreq;
   isValid &= range.fBeginFact >= 0 && range.fBeginFact <= 1;
   isValid &= range.fEndFact >= 0 && range.fEndFact <= 1;
   isValid &= !doesFreqRangeOverlap(range);
   return isValid;
}

bool AtFFTFilter::doesFreqRangeOverlap(const AtFreqRange &newRange)
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

void AtFFTFilter::applyFrequencyCuts(AtPadFFT *pad)
{
   for (const auto &pair : fFactors) {
      auto i = pair.first;
      pad->SetPointIm(i, pad->GetPointIm(i) * pair.second);
      pad->SetPointRe(i, pad->GetPointRe(i) * pair.second);
   }
}

bool operator<(const AtFFTFilter::AtFreqRange &lhs, const AtFFTFilter::AtFreqRange &rhs)
{
   return lhs.fBeginFreq < rhs.fBeginFreq;
}

void AtFFTFilter::DumpFactors()
{
   for (const auto &pair : fFactors)
      std::cout << pair.first << " " << pair.second << std::endl;
}
