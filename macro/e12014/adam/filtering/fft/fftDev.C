#include "TCanvas.h"
#include "TFile.h"
#include "TVirtualFFT.h"

#ifndef __CLING__
#include "../build/include/AtPadFFT.h"
#include "../build/include/AtRawEvent.h"

#endif

#include "../../helper.h"

/***** "Public" functions ******/
void loadRun(int runNum);
void runFFT(int padNum);

/***** "Private" variables *****/

TFile *oFile = nullptr;
TCanvas *cFFT = nullptr;

void fftDev(int runNum = 195)
{
   loadRun(runNum);
}
void loadRun(int runNum)
{
   TString filePath = "/mnt/analysis/e12014/TPC/unpacked/run_%04d.root";
   loadRun(TString::Format(filePath, runNum), "AtRawEvent", "AtRawEventSubtracted", "AtEventH");

   if (dynamic_cast<AtPadFFT *>(rawEventPtr->GetPads().back().get()) == nullptr)
      LOG(error) << "Raw event branch does not contain FFT pads!";
   if (dynamic_cast<AtPadFFT *>(rawEventFilteredPtr->GetPads().back().get()) == nullptr)
      LOG(error) << "Filtered raw event branch does not contain FFT pads!";

   if (oFile != nullptr)
      delete oFile;
   oFile = new TFile(TString::Format("output/traces-%d.root", runNum), "RECREATE");
}

void viewFFT(int eventNum, int padNum)
{
   if (!loadEvent(eventNum))
      return;
   if (!loadPad(padNum))
      return;

   if (cFFT == nullptr) {
      cFFT = new TCanvas("cFFT", "FFT", 900, 900);
      cFFT->Divide(2, 2);
   }

   auto pad = dynamic_cast<AtPadFFT *>(rawEventPtr->GetPad(padNum));
   auto padFiltered = dynamic_cast<AtPadFFT *>(rawEventFilteredPtr->GetPad(padNum));

   TH1F *hFilteredTrace = new TH1F("filteredTrace", "Filtered Trace", 512, 0, 511);
   TH1F *hFilteredMag = new TH1F("filtereMag", "Filtered Magnitude", 512 / 2 + 1, 0, 512 / 2);
   TH1F *hMag = new TH1F("mag", "Magnitude", 512 / 2 + 1, 0, 512 / 2);

   for (int i = 0; i < 512 / 2 + 1; ++i) {
      auto point = pad->GetPoint(i);
      auto mag = std::sqrt(point.first * point.first + point.second * point.second);
      hMag->SetBinContent(i + 1, mag);

      point = padFiltered->GetPoint(i);
      mag = std::sqrt(point.first * point.first + point.second * point.second);
      hFilteredMag->SetBinContent(i + 1, mag);
   }
   for (int i = 0; i < 512; ++i)
      hFilteredTrace->SetBinContent(i + 1, padFiltered->GetADC(i));

   cFFT->cd(1);
   hTrace->Draw();
   cFFT->cd(2);
   hFilteredTrace->Draw();
   cFFT->cd(3);
   hMag->Draw();
   cFFT->cd(4);
   hFilteredMag->Draw();
}
