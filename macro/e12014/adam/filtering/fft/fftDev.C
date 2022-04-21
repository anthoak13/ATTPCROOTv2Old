#include "AtFFTFilter.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TVirtualFFT.h"

#include "../../helper.h"

/***** "Public" functions ******/
void loadRun(int runNum);
void runFFT(int padNum);

/***** "Private" variables *****/
TFile *oFile = nullptr;
TCanvas *cFFT = nullptr;
Int_t fFreqMin = 512;
Int_t fFreqMax = 512;
Int_t fAbsMin = 5;
Int_t fInitTransition = 2;
Int_t fFreqTransition = 5;
Double_t *fact = new Double_t[512];
AtFFTFilter *fftFilter = new AtFFTFilter();

void fftDev(int runNum = 195)
{
   loadRun(runNum);
   fftFilter->Init();
}
void loadRun(int runNum)
{
   TString filePath = "/mnt/analysis/e12014/TPC/unpacked/run_%04d.root";
   loadRun(TString::Format(filePath, runNum), "AtRawEvent", "AtEventH");

   if(oFile != nullptr)
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
   cFFT->cd(1);
   hTrace->Draw();
   std::cout << hTrace->GetNbinsX() << std::endl;

   auto pad = rawEventPtr->GetPad(padNum);
   fftFilter->InitEvent(nullptr);
   cFFT->cd(3);
   fftFilter->Filter(pad);
   return;

   auto n = 512;
   auto fft = TVirtualFFT::FFT(1, &n, "R2CFORWARD M K");
   Double_t *input = new Double_t[n];
   for (int i = 0; i < n; ++i)
      input[i] = hTrace->GetBinContent(i + 1);
   fft->SetPoints(input);
   fft->Transform();

   Double_t *re_full = new Double_t[n];
   Double_t *im_full = new Double_t[n];
   Double_t *mag_full = new Double_t[n];
   for (int i = 0; i < n; ++i) {
      re_full[i] = 0;
      im_full[i] = 0;
      mag_full[i] = 0;
   }
   TH1D *hTrans = new TH1D("trans", "Magnitude", n, 0, n);
   for (int i = 0; i < n; ++i) {
      fft->GetPointComplex(i, re_full[i], im_full[i]);
      re_full[i] /= (double)n;
      im_full[i] /= (double)n;
      mag_full[i] = std::sqrt(std::pow(re_full[i], 2) + std::pow(im_full[i], 2));
      hTrans->SetBinContent(i + 1, mag_full[i]);
   }
   cFFT->cd(3);
   hTrans->Draw();

   // Use the following method to get the full output:
   Double_t *re_full_cut = new Double_t[n];
   Double_t *im_full_cut = new Double_t[n];

   // Add freq cuts

   TH1D *hTransCut = new TH1D("transCut", "Magnitude with cuttoff", n, 0, n);
   for (int i = 0; i < n; ++i) {

      cout << i << " " << fact[i] << endl;
      re_full_cut[i] = fact[i] * re_full[i];
      im_full_cut[i] = fact[i] * im_full[i];
      auto mag = std::sqrt(re_full_cut[i] * re_full_cut[i] + im_full_cut[i] * im_full_cut[i]);

      hTransCut->SetBinContent(i + 1, mag);
   }
   cFFT->cd(4);
   hTransCut->Draw();

   // Now let's make a backward transform:
   TVirtualFFT *fft_back = TVirtualFFT::FFT(1, &n, "C2R M K");
   fft_back->SetPointsComplex(re_full_cut, im_full_cut);
   fft_back->Transform();
   TH1 *hb = nullptr;
   hb = TH1::TransformHisto(fft_back, hb, "Re");
   hb->SetTitle("The backward transform result");
   cFFT->cd(2);
   hb->Draw();
}
