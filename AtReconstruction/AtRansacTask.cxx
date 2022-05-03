#include "AtRansacTask.h"

#include "AtEvent.h"     // for AtEvent
#include "AtLmedsMod.h"  // for AtLmedsMod
#include "AtMlesacMod.h" // for AtMlesacMod
#include "AtPatternEvent.h"
#include "AtRansac.h"    // for AtRansac
#include "AtRansacMod.h" // for AtRansacMod

#include <FairLogger.h>      // for LOG, Logger
#include <FairRootManager.h> // for FairRootManager

#include <TClonesArray.h> // for TClonesArray
#include <TObject.h>      // for TObject

#include <memory> // for allocator

ClassImp(AtRansacTask);

AtRansacTask::AtRansacTask()
   : fInputBranchName("AtEventH"), fOutputBranchName("AtRansac"), kIsPersistence(kFALSE), kIsFullMode(kFALSE),
     kIsReprocess(kFALSE), fPatternEventArray("AtPatternEvent", 1)
{
}

AtRansacTask::~AtRansacTask() = default;

void AtRansacTask::SetPersistence(Bool_t value)
{
   kIsPersistence = value;
}
void AtRansacTask::SetInputBranch(TString branchName)
{
   fInputBranchName = branchName;
}

void AtRansacTask::SetOutputBranch(TString branchName)
{
   fOutputBranchName = branchName;
}

void AtRansacTask::SetModelType(int model)
{
   fRANSACModel = model;
}
void AtRansacTask::SetDistanceThreshold(Float_t threshold)
{
   fRANSACThreshold = threshold;
}
void AtRansacTask::SetFullMode()
{
   kIsFullMode = kTRUE;
}
void AtRansacTask::SetMinHitsLine(Int_t nhits)
{
   fMinHitsLine = nhits;
}
void AtRansacTask::SetTiltAngle(Double_t val)
{
   fTiltAngle = val;
}
void AtRansacTask::SetNumItera(Int_t niterations)
{
   fNumItera = niterations;
}
void AtRansacTask::SetAlgorithm(Int_t val)
{
   fRANSACAlg = val;
}
void AtRansacTask::SetRanSamMode(Int_t mode)
{
   fRandSamplMode = mode;
};
void AtRansacTask::SetIsReprocess(Bool_t value)
{
   kIsReprocess = value;
}
void AtRansacTask::SetChargeThreshold(Double_t value)
{
   fCharThres = value;
}
void AtRansacTask::SetVertexMode(Int_t value)
{
   fVertexMode = value;
}
void AtRansacTask::SetInputBranchName(TString inputName)
{
   fInputBranchName = inputName;
}
void AtRansacTask::SetOutputBranchName(TString outputName)
{
   fOutputBranchName = outputName;
}

InitStatus AtRansacTask::Init()
{

   if (fRANSACAlg == 0)
      fRansacArray = new TClonesArray("AtRANSACN::AtRansac");
   else if (fRANSACAlg == 1)
      fRansacArray = new TClonesArray("AtRansacMod");
   else if (fRANSACAlg == 2)
      fRansacArray = new TClonesArray("AtMlesacMod");
   else if (fRANSACAlg == 3)
      fRansacArray = new TClonesArray("AtLmedsMod");
   else {
      LOG(error) << "Cannot find Ransac algorithm!";
      return kERROR;
   }

   FairRootManager *ioMan = FairRootManager::Instance();
   if (ioMan == nullptr) {
      LOG(error) << "Cannot find RootManager!";
      return kERROR;
   }

   fEventArray = dynamic_cast<TClonesArray *>(ioMan->GetObject(fInputBranchName));
   if (fEventArray == nullptr) {

      LOG(error) << "Cannot find AtEvent array!";
      return kERROR;
   }

   if (fRANSACModel == 0)
      ioMan->Register(fOutputBranchName, "AtTPC", fRansacArray, kIsPersistence);
   else
      ioMan->Register("AtPatternEvent", "AtTPC", &fPatternEventArray, kIsPersistence);

   if (kIsReprocess) {

      ioMan->Register("AtEventH", "AtTPC", fEventArray, kIsPersistence);
   }

   return kSUCCESS;
}

void AtRansacTask::Exec(Option_t *opt)
{

   fRansacArray->Delete();

   if (fEventArray->GetEntriesFast() == 0)
      return;

   fEvent = dynamic_cast<AtEvent *>(fEventArray->At(0));

   fPatternEventArray.Delete();
   auto *patternEvent = (AtPatternEvent *)new (fPatternEventArray[0]) AtPatternEvent();

   LOG(debug) << "Running RANSAC with " << fEvent->GetNumHits() << " hits.";

   if (fRANSACAlg == 0) {
      LOG(debug) << "Running RANSAC algorithm AtRANSACN::AtRansac";
      auto *Ransac = (AtRANSACN::AtRansac *)new ((*fRansacArray)[0]) AtRANSACN::AtRansac();
      Ransac->SetTiltAngle(fTiltAngle);
      if (fRANSACModel != -1)
         Ransac->SetModelType(fRANSACModel);
      Ransac->SetDistanceThreshold(fRANSACThreshold);
      Ransac->SetMinHitsLine(fMinHitsLine);
      if (kIsFullMode)
         Ransac->CalcRANSACFull(fEvent);
      else
         Ransac->CalcRANSAC(fEvent);
   } else {
      LOG(debug) << "Running Unified RANSAC";
      AtRansacMod ransac;
      ransac.SetDistanceThreshold(fRANSACThreshold);
      ransac.SetMinHitsLine(fMinHitsLine);
      ransac.SetNumItera(fNumItera);
      ransac.SetRanSamMode(static_cast<AtRandomSample::SampleMethod>(fRandSamplMode));
      ransac.Solve(fEvent, patternEvent);
      ransac.SetChargeThres(fCharThres);
   }
}
