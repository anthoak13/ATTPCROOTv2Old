#include "AtAuxFilterTask.h"

#include <FairLogger.h>
#include <FairRootManager.h>
#include <FairTask.h>

#include <TString.h>

#include <algorithm>
#include <ostream>
#include <utility>

// Root Classes
#include "AtFilter.h"
#include "AtRawEvent.h"
// AtTPCRoot Classes
#include "AtAuxPad.h"
#include "AtPad.h"

#include <TClonesArray.h>
#include <TObject.h> // for TObject

AtAuxFilterTask::AtAuxFilterTask(AtFilter *filter) : fFilter(filter), fInputEventBranchName("AtRawEvent") {}

AtAuxFilterTask::~AtAuxFilterTask() = default;

void AtAuxFilterTask::AddAuxPad(std::string pad)
{
   auxPads.emplace_back(pad);
}

void AtAuxFilterTask::SetInputBranchName(TString branchName)
{
   fInputEventBranchName = branchName;
}

InitStatus AtAuxFilterTask::Init()
{
   FairRootManager *ioManager = FairRootManager::Instance();

   if (ioManager == nullptr) {
      LOG(ERROR) << "Cannot find RootManager!" << std::endl;
      return kERROR;
   }

   // Get the old data from the io manager
   fInputEventArray = dynamic_cast<TClonesArray *>(ioManager->GetObject(fInputEventBranchName));
   if (fInputEventArray == nullptr) {
      LOG(ERROR) << "AtAuxFilterTask: Cannot find AtRawEvent array " << fInputEventBranchName;
      return kERROR;
   }

   fFilter->Init();

   return kSUCCESS;
}

void AtAuxFilterTask::Exec(Option_t *opt)
{
   if (fInputEventArray->GetEntriesFast() == 0)
      return;

   auto *rawEvent = dynamic_cast<AtRawEvent *>(fInputEventArray->At(0));
   if (!rawEvent->IsGood())
      return;

   fFilter->InitEvent(rawEvent);

   for (auto &auxName : auxPads) {
      auto unfilteredPad = rawEvent->GetAuxPad(auxName);
      if (unfilteredPad == nullptr) {
         LOG(error) << "AtRawEvent does not contain aux pad: " << auxName;
         continue;
      }

      auto itPair = rawEvent->AddAuxPad(auxName + "Filtered");

      auto pad = itPair.first;
      pad->SetRawADC(unfilteredPad->GetRawADC());
      if (unfilteredPad->IsPedestalSubtracted()) {
         pad->SetADC(unfilteredPad->GetADC());
         pad->SetPedestalSubtracted();
      }
      fFilter->Filter(pad);
   }

   rawEvent->SetIsGood(fFilter->IsGoodEvent());
}
