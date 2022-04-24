#ifndef AtEVENT_H
#define AtEVENT_H

#include "AtAuxPad.h"
#include "AtHit.h"

#include <FairLogger.h>

#include <Rtypes.h>
#include <TNamed.h>

#include <algorithm>
#include <array>
#include <map>
#include <memory>
#include <utility>
#include <vector>

class TBuffer;
class TClass;
class TMemberInspector;

class AtRawEvent;

using auxPadVector = std::vector<AtAuxPad>;
using traceArray = std::array<Float_t, 512>;
using hitPtr = std::unique_ptr<AtHit>;
using hitVector = std::vector<AtHit>;

class AtEvent : public TNamed {
private:
   Int_t fEventID;
   Bool_t fIsGood;
   Bool_t fIsInGate;
   ULong_t fTimestamp;
   Double_t fEventCharge = -100;
   Double_t fRhoVariance = 0;

   hitVector fHitArray;
   std::vector<AtAuxPad> fAuxPadArray;
   std::map<Int_t, Int_t> fMultiplicityMap;

   traceArray fMeshSig{};

public:
   AtEvent();
   AtEvent(Int_t eventID, Bool_t isGood, Bool_t isInGate = false, ULong_t timestamp = 0);
   AtEvent(const AtEvent &copy) = default;
   AtEvent(const AtRawEvent &copy);
   ~AtEvent() = default;

   void Clear(Option_t *opt = nullptr) override;

   // Copies everything except the hit array from the passed AtEvent
   void CopyFrom(const AtEvent &event);

   // Adds a new hit to the hit array, and returns a referece to the new hit to be
   // filled. This is done to avoid the create and subsequent copy of a hit and also
   // avoid dealing with the memory managment
   // Takes arguments to any constructor of AtHit, leaving out the first (the hit ID).
   // AtEvent handles the assignment of hit IDs to ensure they are unique within an event.
   template <typename... Ts>
   AtHit &AddHit(Ts &&...params)
   {
      LOG(debug) << "Adding hit with ID " << fHitArray.size() << " to event " << fEventID;
      fHitArray.emplace_back(fHitArray.size(), std::forward<Ts>(params)...);
      return fHitArray.back();
   }

   // Clones the hit and adds it to the event, setting the hitID to the next valid for
   // this event. It is not coppied from the passed hit.
   AtHit &AddHit(const AtHit &hit)
   {
      LOG(debug) << "Adding hit with ID " << fHitArray.size() << " to event " << fEventID;
      fHitArray.emplace_back(hit);
      fHitArray.back().SetHitID(fHitArray.size() - 1);
      return fHitArray.back();
   }

   AtHit &AddHit(AtHit &hit) { return AddHit(const_cast<const AtHit &>(hit)); }

   // Copies passed aux pad into the event's auxiliary pad array
   void AddAuxPad(AtAuxPad auxPad) { fAuxPadArray.push_back(std::move(auxPad)); }

   void SetEventID(Int_t evtid) { fEventID = evtid; }
   void SetTimestamp(ULong_t timestamp) { fTimestamp = timestamp; }
   void SetEventCharge(Double_t Qevent) { fEventCharge = Qevent; }
   void SetRhoVariance(Double_t RhoVariance) { fRhoVariance = RhoVariance; }
   void SetIsGood(Bool_t value) { fIsGood = value; }
   void SetIsInGate(Bool_t value) { fIsInGate = value; }

   void SetMultiplicityMap(std::map<Int_t, Int_t> MultiMap) { fMultiplicityMap = std::move(MultiMap); }
   void SetMeshSignal(const traceArray &mesharray);
   void SetMeshSignal(Int_t idx, Float_t val);

   const AtHit &GetHit(Int_t hitNo) const { return fHitArray.at(hitNo); }
   const hitVector &GetHitArray() const { return fHitArray; }
   const auxPadVector &GetAuxPadArray() const { return fAuxPadArray; }
   Int_t GetEventID() const { return fEventID; }
   Long_t GetTimestamp() const { return fTimestamp; }
   Int_t GetNumHits() const { return fHitArray.size(); }
   Double_t GetEventCharge() const { return fEventCharge; }
   Double_t GetRhoVariance() const { return fRhoVariance; }
   const traceArray &GetMesh() const { return fMeshSig; }
   Int_t GetHitPadMult(Int_t PadNum); // Returns the multiplicity of the pad where this hit belongs to
   const std::map<Int_t, Int_t> &GetMultiMap() { return fMultiplicityMap; }

   Bool_t IsGood() const { return fIsGood; }
   Bool_t IsInGate() const { return fIsInGate; }

   void SortHitArray();
   void SortHitArrayTime();

   ClassDefOverride(AtEvent, 4);
};

#endif
