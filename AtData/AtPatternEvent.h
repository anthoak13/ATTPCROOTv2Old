#ifndef ATPATTERNEVENT_H
#define ATPATTERNEVENT_H

#include "AtTrack.h"

#include <Rtypes.h>
#include <TNamed.h>

#include <algorithm>
#include <vector>
class TBuffer;
class TClass;
class TMemberInspector;

class AtPatternEvent : public TNamed {

public:
   AtPatternEvent();
   ~AtPatternEvent();

   void SetTrackCand(std::vector<AtTrack> tracks);
   void AddTrack(const AtTrack &track) { fTrackCand.push_back(track); }
   void AddTrack(AtTrack &&track) { fTrackCand.push_back(track); }

   std::vector<AtTrack> &GetTrackCand();

private:
   std::vector<AtTrack> fTrackCand; // Candidate tracks

   ClassDef(AtPatternEvent, 1);
};

#endif
