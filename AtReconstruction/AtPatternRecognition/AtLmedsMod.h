/*******************************************************************
// Basic Lmeds Class                                              *
// Author: J.C. Zamora, jczamorac@gmail.com                        *
// University of Sao Paulo, 26-08-2020                             *
********************************************************************/

#ifndef AtLmedsMOD_H
#define AtLmedsMOD_H

#ifdef _OPENMP
#include <omp.h>
#endif

#include "AtRansacMod.h"
#include "AtTrack.h" // for AtTrack

#include <Rtypes.h>   // for Int_t, Double_t, THashConsistencyHolder, ClassDef
#include <TObject.h>  // for TObject
#include <TVector3.h> // for TVector3

#include <stdio.h> // for size_t

#include <algorithm> // for max
#include <utility>   // for pair
#include <vector>    // for vector
class AtEvent;
class TBuffer;
class TClass;
class TMemberInspector;

class AtLmedsMod : public AtRansacMod {
protected:
   std::vector<double> errorsVec;

public:
   AtLmedsMod() = default;
   ~AtLmedsMod() = default;

protected:
   virtual void Reset() override;
   virtual void Solve() override;
   double GetMedian(std::vector<double> errvec);
   ClassDefOverride(AtLmedsMod, 2);
};

#endif
