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

class AtLmedsMod : public AtRansacMod {
protected:
   virtual int evaluateModel(AtTrackModel *model, const std::vector<AtHit> &hitArray) override;

   double GetMedian(std::vector<double> &errvec);
};

#endif
