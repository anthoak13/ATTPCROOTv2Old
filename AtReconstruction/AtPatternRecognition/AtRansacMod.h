/*******************************************************************
// Basic RANSAC Class                                              *
// Author: J.C. Zamora, jczamorac@gmail.com                        *
// University of Sao Paulo, 26-08-2020                             *
********************************************************************/

#ifndef AtRANSACMOD_H
#define AtRANSACMOD_H

#include "AtSampleConsensus.h"

class AtEvent;
class AtPatternEvent;

class AtRansacMod : public AtSampleConsensus {
protected:
   virtual int evaluateModel(AtTrackModel *model, const std::vector<AtHit> &hitsToCheck);
};

#endif
