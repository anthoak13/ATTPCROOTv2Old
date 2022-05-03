#ifndef ATVERTEXFINDER_H
#define ATVERTEXFINDER_H

#include "AtTrack.h"

#include <Math/Point3D.h>

#include <vector>

class TVector3;

class AtVertexFinder {
protected:
   using XYZPoint = ROOT::Math::XYZPoint;
   double fLineDistThreshold;

public:
   XYZPoint FindVertex(const std::vector<AtTrack> &tracks);
   std::vector<AtTrack> FindVertexOnePerTrack(const std::vector<AtTrack> &tracks);

private:
   TVector3 ClosestPoint2Lines(TVector3 d1, TVector3 pt1, TVector3 d2, TVector3 pt2);
};

#endif
