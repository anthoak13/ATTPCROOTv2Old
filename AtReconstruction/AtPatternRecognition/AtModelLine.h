#ifndef ATMODELLINE_H
#define ATMODELLINE_H

#include "AtTrackModel.h"

#include <Math/Point3D.h>
using XYZPoint = ROOT::Math::XYZPoint;

class AtModelLine : public AtTrackModel {
protected:
   XYZPoint fPoint;
   XYZVector fDirection;

public:
   AtModelLine();

   virtual Double_t Fit3D(const std::vector<int> idx, std::vector<double> &fitPar) override;
   virtual void ConstructModel(const std::vector<int> &idx) override;
   virtual Double_t DistanceToModel(const AtHit &hit) override;
};

#endif //#ifndef ATMODELLINE_H
