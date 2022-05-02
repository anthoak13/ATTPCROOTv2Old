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

   virtual void ConstructModel(const std::vector<XYZPoint> &points) override;

   virtual Double_t DistanceToModel(const XYZPoint &point) override;

protected:
   virtual void FitModel(const std::vector<XYZPoint> &points, const std::vector<double> &charge) override;
};

#endif //#ifndef ATMODELLINE_H
