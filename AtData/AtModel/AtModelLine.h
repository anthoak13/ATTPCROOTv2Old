#ifndef ATMODELLINE_H
#define ATMODELLINE_H

#include "AtTrackModel.h"

#include <Math/Point3D.h>
using XYZPoint = ROOT::Math::XYZPoint;

class AtModelLine : public AtTrackModel {
public:
   AtModelLine();

   XYZPoint GetPoint() { return {fModelPar[0], fModelPar[1], fModelPar[2]}; }
   XYZVector GetDirection() { return {fModelPar[3], fModelPar[4], fModelPar[5]}; }

   virtual void ConstructModel(const std::vector<XYZPoint> &points) override;
   virtual Double_t DistanceToModel(const XYZPoint &point) override;
   virtual XYZPoint ClosestPointOnModel(const XYZPoint &point) override;

protected:
   virtual void FitModel(const std::vector<XYZPoint> &points, const std::vector<double> &charge) override;

   ClassDefOverride(AtModelLine, 1)
};

#endif //#ifndef ATMODELLINE_H
