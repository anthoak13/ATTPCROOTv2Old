#ifndef ATPATTERNCIRCLE2D_H
#define ATPATTERNCIRCLE2D_H

#include "AtPattern.h"

#include <Math/Point3D.h>
using XYZPoint = ROOT::Math::XYZPoint;

class AtPatternCircle2D : public AtPattern {
public:
   AtPatternCircle2D();

   XYZPoint GetCenter() { return {fPatternPar[0], fPatternPar[1], 0}; }
   double GetRadius() { return fPatternPar[2]; }

   virtual void DefinePattern(const std::vector<XYZPoint> &points) override;
   virtual Double_t DistanceToPattern(const XYZPoint &point) override;
   virtual XYZPoint ClosestPointOnPattern(const XYZPoint &point) override;
   virtual XYZPoint GetPointAt(double theta) override;

protected:
   virtual void FitPattern(const std::vector<XYZPoint> &points, const std::vector<double> &charge) override;

   ClassDefOverride(AtPatternCircle2D, 1)
};

#endif //#ifndef ATPATTERNCIRCLE2D_H
