#ifndef ATPATTERNLINE_H
#define ATPATTERNLINE_H

#include "AtPattern.h"

#include <Math/Point3D.h>
using XYZPoint = ROOT::Math::XYZPoint;

class AtPatternLine : public AtPattern {
public:
   AtPatternLine();

   XYZPoint GetPoint() { return {fPatternPar[0], fPatternPar[1], fPatternPar[2]}; }
   XYZVector GetDirection() { return {fPatternPar[3], fPatternPar[4], fPatternPar[5]}; }

   virtual void DefinePattern(const std::vector<XYZPoint> &points) override;
   virtual Double_t DistanceToPattern(const XYZPoint &point) override;
   virtual XYZPoint ClosestPointOnPattern(const XYZPoint &point) override;
   virtual XYZPoint GetPointAt(double z) override;

protected:
   virtual void FitPattern(const std::vector<XYZPoint> &points, const std::vector<double> &charge) override;

   ClassDefOverride(AtPatternLine, 1)
};

#endif //#ifndef ATPATTERNLINE_H
