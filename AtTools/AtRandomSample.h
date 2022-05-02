#ifndef ATRANDOMSAMPLE_H
#define ATRANDOMSAMPLE_H

#include <Math/Point3D.h>

#include <vector>
class AtHit;

/**
 * Class to randomly sample AtHits, or vectors of arbitrary types
 */
class AtRandomSample {
public:
   enum class SampleMethod;
   using XYZPoint = ROOT::Math::XYZPoint;

public:
   static std::vector<XYZPoint> SamplePoints(int N, const std::vector<AtHit> &hits, SampleMethod mode);

protected:
   static std::vector<XYZPoint> sampleUniform(int N, const std::vector<AtHit> &hits);
   static std::vector<XYZPoint> sampleGaussian(int N, const std::vector<AtHit> &hits);
   static std::vector<XYZPoint> sampleWeighted(int N, const std::vector<AtHit> &hits);
   static std::vector<XYZPoint> sampleWeightedGaussian(int N, const std::vector<AtHit> &hits);
   static std::vector<double> getPDF(double &avgCharge, const std::vector<AtHit> &hits);
};

std::ostream &operator<<(std::ostream &os, const AtRandomSample::SampleMethod &t);

#endif //#ifndef ATRANDOMSAMPLE_H
