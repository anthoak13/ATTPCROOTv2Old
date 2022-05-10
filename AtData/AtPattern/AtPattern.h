#ifndef ATPATTERN_H
#define ATPATTERN_H

#include "AtHit.h"

#include <Rtypes.h> // for Double_t, Int_t, THashConsistencyHolder, ClassDef
#include <TObject.h>

#include <math.h> // for NAN

#include <algorithm> // for max
#include <vector>    // for vector

class TBuffer;
class TClass;
class TMemberInspector;

/**
 * @defgroup AtPattern Track Patterns
 * @brief Track shapes in 3D space.
 *
 * Collection of classes for describing and drawing track patterns (essentially the shape of a track).
 * These objects are stored as the fit for AtTrack and are used to run sample consensus models (\ref SampleConsensus).
 *
 */
namespace AtPatterns {

/**
 * @brief Describes a shape in 3D space.
 *
 * A base class that describes the pattern of a track. For example a 2D circle or a line
 * These patterns are parameterized by some number (t)
 *
 * @ingroup AtPattern
 */
class AtPattern : public TObject {
protected:
   std::vector<Double_t> fPatternPar; //< Description of pattern
   Double_t fChi2{NAN};               //< How good the pattern is at describing the data
   Int_t fNFree{0};                   //< Degrees of freedom in the fit to the pattern
   const Int_t fNumPoints;            //< Number of 3D points that define the pattern (i.e. size of fIndices)

public:
   AtPattern(Int_t numPoints = 0);
   Double_t FitPattern(const std::vector<AtHit> &pointsToFit, Double_t qThreshold = -1);
   Double_t FitPattern(const std::vector<XYZPoint> &pointsToFit);

   /**
    * @brief Define based on points.
    *
    * Will set-up the pattern using the vector of XYZPoints. Assumes the size of points
    * is equal to fNumPoints
    *
    * @param[in] points 3D points to use when defining the pattern
    */
   virtual void DefinePattern(const std::vector<XYZPoint> &points) = 0;

   /**
    * @brief Closest distance to pattern.
    *
    * @param[in] point Point to get the distance from.
    */
   virtual Double_t DistanceToPattern(const XYZPoint &point) = 0;
   /**
    * @brief Closest point on pattern.
    *
    * @param[in] point Point to get the closest point on the pattern.
    * @return Closest point on pattern
    */
   virtual XYZPoint ClosestPointOnPattern(const XYZPoint &point) = 0;
   /**
    * @brief Point on pattern at t
    *
    * Get the point on the pattern at parameter t. What t physically represents is pattern dependent.
    */
   virtual XYZPoint GetPointAt(double t) = 0;

   /**
    * @brief Number of points to define the pattern.
    *
    * The minimum number of unique points needed to define the pattern.
    */
   Int_t GetNumPoints() const { return fNumPoints; }
   Double_t GetChi2() const { return fChi2; }
   Int_t GetNFree() const { return fNFree; }
   std::vector<double> GetPatternPar() const { return fPatternPar; }
   void SetChi2(double chi2) { fChi2 = chi2; }

protected:
   /**
    * Called by other versions of FitPattern. If pointCharge is not empty does charge weighted fit.
    * Sets fPatternPar, fChi2, and fNFree
    */
   virtual void FitPattern(const std::vector<XYZPoint> &pointsToFit, const std::vector<double> &pointCharge) = 0;

   ClassDef(AtPattern, 1)
};
} // namespace AtPatterns
#endif //#ifndef ATTRACKPATTERN_H
