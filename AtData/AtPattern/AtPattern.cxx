#include "AtPattern.h"

#include "AtPatternCircle2D.h"
#include "AtPatternLine.h"
#include "AtPatternTypes.h"

#include <FairLogger.h>

#include <TMath.h>
#include <TRandom3.h>

ClassImp(AtPattern);
// Defined here so we don't force a recompule

AtPattern::AtPattern(Int_t numPoints) : fNumPoints(numPoints) {}

/**
 * @brief Fit the pattern.
 *
 * Fit the pattern shape using all hits in pointsToFit. Does a charge weighted fit if qThreshold is not equal to -1.
 * Sets fPatternPar parameters, and fChi2.
 *
 * @param[in] points Points in 3D space to fit with charge information
 * @param[in] qThreshold Only fit points that are above this charge threshold.
 * @return Chi-squared of the fit
 */
Double_t AtPattern::FitPattern(const std::vector<AtHit> &pointsToFit, Double_t qThreshold)
{
   std::vector<XYZPoint> points;
   std::vector<double> charge;
   for (const auto &hit : pointsToFit) {
      if (hit.GetCharge() > qThreshold) {
         points.push_back(hit.GetPosition());
         charge.push_back(hit.GetCharge());
      }
   }

   if (qThreshold == -1)
      FitPattern(points);
   else
      FitPattern(points, charge);

   return fChi2;
}

/**
 * @brief Fit the pattern shape.
 *
 * Fit the pattern shape using all points in pointsToFit. Does not weight for charge.
 * Sets fPatternPar parameters, and fChi2.
 *
 * @param[in] pointsToFit Points in 3D space to fit
 * @return Chi-squared of the fit
 */
Double_t AtPattern::FitPattern(const std::vector<XYZPoint> &pointsToFit)
{
   std::vector<double> charge;
   FitPattern(pointsToFit, charge);
   return fChi2;
}

std::unique_ptr<AtPattern> AtPattern::CreatePattern(Type type)
{
   switch (type) {
   case (Type::kLine): return std::make_unique<AtPatternLine>();
   case (Type::kCircle2D): return std::make_unique<AtPatternCircle2D>();
   default: return nullptr;
   }
}
