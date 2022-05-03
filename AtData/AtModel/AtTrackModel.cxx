#include "AtTrackModel.h"

#include <FairLogger.h>

#include <TMath.h>
#include <TRandom3.h>

ClassImp(AtTrackModel);

AtTrackModel::AtTrackModel(Int_t numPoints) : fNumPoints(numPoints) {}

/**
 * @brief Fit the model.
 *
 * Fit the model shape using all hits in pointsToFit. Does a charge weighted fit.
 * Sets fModelPar parameters, and fChi2.
 *
 * @param[in] points Points in 3D space to fit with charge information
 * @return Chi-squared of the fit
 */
Double_t AtTrackModel::FitModel(const std::vector<AtHit> &pointsToFit)
{
   std::vector<XYZPoint> points;
   std::vector<double> charge;
   for (const auto &hit : pointsToFit) {
      points.push_back(hit.GetPosition());
      charge.push_back(hit.GetCharge());
   }
   FitModel(points, charge);
   return fChi2;
}

/**
 * @brief Fit the model shape.
 *
 * Fit the model shape using all points in pointsToFit. Does not weight for charge.
 * Sets fModelPar parameters, and fChi2.
 *
 * @param[in] pointsToFit Points in 3D space to fit
 * @return Chi-squared of the fit
 */
Double_t AtTrackModel::FitModel(const std::vector<XYZPoint> &pointsToFit)
{
   std::vector<double> charge;
   FitModel(pointsToFit, charge);
   return fChi2;
}
