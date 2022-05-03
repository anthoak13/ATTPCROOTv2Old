/**
 * A base class that describes the pattern of a track. For example a 2D circle or a line
 * These patterns are parameterized by some number (t)
 *
 */
#ifndef ATTRACKMODEL_H
#define ATTRACKMODEL_H

#include "AtHit.h"

#include <TObject.h>

class AtTrackModel : public TObject {
protected:
   std::vector<Double_t> fModelPar; //< Description of model
   Double_t fChi2{NAN};             //< How good the model is at describing the data
   Int_t fNFree{0};                 //< Degrees of freedom in the fit to the model
   const Int_t fNumPoints;          //< Number of 3D points that define the model (i.e. size of fIndices)

public:
   AtTrackModel(Int_t numPoints = 0);
   Double_t FitModel(const std::vector<AtHit> &pointsToFit, Double_t qThreshold = -1);
   Double_t FitModel(const std::vector<XYZPoint> &pointsToFit);

   /**
    * @brief Construct a model.
    *
    * Will set-up the model using the vector of XYZPoints. Assumes the size of points
    * is equal to fNumPoints
    *
    * @param[in] points 3D points to use when constructing the model
    */
   virtual void ConstructModel(const std::vector<XYZPoint> &points) = 0;

   /**
    * @brief Closest distance to model.
    *
    * @param[in] point Point to get the distance from.
    */
   virtual Double_t DistanceToModel(const XYZPoint &point) = 0;
   /**
    * @brief Closest point on model.
    *
    * @param[in] point Point to get the closest point on the model.
    * @return Closest point on model
    */
   virtual XYZPoint ClosestPointOnModel(const XYZPoint &point) = 0;
   /**
    * @brief Point on model at t
    *
    * Get the point on the model at parameter t. What t physically represents is pattern dependent.
    */
   virtual XYZPoint GetPointAt(double t) = 0;

   /**
    * @brief Number of points to define the model.
    *
    * The minimum number of unique points needed to define the model.
    */
   Int_t GetNumPoints() const { return fNumPoints; }
   Double_t GetChi2() const { return fChi2; }
   Int_t GetNFree() const { return fNFree; }
   std::vector<double> GetModelPar() const { return fModelPar; }
   void SetChi2(double chi2) { fChi2 = chi2; }

protected:
   /**
    * Called by other versions of FitModel. If pointCharge is not empty does charge weighted fit.
    * Sets fModelPar, fChi2, and fNFree
    */
   virtual void FitModel(const std::vector<XYZPoint> &pointsToFit, const std::vector<double> &pointCharge) = 0;

   ClassDef(AtTrackModel, 1)
};

#endif //#ifndef ATTRACKMODEL_H
