/**
 * A base class that describes a model of a track. For example a 2D circle or a line
 *
 */
#ifndef ATTRACKMODEL_H
#define ATTRACKMODEL_H

#include "AtHit.h"

#include <TObject.h>

class AtTrackModel : public TObject {
public:
   enum class SampleMethod;

protected:
   Double_t fAvgCharge;

   const std::vector<AtHit> *fHitArray; //< Pointer to the hit cloud we are working with
   std::vector<int> fIndices;           //< Vector containing the indcices of fHitArray defining the model
   const Int_t fNumPoints;              //< Number of points that define the model (i.e. size of fIndices)

public:
   AtTrackModel(Int_t numPoints = 0);
   /**
    * @brief Set the hit array (point cloud) for this model.
    */
   void SetHitArray(const std::vector<AtHit> *hitArray);
   /**
    * @brief Indices of fHitArray that define the model
    */
   std::vector<int> GetIndices() const { return fIndices; };

   Int_t GetNumHits() { return fHitArray->size(); }
   /**
    * @brief Construct a model by randomly sampling fHitArray.
    *
    * Will randomly sample fHitArray and set-up the model using the sampled points.
    *
    * @param[in] mode The method to use while sampling
    */
   virtual void ConstructRandomModel(SampleMethod mode);
   /**
    * @brief Construct a model.
    *
    * Will set-up the model using using the passed indices of fHitArray.
    *
    * @param[in] mode The method to use while sampling
    */
   virtual void ConstructModel(const std::vector<int> &idx) = 0;

   /**
    * @brief Closest distance from hit to model.
    *
    * @param[in] hit AtHit to get the distance from.
    */
   virtual Double_t DistanceToModel(const AtHit &hit) = 0;
   /**
    * @brief Closest distance from hit to model.
    *
    * @param[in] hitIndex Index of AtHit in fHitArray to get the distance from.
    */
   virtual Double_t DistanceToModel(Int_t hitIndex) { return DistanceToModel(fHitArray->at(hitIndex)); }
   /**
    * @brief Fit the model shape.
    *
    * Fit the model shape using all points in idx.
    *
    * @param[in] idx Indices of fHitArray to fit3d
    * @param[out] fitPar Fit parameters for this model
    * @return Chi-squared of the fit
    */
   virtual Double_t Fit3D(const std::vector<int> idx, std::vector<double> &fitPar) = 0;

protected:
   std::vector<int> sampleModelPoints(SampleMethod mode);
   std::vector<int> sampleUniform();
   std::vector<int> sampleGaussian();
   std::vector<int> sampleWeighted();
   std::vector<int> sampleWeightedGaussian();
   std::vector<double> getPDF(); //< Returns the PDF weighted by charge for fHitArray.

   /**
    * Set the internal model parameters used by other functions based on fIndices
    */
   virtual void Reset();
};

#endif //#ifndef ATTRACKMODEL_H
