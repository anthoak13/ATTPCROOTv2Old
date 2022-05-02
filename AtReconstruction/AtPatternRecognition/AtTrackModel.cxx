#include "AtTrackModel.h"

#include <FairLogger.h>

#include <TMath.h>
#include <TRandom3.h>
enum class AtTrackModel::SampleMethod { kUniform = 0, kGaussian = 1, kWeighted = 2, kWeightedGaussian = 3 };

AtTrackModel::AtTrackModel(Int_t numPoints) : fNumPoints(numPoints) {}

void AtTrackModel::SetHitArray(const std::vector<AtHit> *hitArray)
{
   Reset();
   LOG(debug) << "Setting hit array with " << hitArray->size();
   fHitArray = hitArray;
}

void AtTrackModel::ConstructRandomModel(SampleMethod mode)
{
   ConstructModel(sampleModelPoints(mode));
}
/**
 * @brief Sample points required to define a model.
 *
 * Sample N points, where N is the minimum number of points required to define
 * the model being used (right now, that's a linear model so N = 2)
 *
 * @todo Generalize the functions being called so they return and arbitrary number
 * of points based on the requirements of the model
 *
 * @param[in] mode flag instructing what sampling algorithm to use
 * @return vector of indices of the points describing the model
 */
std::vector<int> AtTrackModel::sampleModelPoints(SampleMethod mode)
{
   switch (mode) {
   case (SampleMethod::kUniform): return sampleUniform();
   case (SampleMethod::kGaussian): return sampleGaussian();
   case (SampleMethod::kWeighted): return sampleWeighted();
   case (SampleMethod::kWeightedGaussian): return sampleWeightedGaussian();
   default: return {};
   }
}

/**
 * @brief Sample two random points from indX
 */
std::vector<int> AtTrackModel::sampleUniform()
{
   //-------Uniform sampling
   int ind1 = gRandom->Uniform(0, fHitArray->size());
   int ind2;
   do {
      ind2 = gRandom->Uniform(0, fHitArray->size());
   } while (ind1 == ind2);

   return {ind1, ind2};
}

/**
 * @brief Sample two points from indX
 *
 * The first point is sampled randomly. The second points is sampled according to
 * a gaussian distribition around the first point with a sigma of 30. If in 20 samples it
 * doesn't find a point close enough, it defaults to a uniform sample.
 *
 * @TODO We should probably set sigma so it scales with the size of the pad plane or make it
 * a tunable parameter.
 */
std::vector<int> AtTrackModel::sampleGaussian()
{
   //--------Gaussian sampling
   double sigma = 30.0;
   double y = 0;
   double gauss = 0;
   int counter = 0;
   int p1 = gRandom->Uniform(0, fHitArray->size());
   int p2;
   auto &P1 = fHitArray->at(p1).GetPosition();

   do {
      p2 = gRandom->Uniform(0, fHitArray->size());
      auto &P2 = fHitArray->at(p2).GetPosition();
      auto dist = std::sqrt((P2 - P1).Mag2());

      gauss = 1.0 * exp(-1.0 * pow(dist / sigma, 2.0));
      y = (gRandom->Uniform(0, 1));
      counter++;
      if (counter > 20 && p2 != p1)
         break;
   } while (p2 == p1 || y > gauss);

   return {p1, p2};
}
/**
 * @ brief Sample two points based on charge
 */
std::vector<int> AtTrackModel::sampleWeighted()
{
   //-------Weighted sampling
   auto Proba = getPDF();
   bool validSecondPoint = false;
   int counter = 0;
   int p1 = gRandom->Uniform(0, fHitArray->size());
   int p2 = p1;

   do {
      validSecondPoint = false;
      counter++;
      if (counter > 30 && p2 != p1)
         break;

      p2 = gRandom->Uniform(0, fHitArray->size());
      double TwiceAvCharge = 2 * fAvgCharge;

      if (Proba.size() == fHitArray->size())
         validSecondPoint = (Proba[p2] >= gRandom->Uniform(0, TwiceAvCharge));
      else
         validSecondPoint = true;

   } while (p2 == p1 || validSecondPoint == false);

   return {p1, p2};
}

/**
 * @brief Sample two points from indX
 *
 * The first point is sampled randomly by charge. The second points is sampled according to
 * a gaussian distribition around the first point with a sigma of 30. If in 20 samples it
 * doesn't find a point close enough, it defaults to a charge weighted sample.
 *
 * @TODO We should probably set sigma so it scales with the size of the pad plane or make it
 * a tunable parameter.
 */
std::vector<int> AtTrackModel::sampleWeightedGaussian()
{
   //-------Weighted sampling + Gauss dist.
   auto Proba = getPDF();
   bool validSecondPoint = false;
   double sigma = 30.0;
   double y = 0;
   double gauss = 0;
   int counter = 0;

   int p1 = gRandom->Uniform(0, fHitArray->size());
   int p2 = p1;
   auto &P1 = fHitArray->at(p1).GetPosition();

   do {
      p2 = gRandom->Uniform(0, fHitArray->size());
      auto &P2 = fHitArray->at(p2).GetPosition();
      auto dist = std::sqrt((P2 - P1).Mag2());

      gauss = 1.0 * exp(-1.0 * pow(dist / sigma, 2));
      y = (gRandom->Uniform(0, 1));

      counter++;
      if (counter > 30 && p2 != p1)
         break;

      validSecondPoint = false;
      double TwiceAvCharge = 2 * fAvgCharge;

      if (Proba.size() == fHitArray->size())
         validSecondPoint = (Proba[p2] >= gRandom->Uniform(0, TwiceAvCharge));
      else
         validSecondPoint = true;

   } while (p2 == p1 || validSecondPoint == false || y > gauss);

   return {p1, p2};
}

std::vector<double> AtTrackModel::getPDF()
{
   double Tcharge = 0;
   for (const auto &hit : *fHitArray)
      Tcharge += hit.GetCharge();

   fAvgCharge = Tcharge / fHitArray->size();
   std::vector<double> w;
   if (Tcharge > 0)
      for (const auto &hit : *fHitArray)
         w.push_back(hit.GetCharge() / Tcharge);

   return w;
}

void AtTrackModel::Reset()
{
   fHitArray = nullptr;
   fIndices.clear();
   fAvgCharge = 0;
}
