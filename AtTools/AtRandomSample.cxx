#include "AtRandomSample.h"

#include "AtHit.h"

#include <FairLogger.h>

#include <TRandom3.h>

#include <numeric>

enum class AtRandomSample::SampleMethod { kUniform = 0, kGaussian = 1, kWeighted = 2, kWeightedGaussian = 3 };
std::ostream &operator<<(std::ostream &os, const AtRandomSample::SampleMethod &t)
{
   using method = AtRandomSample::SampleMethod;
   switch (t) {
   case (method::kUniform): os << "SampleMethod::kUniform"; break;
   case (method::kGaussian): os << "SampleMethod::kGaussian"; break;
   case (method::kWeighted): os << "SampleMethod::kWeighted"; break;
   case (method::kWeightedGaussian): os << "SampleMethod::kWeightedGaussian"; break;
   default: os << "SampleMethod::Other";
   }
   return os;
}

/**
 * @brief Sample points from hits
 *
 * Sample N unique points from hits using the method mode
 *
 * @todo Generalize the functions being called so they return and arbitrary number
 * of points based on the requirements of the model
 *
 * @param[in] mode flag instructing what sampling algorithm to use
 * @return vector of sampled points
 */
std::vector<AtRandomSample::XYZPoint>
AtRandomSample::SamplePoints(int N, const std::vector<AtHit> &hits, SampleMethod mode = SampleMethod::kUniform)
{
   switch (mode) {
   case (SampleMethod::kUniform): return sampleUniform(N, hits);
   case (SampleMethod::kGaussian): return sampleGaussian(N, hits);
   case (SampleMethod::kWeighted): return sampleWeighted(N, hits);
   case (SampleMethod::kWeightedGaussian): return sampleWeightedGaussian(N, hits);
   default: LOG(error) << "Invalid sample method passed " << mode; return {};
   }
}
/**
 * @brief Sample random points
 */
std::vector<AtRandomSample::XYZPoint> AtRandomSample::sampleUniform(int N, const std::vector<AtHit> &hits)
{
   //-------Uniform sampling
   int ind1 = gRandom->Uniform(0, hits.size());
   int ind2;
   do {
      ind2 = gRandom->Uniform(0, hits.size());
   } while (ind1 == ind2);

   return {hits.at(ind1).GetPosition(), hits.at(ind2).GetPosition()};
}

/**
 * @brief Sample points from
 *
 * The first point is sampled randomly. The remaining points are sampled according to
 * a gaussian distribition around the first point with a sigma of 30. If in 20 samples it
 * doesn't find a point close enough, it defaults to a uniform sample.
 *
 * @TODO We should probably set sigma so it scales with the size of the pad plane or make it
 * a tunable parameter.
 */
std::vector<AtRandomSample::XYZPoint> AtRandomSample::sampleGaussian(int numPoints, const std::vector<AtHit> &hits)
{
   //--------Gaussian sampling
   double sigma = 30.0;
   double y = 0;
   double gauss = 0;
   int counter = 0;
   int p1 = gRandom->Uniform(0, hits.size());
   int p2;
   auto &P1 = hits.at(p1).GetPosition();

   do {
      p2 = gRandom->Uniform(0, hits.size());
      auto &P2 = hits.at(p2).GetPosition();
      auto dist = std::sqrt((P2 - P1).Mag2());

      gauss = 1.0 * exp(-1.0 * pow(dist / sigma, 2.0));
      y = (gRandom->Uniform(0, 1));
      counter++;
      if (counter > 20 && p2 != p1) {
         LOG(info) << "Gaussian sampling defaulting to uniform";
         break;
      }
   } while (p2 == p1 || y > gauss);

   return {hits.at(p1).GetPosition(), hits.at(p2).GetPosition()};
}

/**
 * @ brief Sample two points based on charge
 */
std::vector<AtRandomSample::XYZPoint> AtRandomSample::sampleWeighted(int numPoints, const std::vector<AtHit> &hits)
{
   //-------Weighted sampling
   double avgCharge = 0;
   auto Proba = getPDF(avgCharge, hits);
   bool validSecondPoint = false;
   int counter = 0;
   int p1 = gRandom->Uniform(0, hits.size());
   int p2 = p1;

   do {
      validSecondPoint = false;
      counter++;
      if (counter > 30 && p2 != p1)
         break;

      p2 = gRandom->Uniform(0, hits.size());
      double TwiceAvCharge = 2 * avgCharge;

      if (Proba.size() == hits.size())
         validSecondPoint = (Proba[p2] >= gRandom->Uniform(0, TwiceAvCharge));
      else
         validSecondPoint = true;

   } while (p2 == p1 || validSecondPoint == false);

   return {hits.at(p1).GetPosition(), hits.at(p2).GetPosition()};
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
std::vector<AtRandomSample::XYZPoint> AtRandomSample::sampleWeightedGaussian(int N, const std::vector<AtHit> &hits)
{
   //-------Weighted sampling + Gauss dist.
   double avgCharge = 0;
   auto Proba = getPDF(avgCharge, hits);
   bool validSecondPoint = false;
   double sigma = 30.0;
   double y = 0;
   double gauss = 0;
   int counter = 0;

   int p1 = gRandom->Uniform(0, hits.size());
   int p2 = p1;
   auto &P1 = hits.at(p1).GetPosition();

   do {
      p2 = gRandom->Uniform(0, hits.size());
      auto &P2 = hits.at(p2).GetPosition();
      auto dist = std::sqrt((P2 - P1).Mag2());

      gauss = 1.0 * exp(-1.0 * pow(dist / sigma, 2));
      y = (gRandom->Uniform(0, 1));

      counter++;
      if (counter > 30 && p2 != p1)
         break;

      validSecondPoint = false;
      double TwiceAvCharge = 2 * avgCharge;

      if (Proba.size() == hits.size())
         validSecondPoint = (Proba[p2] >= gRandom->Uniform(0, TwiceAvCharge));
      else
         validSecondPoint = true;

   } while (p2 == p1 || validSecondPoint == false || y > gauss);

   return {hits.at(p1).GetPosition(), hits.at(p2).GetPosition()};
}

std::vector<double> AtRandomSample::getPDF(double &avgCharge, const std::vector<AtHit> &hits)
{
   double Tcharge =
      std::accumulate(hits.begin(), hits.end(), 0, [](double sum, auto &a) { return sum + a.GetCharge(); });

   avgCharge = Tcharge / hits.size();
   std::vector<double> w;
   if (Tcharge > 0)
      for (const auto &hit : hits)
         w.push_back(hit.GetCharge() / Tcharge);

   return w;
}
