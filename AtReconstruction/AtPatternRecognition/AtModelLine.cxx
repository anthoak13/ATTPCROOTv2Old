#include "AtModelLine.h"

#include <FairLogger.h>

#include <TMath.h>
AtModelLine::AtModelLine() : AtTrackModel(2) {}

Double_t AtModelLine::DistanceToModel(const AtHit &hit)
{
   auto pos = hit.GetPosition();
   auto vec = fPoint - pos;
   auto nD = fDirection.Cross(vec);
   double dist2 = nD.Mag2() / fDirection.Mag2();
   return std::sqrt(dist2);
}

void AtModelLine::ConstructModel(const std::vector<int> &idx)
{
   if (idx.size() != fNumPoints)
      LOG(error) << "Trying to create model with wrong number of points " << idx.size();
   fIndices = idx;
   auto ind1 = idx.at(0);
   auto ind2 = idx.at(1);
   auto P1 = fHitArray->at(ind1).GetPosition();
   auto P2 = fHitArray->at(ind2).GetPosition();

   fPoint = P1;
   fDirection = P2 - P1;
}

Double_t AtModelLine::Fit3D(const std::vector<int> idx, std::vector<double> &fitPar)
{
   //------3D Line Regression
   //----- adapted from: http://fr.scribd.com/doc/31477970/Regressions-et-trajectoires-3D
   int R, C;
   double Q;
   double Xm, Ym, Zm;
   double Xh, Yh, Zh;
   double a, b;
   double Sxx, Sxy, Syy, Sxz, Szz, Syz;
   double theta;
   double K11, K22, K12, K10, K01, K00;
   double c0, c1, c2;
   double p, q, r, dm2;
   double rho, phi;

   Q = Xm = Ym = Zm = 0.;
   double total_charge = 0;
   Sxx = Syy = Szz = Sxy = Sxz = Syz = 0.;

   for (auto i : idx) {
      auto hitQ = fHitArray->at(i).GetCharge();
      auto &pos = fHitArray->at(i).GetPosition();
      Q += hitQ / 10.;
      Xm += pos.X() * hitQ / 10.;
      Ym += pos.Y() * hitQ / 10.;
      Zm += pos.Z() * hitQ / 10.;
      Sxx += pos.X() * pos.X() * hitQ / 10.;
      Syy += pos.Y() * pos.Y() * hitQ / 10.;
      Szz += pos.Z() * pos.Z() * hitQ / 10.;
      Sxy += pos.X() * pos.Y() * hitQ / 10.;
      Sxz += pos.X() * pos.Z() * hitQ / 10.;
      Syz += pos.Y() * pos.Z() * hitQ / 10.;
   }

   Xm /= Q;
   Ym /= Q;
   Zm /= Q;
   Sxx /= Q;
   Syy /= Q;
   Szz /= Q;
   Sxy /= Q;
   Sxz /= Q;
   Syz /= Q;
   Sxx -= (Xm * Xm);
   Syy -= (Ym * Ym);
   Szz -= (Zm * Zm);
   Sxy -= (Xm * Ym);
   Sxz -= (Xm * Zm);
   Syz -= (Ym * Zm);

   theta = 0.5 * atan((2. * Sxy) / (Sxx - Syy));

   K11 = (Syy + Szz) * pow(cos(theta), 2) + (Sxx + Szz) * pow(sin(theta), 2) - 2. * Sxy * cos(theta) * sin(theta);
   K22 = (Syy + Szz) * pow(sin(theta), 2) + (Sxx + Szz) * pow(cos(theta), 2) + 2. * Sxy * cos(theta) * sin(theta);
   // K12 = -Sxy * (pow(cos(theta), 2) - pow(sin(theta), 2)) + (Sxx - Syy) * cos(theta) * sin(theta);
   K10 = Sxz * cos(theta) + Syz * sin(theta);
   K01 = -Sxz * sin(theta) + Syz * cos(theta);
   K00 = Sxx + Syy;

   c2 = -K00 - K11 - K22;
   c1 = K00 * K11 + K00 * K22 + K11 * K22 - K01 * K01 - K10 * K10;
   c0 = K01 * K01 * K11 + K10 * K10 * K22 - K00 * K11 * K22;

   p = c1 - pow(c2, 2) / 3.;
   q = 2. * pow(c2, 3) / 27. - c1 * c2 / 3. + c0;
   r = pow(q / 2., 2) + pow(p, 3) / 27.;

   if (r > 0)
      dm2 = -c2 / 3. + pow(-q / 2. + sqrt(r), 1. / 3.) + pow(-q / 2. - sqrt(r), 1. / 3.);
   else {
      rho = sqrt(-pow(p, 3) / 27.);
      phi = acos(-q / (2. * rho));
      dm2 = std::min(-c2 / 3. + 2. * pow(rho, 1. / 3.) * cos(phi / 3.),
                     std::min(-c2 / 3. + 2. * pow(rho, 1. / 3.) * cos((phi + 2. * TMath::Pi()) / 3.),
                              -c2 / 3. + 2. * pow(rho, 1. / 3.) * cos((phi + 4. * TMath::Pi()) / 3.)));
   }

   a = -K10 * cos(theta) / (K11 - dm2) + K01 * sin(theta) / (K22 - dm2);
   b = -K10 * sin(theta) / (K11 - dm2) - K01 * cos(theta) / (K22 - dm2);

   Xh = ((1. + b * b) * Xm - a * b * Ym + a * Zm) / (1. + a * a + b * b);
   Yh = ((1. + a * a) * Ym - a * b * Xm + b * Zm) / (1. + a * a + b * b);
   Zh = ((a * a + b * b) * Zm + a * Xm + b * Ym) / (1. + a * a + b * b);

   // First 3 are point1. Second 3 are point 2
   fitPar = {Xm, Ym, Zm, Xh, Yh, Zh};
   return (fabs(dm2 / Q));
}
