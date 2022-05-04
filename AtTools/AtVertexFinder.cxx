#include "AtVertexFinder.h"

#include <TVector3.h>

XYZPoint AtVertexFinder::FindVertex(const std::vector<AtTrack> &tracks)
{
   XYZPoint fVertex;
   // Assumes the minimum distance between two lines, with respect a given threshold, the first vertex candidate. Then
   // evaluates the distance of each remaining line with respect to the others (vertex) to decide the particles of the
   // reaction.
   // std::cout<<" New find vertex call "<<std::endl;

   Double_t mad = 10; // Minimum approach distance. This is the minimum distance between the lines in 3D. Must be bigger
                      // than fLineDistThreshold
   // ROOT::Math::XYZVector c_1(-1000,-1000,-1000);
   // ROOT::Math::XYZVector c_2(-1000,-1000,-1000);
   TVector3 c_1(-1000, -1000, -1000);
   TVector3 c_2(-1000, -1000, -1000);
   // std::vector<AtTrack*> *TrackCand;

   // Current  parametrization
   // x = p[0] + p[1]*t;
   // y = p[2] + p[3]*t;
   // z = p[4] + p[5]*t;
   //  (x,y,z) = (p[0],p[2],p[4]) + (p[1],p[3],p[5])*t

   // Vector of the beam determined from the experimental data
   TVector3 BeamDir(0., 0., 1.0);
   TVector3 BeamPoint(0., 0., 500.0);

   std::vector<Bool_t> IsFilled(tracks.size(), false);

   // Test each line against the others to find a vertex candidate
   for (Int_t i = 0; i < int(tracks.size()) - 1; i++) {

      const AtTrack &track = tracks.at(i);
      std::vector<Double_t> p = track.GetFitPar();

      if (p.size() == 0)
         continue;

      TVector3 p1(p[0], p[1], p[2]); // p1
      TVector3 e1(p[3], p[4], p[5]); // d1

      // Loop through every other track
      for (Int_t j = i + 1; j < tracks.size(); j++) {
         const AtTrack &track_f = tracks.at(j);
         std::vector<Double_t> p_f = track_f.GetFitPar();
         if (track.GetTrackID() == track_f.GetTrackID())
            continue;

         if (p_f.size() == 0)
            continue;

         TVector3 p2(p_f[0], p_f[1], p_f[2]); // p2
         TVector3 e2(p_f[3], p_f[4], p_f[5]); // d2
         double angle = e1.Angle(e2) * 180. / 3.1415;

         TVector3 n = e1.Cross(e2);
         double sdist = fabs(n.Dot(p1 - p2) / n.Mag());
         TVector3 vertexbuff = ClosestPoint2Lines(e1, p1, e2, p2);

         if (sdist < fLineDistThreshold) {
            TVector3 meanVer = vertexbuff;
            double radius = sqrt(vertexbuff.X() * vertexbuff.X() + vertexbuff.Y() * vertexbuff.Y());

            if (vertexbuff.Z() > 0 && vertexbuff.Z() < 1000 && radius < 25.0) {

               fVertex.SetXYZ(vertexbuff.X(), vertexbuff.Y(), vertexbuff.Z());

               // condition for almost parallel tracks and vertex out of the chamber
            } else if (angle < 10 && angle > 170) {

               TVector3 vertexbuff1 = ClosestPoint2Lines(e1, p1, BeamDir, BeamPoint);
               TVector3 vertexbuff2 = ClosestPoint2Lines(e2, p2, BeamDir, BeamPoint);
               TVector3 vertexbuffav = 0.5 * (vertexbuff1 + vertexbuff2);
               TVector3 vertexbuffdif = vertexbuff1 - vertexbuff2;

               if (vertexbuffdif.Mag() > 2 * fLineDistThreshold)
                  continue;

               fVertex.SetXYZ(vertexbuffav.X(), vertexbuffav.Y(), vertexbuffav.Z());
            }

         } // if fLineDistThreshold

      } // End of track_f

   } // Loop over the tracks
   return fVertex;
}

std::vector<AtTrack> AtVertexFinder::FindVertexOnePerTrack(const std::vector<AtTrack> &tracks)
{
   std::vector<AtTrack> fTrackCand;

   // Assumes the minimum distance between two lines, with respect a given threshold, the first vertex candidate. Then
   // evaluates the distance of each remaining line with respect to the others (vertex) to decide the particles of the
   // reaction.
   // std::cout<<" New find vertex call "<<std::endl;

   // Minimum approach distance. This is the minimum distance between the lines in 3D. Must be bigger
   // than fLineDistThreshold

   Double_t mad = 10;

   TVector3 c_1(-1000, -1000, -1000);
   TVector3 c_2(-1000, -1000, -1000);

   // Current  parametrization
   // x = p[0] + p[1]*t;
   // y = p[2] + p[3]*t;
   // z = p[4] + p[5]*t;
   //  (x,y,z) = (p[0],p[2],p[4]) + (p[1],p[3],p[5])*t

   // Vector of the beam determined from the experimental data
   TVector3 BeamDir(0., 0., 1.0);
   TVector3 BeamPoint(0., 0., 500.0);

   std::vector<Bool_t> IsFilled(tracks.size(), false);

   // Test each line against the others to find a vertex candidate
   for (auto track : tracks) {

      std::vector<Double_t> p = track.GetFitPar();

      if (p.size() == 0)
         continue;

      TVector3 p1(p[0], p[1], p[2]); // p1
      TVector3 e1(p[3], p[4], p[5]); // d1

      TVector3 n = e1.Cross(BeamDir);
      double sdist = fabs(n.Dot(p1 - BeamPoint) / n.Mag());

      TVector3 vertexbuff = ClosestPoint2Lines(e1, p1, BeamDir, BeamPoint);

      if (sdist < fLineDistThreshold) {
         TVector3 meanVer = vertexbuff;
         double radius = sqrt(vertexbuff.X() * vertexbuff.X() + vertexbuff.Y() * vertexbuff.Y());

         // if(vertexbuff.Z()>0 && vertexbuff.Z()<1000 && radius<25.0){
         if (radius < 25.0) {
            track.SetTrackVertex(XYZPoint(vertexbuff));

            if (!IsFilled.at(track.GetTrackID())) {
               IsFilled[track.GetTrackID()] = true;
               fTrackCand.push_back(track);
            }

            // condition for almost parallel tracks and vertex out of the chamber
         }

      } // if fLineDistThreshold

   } // Loop over the tracks

   if (fTrackCand.size() > 5)
      fTrackCand.resize(5); // Truncate the maximum number of lines to 5
   return fTrackCand;
}

TVector3 AtVertexFinder::ClosestPoint2Lines(TVector3 d1, TVector3 pt1, TVector3 d2, TVector3 pt2)
{
   TVector3 n1 = d1.Cross(d2.Cross(d1));
   TVector3 n2 = d2.Cross(d1.Cross(d2));
   double t1 = (pt2 - pt1).Dot(n2) / (d1.Dot(n2));
   double t2 = (pt1 - pt2).Dot(n1) / (d2.Dot(n1));
   TVector3 c1 = pt1 + t1 * d1;
   TVector3 c2 = pt2 + t2 * d2;
   TVector3 meanpoint = 0.5 * (c1 + c2);

   return meanpoint;
}
