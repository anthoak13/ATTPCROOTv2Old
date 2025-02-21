/********************************************************************************
 *    Copyright (C) 2014 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH    *
 *                                                                              *
 *              This software is distributed under the terms of the             *
 *         GNU Lesser General Public Licence version 3 (LGPL) version 3,        *
 *                  copied verbatim in the file "LICENSE"                       *
 ********************************************************************************/
// in root all sizes are given in cm

#include "TFile.h"
#include "TGeoCompositeShape.h"
#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMatrix.h"
#include "TGeoMedium.h"
#include "TGeoPgon.h"
#include "TGeoVolume.h"
#include "TList.h"
#include "TROOT.h"
#include "TString.h"
#include "TSystem.h"

#include <iostream>

// Name of geometry version and output file
const TString geoVersion = "AT3PC_3H_5torr";
const TString FileName = geoVersion + ".root";
const TString FileName1 = geoVersion + "_geomanager.root";

// Names of the different used materials which are used to build the modules
// The materials are defined in the global media.geo file
const TString MediumGas = "ATTPCIsoButane_50";
const TString CylinderVolumeMedium = "steel";
const TString MediumVacuum = "vacuum4";
const TString CellMedium = "3H_5torr";
const TString CellMaterial = "mylar";

// Distance of the center of the first detector layer [cm];
const Float_t First_Z_Position = 10;
const Float_t Z_Distance = 10;

const Float_t tpc_diameter = 50.;
const Float_t drift_length = 100.;
const Float_t cell_thickness = 0.0005;

// some global variables
TGeoManager *gGeoMan = new TGeoManager("ATTPC", "ATTPC");
;                     // Pointer to TGeoManager instance
TGeoVolume *gModules; // Global storage for module types

// Forward declarations
void create_materials_from_media_file();
TGeoVolume *create_detector();
void position_detector();
void add_alignable_volumes();

void AT3PC_3H_5torr()
{
   // Load the necessary FairRoot libraries
   // gROOT->LoadMacro("$VMCWORKDIR/gconfig/basiclibs.C");
   // basiclibs();
   // gSystem->Load("libGeoBase");
   // gSystem->Load("libParBase");
   // gSystem->Load("libBase");

   // Load needed material definition from media.geo file
   create_materials_from_media_file();

   // Get the GeoManager for later usage
   gGeoMan = (TGeoManager *)gROOT->FindObject("FAIRGeom");
   gGeoMan->SetVisLevel(7);

   // Create the top volume

   TGeoVolume *top = new TGeoVolumeAssembly("TOP");
   gGeoMan->SetTopVolume(top);

   TGeoMedium *gas = gGeoMan->GetMedium(MediumVacuum);
   TGeoVolume *tpcvac = new TGeoVolumeAssembly(geoVersion);
   tpcvac->SetMedium(gas);
   top->AddNode(tpcvac, 1);

   gModules = create_detector();

   // position_detector();

   cout << "Voxelizing." << endl;
   top->Voxelize("");
   gGeoMan->CloseGeometry();

   // add_alignable_volumes();

   gGeoMan->CheckOverlaps(0.001);
   gGeoMan->PrintOverlaps();
   gGeoMan->Test();

   TFile *outfile = new TFile(FileName, "RECREATE");
   top->Write();
   outfile->Close();

   TFile *outfile1 = new TFile(FileName1, "RECREATE");
   gGeoMan->Write();
   outfile1->Close();

   top->Draw("ogl");
   // top->Raytrace();
}

void create_materials_from_media_file()
{
   // Use the FairRoot geometry interface to load the media which are already defined
   FairGeoLoader *geoLoad = new FairGeoLoader("TGeo", "FairGeoLoader");
   FairGeoInterface *geoFace = geoLoad->getGeoInterface();
   TString geoPath = gSystem->Getenv("VMCWORKDIR");
   TString geoFile = geoPath + "/geometry/media.geo";
   geoFace->setMediaFile(geoFile);
   geoFace->readMedia();

   // Read the required media and create them in the GeoManager
   FairGeoMedia *geoMedia = geoFace->getMedia();
   FairGeoBuilder *geoBuild = geoLoad->getGeoBuilder();

   FairGeoMedium *isobutane = geoMedia->getMedium("ATTPCIsoButane_50");
   FairGeoMedium *steel = geoMedia->getMedium("steel");
   FairGeoMedium *vacuum4 = geoMedia->getMedium("vacuum4");
   FairGeoMedium *mylar = geoMedia->getMedium("mylar");
   FairGeoMedium *tritium = geoMedia->getMedium("3H_5torr");

   // include check if all media are found

   geoBuild->createMedium(isobutane);
   geoBuild->createMedium(steel);
   geoBuild->createMedium(vacuum4);
   geoBuild->createMedium(mylar);
   geoBuild->createMedium(tritium);
}

TGeoVolume *create_detector()
{

   // needed materials
   TGeoMedium *OuterCylinder = gGeoMan->GetMedium(CylinderVolumeMedium);
   TGeoMedium *gas = gGeoMan->GetMedium(MediumGas);
   TGeoMedium *mylar = gGeoMan->GetMedium(CellMaterial);
   TGeoMedium *tritiumgas = gGeoMan->GetMedium(CellMedium);

   TGeoVolume *chamber_volume =
      gGeoManager->MakeTube("chamber_volume", OuterCylinder, tpc_diameter / 2, (tpc_diameter + 10.0) / 2, 120.0 / 2);
   gGeoMan->GetVolume(geoVersion)->AddNode(chamber_volume, 1, new TGeoTranslation(0, 0, 120.0 / 2));

   TGeoVolume *drift_volume =
      gGeoManager->MakeTube("drift_volume", gas, (5.0 + cell_thickness) / 2, tpc_diameter / 2, drift_length / 2);
   gGeoMan->GetVolume(geoVersion)->AddNode(drift_volume, 1, new TGeoTranslation(0, 0, drift_length / 2));
   drift_volume->SetTransparency(80);

   TGeoVolume *mylar_cell =
      gGeoManager->MakeTube("mylar_cell", mylar, 5.0 / 2, (5.0 + cell_thickness) / 2, drift_length / 2);
   mylar_cell->SetLineColor(kBlue);
   gGeoMan->GetVolume(geoVersion)->AddNode(mylar_cell, 1, new TGeoTranslation(0, 0, drift_length / 2));
   mylar_cell->SetTransparency(80);

   TGeoVolume *cell = gGeoManager->MakeTube("cell", tritiumgas, 0, 5.0 / 2, drift_length / 2);
   cell->SetLineColor(kRed);
   gGeoMan->GetVolume(geoVersion)->AddNode(cell, 1, new TGeoTranslation(0, 0, drift_length / 2));
   cell->SetTransparency(10);

   return drift_volume;
}

void position_detector()
{

   /* TGeoTranslation* det_trans=NULL;

    Int_t numDets=0;
    for (Int_t detectorPlanes = 0; detectorPlanes < 40; detectorPlanes++) {
      det_trans
        = new TGeoTranslation("", 0., 0., First_Z_Position+(numDets*Z_Distance));
      gGeoMan->GetVolume(geoVersion)->AddNode(gModules, numDets, det_trans);
      numDets++;

    }*/
}

void add_alignable_volumes()
{

   /* TString volPath;
    TString symName;
    TString detStr   = "Tutorial4/det";
    TString volStr   = "/TOP_1/tutorial4_1/tut4_det_";

    for (Int_t detectorPlanes = 0; detectorPlanes < 40; detectorPlanes++) {

      volPath  = volStr;
      volPath += detectorPlanes;

      symName  = detStr;
      symName += Form("%02d",detectorPlanes);

      cout<<"Path: "<<volPath<<", "<<symName<<endl;
  //    gGeoMan->cd(volPath);

      gGeoMan->SetAlignableEntry(symName.Data(),volPath.Data());

    }
      cout<<"Nr of alignable objects: "<<gGeoMan->GetNAlignable()<<endl;*/
}
