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
const TString geoVersion = "MAGNEX_in_development";
const TString outputGeoFileName = geoVersion + ".root";
const TString outputGeoManagerFileName = geoVersion + "_geomanager.root";

// Names of the different used materials which are used to build the modules
// The materials are defined in the global media.geo file
const TString MediumVacuum = "vacuum4";
const TString MediumGas = "H_1bar";

TGeoManager *gGeoMan = new TGeoManager("MAGNEX", "MAGNEX");   // Pointer to TGeoManager instance
TGeoVolume *gModules;                                         // Global storage for module types


const Float_t xSize = 30.;   // cm
const Float_t ySize = 18.5;  // cm
const Float_t zSize = 10.87; // cm

void create_materials_from_media_file();
TGeoVolume *create_detector();

void MAGNEX_geometry()
{

   gSystem->Load("libGeoBase");
   gSystem->Load("libParBase");
   gSystem->Load("libBase");

   create_materials_from_media_file();

   // Get the GeoManager for later usage
   gGeoMan = (TGeoManager *)gROOT->FindObject("FAIRGeom");
   gGeoMan->SetVisLevel(7);

   TGeoVolume *top = new TGeoVolumeAssembly("TOP");
   gGeoMan->SetTopVolume(top);

   TGeoMedium *vac = gGeoMan->GetMedium(MediumVacuum);
   TGeoVolume *tpcvac = new TGeoVolumeAssembly(geoVersion);
   tpcvac->SetMedium(vac);
   top->AddNode(tpcvac, 1);

   gModules = create_detector();

   std::cout << "Voxelizing." << std::endl;
   top->Voxelize("");
   gGeoMan->CloseGeometry();

   gGeoMan->CheckOverlaps(0.001);
   gGeoMan->PrintOverlaps();
   gGeoMan->Test();

   TFile *outGeoFile = new TFile(outputGeoFileName, "RECREATE");
   top->Write();
   outGeoFile->Close();

   TFile *outGeoManFile = new TFile(outputGeoManagerFileName, "RECREATE");
   gGeoMan->Write();
   outGeoManFile->Close();

   top->Draw("ogl");
}


void create_materials_from_media_file()
{
   // Use the FairRoot geometry interface to load the media which are already defined
   FairGeoLoader *geoLoad = new FairGeoLoader("TGeo", "FairGeoLoader");
   FairGeoInterface *geoFace = geoLoad->getGeoInterface();
   TString attpcrootPath = gSystem->Getenv("VMCWORKDIR");
   TString geoFile = attpcrootPath + "/geometry/media.geo";
   geoFace->setMediaFile(geoFile);
   geoFace->readMedia();

   // Read the required media and create them in the GeoManager
   FairGeoMedia *geoMedia = geoFace->getMedia();
   FairGeoBuilder *geoBuild = geoLoad->getGeoBuilder();

   FairGeoMedium *vacuum = geoMedia->getMedium(MediumVacuum);
   FairGeoMedium *gas = geoMedia->getMedium(MediumGas);

   // Build the required media

   geoBuild->createMedium(vacuum);
   geoBuild->createMedium(gas);
}

TGeoVolume *create_detector()
{
   TGeoMedium *gas = gGeoMan->GetMedium(MediumGas);

   TGeoVolume *drift_volume = gGeoManager->MakeBox("drift_volume", gas, xSize / 2., ySize / 2., zSize / 2.);
   gGeoMan->GetVolume(geoVersion)->AddNode(drift_volume, 1, new TGeoTranslation(xSize / 2., ySize / 2., zSize / 2.));
   drift_volume->SetTransparency(80);

   return drift_volume;
}

