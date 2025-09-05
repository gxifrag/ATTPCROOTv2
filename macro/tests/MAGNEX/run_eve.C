
void run_eve()
{
   TString InputDataPath = "/home/aurio/research/NUMEN/NUMEN_output/reconstructed_005.root";
   TString OutputDataPath = "/home/aurio/research/NUMEN/NUMEN_output/output.reco_display.root";
   std::cout << "Opening: " << InputDataPath << std::endl;

   TString attpcrootPath = gSystem->Getenv("VMCWORKDIR");
   TString geoFile = "MAGNEX_in_development_geomanager.root";
   TString mapFile = "e12014_pad_mapping.xml";
   TString GeoDataPath = attpcrootPath + "/geometry/" + geoFile;
   TString mapDir = attpcrootPath + "/scripts/" + mapFile;

   FairRunAna *fRun = new FairRunAna();
   FairRootFileSink *sink = new FairRootFileSink(OutputDataPath);
   FairFileSource *source = new FairFileSource(InputDataPath);
   fRun->SetSource(source);
   fRun->SetSink(sink);
   fRun->SetGeomFile(GeoDataPath);

   FairRuntimeDb *rtdb = fRun->GetRuntimeDb();
   FairParRootFileIo *parIo = new FairParRootFileIo();
   rtdb->setFirstInput(parIo);

   TString parFile = "/home/aurio/research/NUMEN/NUMEN_data/testParFile.txt";
   auto fMap = std::make_shared<AtMAGNEXMap>(parFile);
   fMap->ParseXMLMap(mapDir.Data());
   fMap->SetPadPlanePlane("XZ");
   std::cout << "Current PadPlane plane: " << fMap->GetPadPlanePlane() << "." << std::endl;
   AtViewerManager *eveMan = new AtViewerManager(fMap);

   auto tabMain = std::make_unique<AtTabMAGNEX>();
   tabMain->SetMultiHit(100);

   eveMan->AddTab(std::move(tabMain));

   eveMan->Init();

   std::cout << " Finished init. " << std::endl;
}
