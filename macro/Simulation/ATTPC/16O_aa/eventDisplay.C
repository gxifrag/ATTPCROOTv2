void eventDisplay()
{
  //-----User Settings:-----------------------------------------------
  TString InputDataFile = "/home/georgina/fair_install/ATTPCROOTv2_KF/macro/Simulation/ATTPC/16O_aa/data/attpcsim.root";
  TString ParFile = "/home/georgina/fair_install/ATTPCROOTv2_KF/macro/Simulation/ATTPC/16O_aa/data/attpcpar.root";
  TString OutFile = "/home/georgina/fair_install/ATTPCROOTv2_KF/macro/Simulation/ATTPC/16O_aa/data/attpctest.root";

  // -----   Reconstruction run   -------------------------------------------
  
   FairRunAna *fRun = new FairRunAna();
   FairRootFileSink *sink = new FairRootFileSink(OutFile);
   FairFileSource *source = new FairFileSource(InputDataFile);
   fRun->SetSource(source);
   fRun->SetSink(sink);
   // fRun->SetGeomFile(GeoDataPath);

   FairRuntimeDb *rtdb = fRun->GetRuntimeDb();
   FairParRootFileIo *parInput1 = new FairParRootFileIo();
   parInput1->open(ParFile.Data());
   rtdb->setFirstInput(parInput1);

   FairEventManager *fMan = new FairEventManager();

   //----------------------Traks and points -------------------------------------
   // FairMCTracks *Track = new FairMCTracks("Monte-Carlo Tracks");
   FairMCPointDraw *AtTpcPoints = new FairMCPointDraw("AtTpcPoint", kBlue, kFullSquare);

   // fMan->AddTask(Track);
   fMan->AddTask(AtTpcPoints);

   fMan->Init();
}


