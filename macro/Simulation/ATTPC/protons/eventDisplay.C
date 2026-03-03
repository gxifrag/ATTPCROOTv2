void eventDisplay()
{
   //-----User Settings:-----------------------------------------------
   TString InputDataFile = "/home/georgina/fair_install/ATTPCROOTv2_KF/macro/Simulation/ATTPC/protons/data/protonssim_3T_H60torr_60MeV_theta30.root";
   TString ParFile = "/home/georgina/fair_install/ATTPCROOTv2_KF/macro/Simulation/ATTPC/protons/data/protonspar_3T_H60torr_60MeV_theta30.root";
   TString OutputDataFile = "/home/georgina/fair_install/ATTPCROOTv2_KF/macro/Simulation/ATTPC/protons/data/protonsctest_3T_H60torr_60MeV_theta30.root";



   // -----   Reconstruction run   -------------------------------------------

   FairRunAna *fRun = new FairRunAna();
   FairRootFileSink *sink = new FairRootFileSink(OutputDataFile);
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
   //FairMCTracks *Track = new FairMCTracks("Monte-Carlo Tracks");
   FairMCPointDraw *AtTpcPoints = new FairMCPointDraw("AtTpcPoint", kBlue, kFullSquare);

   //fMan->AddTask(Track);
   fMan->AddTask(AtTpcPoints);

   fMan->Init();
}
