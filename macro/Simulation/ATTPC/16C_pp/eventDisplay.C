void eventDisplay()
{
   //-----User Settings:-----------------------------------------------
   TString InputDataFile = "/home/georgina/fair_install/ATTPCROOTv2_KF/macro/Simulation/ATTPC/16C_pp/data/attpcsim_Bfield.root";
   TString ParFile = "/home/georgina/fair_install/ATTPCROOTv2_KF/macro/Simulation/ATTPC/16C_pp/data/attpcpar_Bfield.root";
   TString OutputDataFile = "/home/georgina/fair_install/ATTPCROOTv2_KF/macro/Simulation/ATTPC/16C_pp/data/attpctest_Bfield.root";

   /*TString InputDataFile = "/home/georgina/fair_install/ATTPCROOTv2_KF/macro/Simulation/ATTPC/16C_pp/data/attpcsim_proton_gun.root";
   TString ParFile = "/home/georgina/fair_install/ATTPCROOTv2_KF/macro/Simulation/ATTPC/16C_pp/data/attpcpar_proton_gun.root";  
   TString OutputDataFile = "/home/georgina/fair_install/ATTPCROOTv2_KF/macro/Simulation/ATTPC/16C_pp/data/attpctest_proton_gun.root";*/

/*   TString InputDataFile = "/home/georgina/fair_install/ATTPCROOTv2_KF/macro/Simulation/ATTPC/16C_pp/data/attpcsim.root";
   TString ParFile = "/home/georgina/fair_install/ATTPCROOTv2_KF/macro/Simulation/ATTPC/16C_pp/data/attpcpar.root";
   TString OutputDataFile = "/home/georgina/fair_install/ATTPCROOTv2_KF/macro/Simulation/ATTPC/16C_pp/data/attpctest.root";
*/

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
   // FairMCTracks *Track = new FairMCTracks("Monte-Carlo Tracks");
   FairMCPointDraw *AtTpcPoints = new FairMCPointDraw("AtTpcPoint", kBlue, kFullSquare);

   // fMan->AddTask(Track);
   fMan->AddTask(AtTpcPoints);

   fMan->Init();
}
