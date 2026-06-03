void run_sim(Int_t nEvents = 10, TString mcEngine = "TGeant4")
{

   TString dir = getenv("VMCWORKDIR");

   // Output file name
   TString outFile = "./data/attpcsim_Bfield.root";

   // Parameter file name
   TString parFile = "./data/attpcpar_Bfield.root";

   // -----   Timer   --------------------------------------------------------
   TStopwatch timer;
   timer.Start();
   // ------------------------------------------------------------------------

   // gSystem->Load("libAtGen.so");

   // -----   Create simulation run   ----------------------------------------
   FairRunSim *run = new FairRunSim();
   run->SetName(mcEngine);      // Transport engine
   run->SetOutputFile(outFile); // Output file
   FairRuntimeDb *rtdb = run->GetRuntimeDb();
   // ------------------------------------------------------------------------

   // -----   Create media   -------------------------------------------------
   run->SetMaterials("media.geo"); // Materials
   // ------------------------------------------------------------------------

   // -----   Create geometry   ----------------------------------------------

   FairModule *cave = new AtCave("CAVE");
   cave->SetGeometryFileName("cave.geo");
   run->AddModule(cave);

   // FairModule* magnet = new AtMagnet("Magnet");
   // run->AddModule(magnet);

   /*FairModule* pipe = new AtPipe("Pipe");
   run->AddModule(pipe);*/

   FairDetector *ATTPC = new AtTpc("ATTPC", kTRUE);
   ATTPC->SetGeometryFileName("/home/georgina/fair_install/ATTPCROOTv2_KF/geometry/ATTPC_H300torr.root");
   // ATTPC->SetModifyGeometry(kTRUE);
   run->AddModule(ATTPC);

   // ------------------------------------------------------------------------

   // -----   Magnetic field   -------------------------------------------
   // Constant Field
   AtConstField *fMagField = new AtConstField();
   fMagField->SetField(0., 0., 10.);                      // values are in kG
   fMagField->SetFieldRegion(-50, 50, -50, 50, -10, 230); // values are in cm
                                                          //  (xmin,xmax,ymin,ymax,zmin,zmax)
   // run->SetField(fMagField);
   //  --------------------------------------------------------------------

   // -----   Create PrimaryGenerator   --------------------------------------

   Int_t pdgCode = 2212;   // Proton
   Int_t multiplicity = 1; // 1 proton por evento para ver las tracks claras

   /*
   FairBoxGenerator* boxGen = new FairBoxGenerator(pdgCode, multiplicity);

   boxGen->SetPRange(4e6, 4e6);    // Momentum (GeV/c)
   boxGen->SetPhiRange(0., 360.);   // Cobertura azimutal completa
   boxGen->SetThetaRange(0, 90.); // Cobertura polar completa
   boxGen->SetXYZ(0., 0., 0.);     // Origen (centro, desplazado ligeramente en Z)

   FairPrimaryGenerator* primGen = new FairPrimaryGenerator();
   primGen->AddGenerator(boxGen);
   run->SetGenerator(primGen);
   */

   FairPrimaryGenerator *primGen = new FairPrimaryGenerator();
   auto ionGen =
      new FairIonGenerator(9, 25, 9, 1, 0., 0., 10, 0., 0.,
                           -50.); // Z, A, Q, Multiplicity, Px, Py, Pz, Excitation Energy, Mass, Nominal Energy
   primGen->AddGenerator(ionGen);
   run->SetGenerator(primGen);

   // ------------------------------------------------------------------------

   //---Store the visualiztion info of the tracks, this make the output file very large!!
   //--- Use it only to display but not for production!
   run->SetStoreTraj(kTRUE);

   // -----   Initialize simulation run   ------------------------------------
   run->Init();
   // ------------------------------------------------------------------------

   // -----   Runtime database   ---------------------------------------------

   Bool_t kParameterMerged = kTRUE;
   FairParRootFileIo *parOut = new FairParRootFileIo(kParameterMerged);
   parOut->open(parFile.Data());
   rtdb->setOutput(parOut);
   rtdb->saveOutput();
   rtdb->closeOutput();
   rtdb->print();
   // ------------------------------------------------------------------------

   // -----   Start run   ----------------------------------------------------
   run->Run(nEvents);

   // You can export your ROOT geometry ot a separate file
   run->CreateGeometryFile("./data/geofile_full.root");
   // ------------------------------------------------------------------------

   // -----   Finish   -------------------------------------------------------
   timer.Stop();
   Double_t rtime = timer.RealTime();
   Double_t ctime = timer.CpuTime();
   cout << endl << endl;
   cout << "Macro finished succesfully." << endl;
   cout << "Output file is " << outFile << endl;
   cout << "Parameter file is " << parFile << endl;
   cout << "Real time " << rtime << " s, CPU time " << ctime << "s" << endl << endl;
   // ------------------------------------------------------------------------
}
