void run_sim_protons(Int_t nEvents = 10000, TString mcEngine = "TGeant4")
{

   TString dir = getenv("VMCWORKDIR");

   // Output file name
   TString outFile = "./data/protonssim_3T_H300torr_40MeV_theta30.root";


   // Parameter file name
   TString parFile = "./data/protonspar_3T_H300torr_40MeV_theta30.root";

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
   //ATTPC->SetModifyGeometry(kTRUE);
   run->AddModule(ATTPC);

   // -----   Magnetic field   -------------------------------------------
   // Constant Field
   AtConstField *fMagField = new AtConstField();
   fMagField->SetField(0., 0., 30.);                     // values are in 20 kG
   fMagField->SetFieldRegion(-50, 50, -50, 50, -10, 230); // values are in cm
                                                          //  (xmin,xmax,ymin,ymax,zmin,zmax)
   run->SetField(fMagField);
   
   // -----   Create PrimaryGenerator   --------------------------------------

   Int_t pdgCode = 2212; // Proton
   //Int_t pdgCode = 211; // Pion
   Int_t multiplicity = 1; // 1 proton por evento para ver las tracks claras
   Double_t pMin = 0.04; // Momentum mínimo en GeV/c
   Double_t pMax = 0.04; // Momentum máximo en GeV/c

   FairBoxGenerator* boxGen = new FairBoxGenerator(pdgCode, multiplicity);
   
   boxGen->SetPRange(pMin, pMax);    // Momentum (GeV/c)
   boxGen->SetPhiRange(0., 360.);   // Cobertura azimutal completa
   boxGen->SetThetaRange(30., 30.); // Cobertura polar completa
   boxGen->SetXYZ(0., 0., 0.);     // Origen (centro, desplazado ligeramente en Z)

   FairPrimaryGenerator* primGen = new FairPrimaryGenerator();
   primGen->AddGenerator(boxGen);
   run->SetGenerator(primGen);

   // ====================================================================
   // FIX: Manually initialize AtVertexPropagator for the AtTpc detector
   // ====================================================================
   AtVertexPropagator* vertexProp = AtVertexPropagator::Instance();
   
   // Set the mass of the beam particle (Proton mass is ~1.0078 u)
   Double_t massAmu = 1.007825; 
   Double_t massGeV = 0.938272;

   vertexProp->SetBeamMass(massAmu); 
   //vertexProp->SetBeamMass(0.13957); // Pion mass in GeV/c^2
   
   // Set Nominal Energy (matching your BoxGen momentum)
   Double_t pAvg = (pMin + pMax) / 2.0; 
   Double_t eTotal = TMath::Sqrt(pAvg*pAvg + massGeV*massGeV); // E = sqrt(p^2 + m^2)
   Double_t nominalEnergy = eTotal - massGeV;
   vertexProp->SetBeamNomE(nominalEnergy);
   
   // Prevent the "reactionOccursHere" logic from triggering randomly 
   // by setting a very high energy loss threshold, effectively treating 
   // this as a pure tracking event, not a beam-reaction event.
   vertexProp->SetRndELoss(1e9); 
   // ====================================================================


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
