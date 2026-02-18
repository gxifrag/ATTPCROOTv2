void proton_particle_gun(Int_t nEvents = 100, TString mcEngine = "TGeant4")
{
   TString outFile = "./data/attpcsim_proton_gun.root";
   TString parFile = "./data/attpcpar_proton_gun.root";

      // -----   Timer   --------------------------------------------------------
   TStopwatch timer;
   timer.Start();

   FairRunSim *run = new FairRunSim();
   run->SetName("TGeant4");      
   run->SetOutputFile(outFile); 
   FairRuntimeDb *rtdb = run->GetRuntimeDb();

   run->SetMaterials("media.geo"); 

   FairModule *cave = new AtCave("CAVE");
   cave->SetGeometryFileName("cave.geo");
   run->AddModule(cave);

   FairDetector *ATTPC = new AtTpc("ATTPC", kTRUE);
   // USAMOS LA RUTA ABSOLUTA QUE SÍ TE FUNCIONA
   ATTPC->SetGeometryFileName("/home/georgina/fair_install/ATTPCROOTv2_KF/geometry/ATTPC_H300torr.root");
   run->AddModule(ATTPC);

   AtConstField *fMagField = new AtConstField();
   fMagField->SetField(0., 0., 28.5); // 2.85 Tesla
   fMagField->SetFieldRegion(-50, 50, -50, 50, -10, 230); 
   run->SetField(fMagField);

   FairPrimaryGenerator *primGen = new FairPrimaryGenerator();

  // Definición del Protón (Z=1, A=1)
   Int_t z = 1;  // Número Atómico
   Int_t a = 1;  // Número de Masa
   Int_t q = 1;  // Estado de carga (1 para protón)
   Int_t m = 1;  // Multiplicidad (cuántos protones por evento)

   // Momento (GeV/c)
   Double_t px = 0.000; 
   Double_t py = 0.000; 
   Double_t pz = 0.400; // Por ejemplo, 400 MeV/c
   
   Double_t BExcEner = 0.0;       // Energía de excitación
   Double_t Bmass    = 1.007825;  // Masa en uma
   Double_t NomEnergy = 10.0;     // Energía nominal en MeV/u (opcional para el cálculo interno)

   /*AtTPCIonGenerator *ionGen = new AtTPCIonGenerator("Ion", z, a, q, m, px, py, pz, BExcEner, Bmass, NomEnergy);   
   primGen->AddGenerator(ionGen);
   run->SetGenerator(primGen);*/


   FairBoxGenerator* boxGen = new FairBoxGenerator(2212, 1);
   boxGen->SetXYZ(0, 0, 10);
   boxGen->SetPRange(0.4, 0.6);
   boxGen->SetThetaRange(0, 5);
   primGen->AddGenerator(boxGen);
   run->SetGenerator(primGen);

   run->Init();

   // -----   Runtime database   ---------------------------------------------

   Bool_t kParameterMerged = kTRUE;
   FairParRootFileIo *parOut = new FairParRootFileIo(kParameterMerged);
   parOut->open(parFile.Data());
   rtdb->setOutput(parOut);
   rtdb->saveOutput();
   //rtdb->print();
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
