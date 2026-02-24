void pion_box_sim(Int_t nEvents = 100, TString mcEngine = "TGeant4")
{
   // 1. Definición de Archivos
   TString dir = getenv("VMCWORKDIR");
   TString outFile = "./output_pions_box.root";
   TString parFile = "./params_pions_box.root";

   // 2. Timer
   TStopwatch timer;
   timer.Start();

   // 3. Cargar librerías (Por seguridad)
   gSystem->Load("libFairRoot");
   gSystem->Load("libAtTpcRoot");

   // 4. Run Manager
   FairRunSim *run = new FairRunSim();
   run->SetName(mcEngine);      
   run->SetOutputFile(outFile); 
   FairRuntimeDb *rtdb = run->GetRuntimeDb();

   // 5. Materiales
   run->SetMaterials("media.geo"); 

   // 6. Geometría (USANDO TU RUTA EXACTA QUE FUNCIONA)
   FairModule *cave = new AtCave("CAVE");
   cave->SetGeometryFileName("cave.geo");
   run->AddModule(cave);

   FairDetector *ATTPC = new AtTpc("ATTPC", kTRUE);
   // ¡Esta es la línea clave que hace que funcione!
   ATTPC->SetGeometryFileName("/home/georgina/fair_install/ATTPCROOTv2_KF/geometry/ATTPC_H300torr.root");
   run->AddModule(ATTPC);

   // 7. Campo Magnético
   AtConstField *fMagField = new AtConstField();
   fMagField->SetField(0., 0., 20.0);                     // 20 kG (2 Tesla) para piones
   fMagField->SetFieldRegion(-50, 50, -50, 50, -10, 230); // Región en cm
   run->SetField(fMagField);

   // 8. GENERADOR DE PIONES (Aquí cambiamos la física)
   FairPrimaryGenerator *primGen = new FairPrimaryGenerator();

   // Configuración del Box Generator
   // PDG -211 = Pi- (Negativo)
   // Multiplicidad = 1 (1 pión por evento)
   FairBoxGenerator* boxGen = new FairBoxGenerator(-211, 1);
   
   // Rangos para visualizar espirales
   boxGen->SetPRange(0.05, 0.2);    // Momento: 50 a 200 MeV/c 
   boxGen->SetPhiRange(0., 360.);   // 360 grados
   boxGen->SetThetaRange(0., 90.);  // Hacia adelante
   boxGen->SetXYZ(0., 0., 10.);     // Salen desde z=10cm

   primGen->AddGenerator(boxGen);
   run->SetGenerator(primGen);

   // 9. Inicializar y Correr
   // Guardar trayectorias para verlas en visualización
   run->SetStoreTraj(kTRUE);
   
   run->Init();

   // Guardar Parámetros
   Bool_t kParameterMerged = kTRUE;
   FairParRootFileIo *parOut = new FairParRootFileIo(kParameterMerged);
   parOut->open(parFile.Data());
   rtdb->setOutput(parOut);
   rtdb->saveOutput();
   rtdb->print();

   // Correr eventos
   run->Run(nEvents);

   // Finalizar
   timer.Stop();
   Double_t rtime = timer.RealTime();
   Double_t ctime = timer.CpuTime();
   cout << endl << endl;
   cout << ">>> Simulacion de Piones terminada con exito." << endl;
   cout << ">>> Output: " << outFile << endl;
   cout << ">>> Real time " << rtime << " s, CPU time " << ctime << "s" << endl << endl;
}
