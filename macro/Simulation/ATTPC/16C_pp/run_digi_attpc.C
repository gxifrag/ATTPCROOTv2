bool reduceFunc(AtRawEvent *evt);

void run_digi_attpc()
{
   TString inOutDir = "/home/georgina/fair_install/ATTPCROOTv2_KF/macro/Simulation/ATTPC/16C_pp/data/";
   TString outputFile = inOutDir + "output_digi.root";
   TString scriptfile = "Lookup20150611.xml";
   TString paramFile = "ATTPC.e20009_sim.par";

   TString dir = getenv("VMCWORKDIR");

   // TString mcFile = "./data/sim_attpc.root";
   TString mcFile = inOutDir + "attpcsim.root";

   // Create the full parameter file paths
   TString digiParFile = dir + "/parameters/" + paramFile;
   TString mapParFile = dir + "/scripts/" + scriptfile;

   // -----   Timer   --------------------------------------------------------
   TStopwatch timer;

   // ------------------------------------------------------------------------

   // __ Run ____________________________________________
   FairRunAna *fRun = new FairRunAna();
   //leer eventos de un ROOT file
   FairFileSource *source = new FairFileSource(mcFile); 
   fRun->SetSource(source);
   fRun->SetOutputFile(outputFile);

   FairRuntimeDb *rtdb = fRun->GetRuntimeDb();
   FairParAsciiFileIo *parIo1 = new FairParAsciiFileIo();
   parIo1->open(digiParFile.Data(), "in");
   rtdb->setFirstInput(parIo1);

   // Create the detector map to pass to the simulation
   auto mapping = std::make_shared<AtTpcMap>();
   mapping->ParseXMLMap(mapParFile.Data());
   mapping->GeneratePadPlane();
   //mapping->ParseInhibitMap("./data/inhibit.txt", AtMap::InhibitType::kTotal);

   //Clusterize → Pulse → PSA → PRA
   AtClusterizeTask *clusterizer = new AtClusterizeTask(); //agrupa señales vecinas
   clusterizer->SetPersistence(kFALSE);

   AtPulseTask *pulse = new AtPulseTask(std::make_shared<AtPulse>(mapping)); //convierte las señales en pulsos físicos
   pulse->SetPersistence(kTRUE);
   pulse->SetSaveMCInfo();

   auto psa = std::make_unique<AtPSAMax>(); //Extrae posición/tiempo/carga → hits.
   psa->SetThreshold(0);

   // Create PSA task
   AtPSAtask *psaTask = new AtPSAtask(std::move(psa));
   psaTask->SetPersistence(kTRUE);

   AtPRAtask *praTask = new AtPRAtask();
   praTask->SetPersistence(kTRUE); //Reconstrucción de trayectorias.
   

   fRun->AddTask(clusterizer);
   fRun->AddTask(pulse);
   fRun->AddTask(psaTask);
   fRun->AddTask(praTask);

   //  __ Init and run ___________________________________
   //loop de eventos
   fRun->Init();

   timer.Start();
   fRun->Run(0, 1000);
   timer.Stop();


   std::cout << std::endl << std::endl;
   std::cout << "Macro finished succesfully." << std::endl << std::endl;
   // -----   Finish   -------------------------------------------------------

   Double_t rtime = timer.RealTime();
   Double_t ctime = timer.CpuTime();
   cout << endl;
   cout << "Real time " << rtime << " s, CPU time " << ctime << " s" << endl;
   cout << endl;
   // ------------------------------------------------------------------------
}

bool reduceFunc(AtRawEvent *evt) //se evalua justo antes de ejecutar las tasks
{
   return (evt->GetNumPads() > 0) && evt->IsGood();
}
