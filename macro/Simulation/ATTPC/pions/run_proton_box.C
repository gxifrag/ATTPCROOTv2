void run_proton_box(Int_t nEvents = 1000)
{
    TString outFile = "proton_box.root";
    TString parFile = "proton_box_params.root";

    // Timer
    TStopwatch timer;
    timer.Start();

    // ------------------------------------------------------------------------
    // FairRunSim
    // ------------------------------------------------------------------------
    // -----   Create simulation run   ----------------------------------------
   FairRunSim *run = new FairRunSim();
   run->SetName(mcEngine);      // Transport engine
   run->SetOutputFile(outFile); // Output file
   FairRuntimeDb *rtdb = run->GetRuntimeDb();
   // ------------------------------------------------------------------------

   // -----   Create media   -------------------------------------------------
   run->SetMaterials("media.geo"); // Materials
    // ------------------------------------------------------------------------
    // Geometría (ejemplo vacío)
    // ------------------------------------------------------------------------
    FairModule *cave = new FairCave("CAVE");
    cave->SetGeometryFileName("cave.geo");
    run->AddModule(cave);

    // ------------------------------------------------------------------------
    // Campo magnético uniforme
    // ------------------------------------------------------------------------
    Double_t Bx = 0.0;   // Tesla
    Double_t By = 0.0;
    Double_t Bz = 1.0;   // Campo en z

    FairConstField *magField = new FairConstField();
    magField->SetField(Bx, By, Bz);

    // Región donde actúa el campo (cm)
    magField->SetFieldRegion(-100, 100,
                             -100, 100,
                             -100, 100);

    run->SetField(magField);

    // ------------------------------------------------------------------------
    // Generador primario
    // ------------------------------------------------------------------------
    FairPrimaryGenerator *primGen = new FairPrimaryGenerator();

    // PDG: protón = 2212
    Int_t pdgId = 2212;

    // Generador tipo caja
    FairBoxGenerator *boxGen = new FairBoxGenerator(pdgId, 1);

    // Distribución en posición (cm)
    boxGen->SetXYZ(0.0, 0.0, 0.0);

    // Rango de ángulos
    boxGen->SetThetaRange(0., 10.);   // grados
    boxGen->SetPhiRange(0., 360.);    // grados

    // Energía o momento
    boxGen->SetPRange(0.5, 1.0);  // GeV/c

    primGen->AddGenerator(boxGen);
    run->SetGenerator(primGen);

    // ------------------------------------------------------------------------
    // Inicialización
    // ------------------------------------------------------------------------
    run->Init();

    // ------------------------------------------------------------------------
    // Parámetros de salida
    // ------------------------------------------------------------------------
    FairRuntimeDb *rtdb = run->GetRuntimeDb();
    FairParRootFileIo *parOut = new FairParRootFileIo(kTRUE);
    parOut->open(parFile);
    rtdb->setOutput(parOut);
    rtdb->saveOutput();

    // ------------------------------------------------------------------------
    // Simulación
    // ------------------------------------------------------------------------
    run->Run(nEvents);

    // ------------------------------------------------------------------------
    // Final
    // ------------------------------------------------------------------------
    timer.Stop();
    std::cout << "Tiempo real: " << timer.RealTime() << " s" << std::endl;
    std::cout << "Tiempo CPU:  " << timer.CpuTime()  << " s" << std::endl;
}

