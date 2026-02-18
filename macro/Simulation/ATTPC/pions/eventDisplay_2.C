void eventDisplay_2()
{
    FairRunAna* fRun = new FairRunAna();

    FairRuntimeDb* rtdb = fRun->GetRuntimeDb();
    FairParRootFileIo* parIo1 = new FairParRootFileIo();
    parIo1->open("./data/attpcpar_Bfield.root");
    rtdb->setFirstInput(parIo1);
    rtdb->print();

    // fRun->SetInputFile("sim.root");
    // fRun->SetOutputFile("test.root");

    FairRootFileSink* sink = new FairRootFileSink("./data/attpcpar_Bfield.root");
    FairFileSource* source = new FairFileSource("./data/attpctest_Bfield.root");
    fRun->SetSource(source);
    fRun->SetSink(sink);
    // fRun->SetGeomFile(GeoDataPath);

    FairEventManager* fMan = new FairEventManager();
    FairMCTracksDraw* Track = new FairMCTracksDraw();
    FairMCPointDraw* GTPCPoints = new FairMCPointDraw("GTPCPoint", kOrange, kFullSquare);

    fMan->AddTask(Track);
    fMan->AddTask(GTPCPoints);

    fMan->Init();
}
