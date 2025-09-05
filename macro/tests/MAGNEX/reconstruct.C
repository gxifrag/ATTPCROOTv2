void reconstruct(int runNumber = 210)
{
   // Load the library for unpacking and reconstruction
   gSystem->Load("libAtReconstruction.so");

   FairLogger *logger = FairLogger::GetLogger();
   logger->SetLogVerbosityLevel("high");

   TStopwatch timer;
   timer.Start();

   // Set the input/output directories
   TString inputDir = "/home/aurio/research/NUMEN/NUMEN_data/";
   TString outDir = "/home/aurio/research/NUMEN/NUMEN_output/";

   // Set the in/out files
   TString inputFile = inputDir + "merg_005.root";
   TString SiCFile = inputDir + "sic_005.root";
   TString outputFile = outDir + "reconstructed_005.root";

   // Set the mapping for the TPC
   TString mapFile = "e12014_pad_mapping.xml"; //"Lookup20150611.xml";
   TString parameterFile = "ATTPC.e12014.par";
   TString planeMapFile = inputDir + "channel2pad_2.txt";
   TString parPadFileName = "/home/aurio/research/NUMEN/NUMEN_data/testParFile.txt";

   // Set directories
   TString dir = gSystem->Getenv("VMCWORKDIR");
   TString mapDir = dir + "/scripts/" + mapFile;
   TString geomDir = dir + "/geometry/";
   gSystem->Setenv("GEOMPATH", geomDir.Data());
   TString digiParFile = dir + "/parameters/" + parameterFile;
   TString geoManFile = dir + "/geometry/ATTPC_H1bar.root";

   // Create a run
   FairRunAna *run = new FairRunAna();
   run->SetSink(new FairRootFileSink(outputFile));
   run->SetGeomFile(geoManFile);

   // Set the parameter file
   FairRuntimeDb *rtdb = run->GetRuntimeDb();
   FairParAsciiFileIo *parIo1 = new FairParAsciiFileIo();

   std::cout << "Setting par file: " << digiParFile << std::endl;
   parIo1->open(digiParFile.Data(), "in");
   rtdb->setFirstInput(parIo1);
   std::cout << "Getting containers..." << std::endl;
   // We must get the container before initializing a run
   rtdb->getContainer("AtDigiPar");

   // Create the detector map
   auto fAtMapPtr = std::make_shared<AtTpcMap>();
   fAtMapPtr->ParseXMLMap(mapDir.Data());
   fAtMapPtr->GeneratePadPlane();

   auto *parser = new AtMAGNEXParserAndClusterTask(inputFile, SiCFile, planeMapFile, parPadFileName, "AtHitClusterEventH");
   parser->SetPersistence(kTRUE);
   parser->SetVerbose(kFALSE);
   //parser->SetWindowSize(2000000);             // 2000000 ps is the default value.
   //parser->SetDriftVelocity(9.25);             // 9.25 cm/us is the default value (completely arbitrary choice).
   parser->SetStripCluster(3);                 // 1 is the default value. Controls the maximum separation between hits in columns (or strip number) to be considered to be in the same cluster.
   //parser->SetTimeCluster(100000);             // 100000 ps is the default value. Controls the maximum separation between hits in TS to be considered as part of the same cluster.

   /*auto sac = std::make_unique<SampleConsensus::AtSampleConsensus>(
      SampleConsensus::Estimators::kRANSAC, AtPatterns::PatternType::kLine, RandomSample::SampleMethod::kUniform);
   auto sacTask = new AtSampleConsensusTask(std::move(sac));
   sacTask->SetPersistence(true);*/

   AtMAGNEXRansacTask *ransacTask = new AtMAGNEXRansacTask();
   ransacTask->SetPersistence(kTRUE);
   ransacTask->SetVerbose(kFALSE);
   ransacTask->SetDistanceThreshold(5.0);
   ransacTask->SetMinHitsLine(5);
   // in AtRansacTask parttern type set to line : auto patternType = AtPatterns::PatternType::kLine;
   ransacTask->SetAlgorithm(1); // 1=Homemade Ransac (default); 2=Homemade Mlesac; 3=Homemade Lmeds;
   ransacTask->SetRanSamMode(
      0); // SampleMethod { kUniform = 0, kChargeWeighted = 1, kGaussian = 2, kWeightedGaussian = 3, kWeightedY = 4 };
   ransacTask->SetChargeThreshold(0);

   run->AddTask(parser);
   // run->AddTask(sacTask);
   run->AddTask(ransacTask);

   std::cout << "***** Starting Init ******" << std::endl;
   run->Init();
   std::cout << "***** Ending Init ******" << std::endl;

   std::cout << "starting run" << std::endl;
   run->Run(0, 4000000);

   std::cout << std::endl << std::endl;
   std::cout << "Done unpacking events" << std::endl << std::endl;
   std::cout << "- Output file : " << outputFile << std::endl << std::endl;
   // -----   Finish   -------------------------------------------------------
   timer.Stop();
   Double_t rtime = timer.RealTime();
   Double_t ctime = timer.CpuTime();
   cout << endl << endl;
   cout << "Real time " << rtime << " s, CPU time " << ctime << " s" << endl;
   cout << endl;
   // ------------------------------------------------------------------------

   return 0;
}
