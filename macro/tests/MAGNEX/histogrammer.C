void histogrammer()
{
   using XYZPoint = ROOT::Math::XYZPoint;
   using XYZVector = ROOT::Math::XYZVector;

   // Define the histograms to be filled.
   TH1F *xHist = new TH1F("xHist", "xHist", 60, 0, 300);
   TH1F *yHist = new TH1F("yHist", "yHist", 100, 0, 285);
   TH1F *zHist = new TH1F("zHist", "zHist", 21, -10, 10);

   TH1F *thetaHist = new TH1F("thetaHist", "thetaHist", 180, -90, 90);
   TH1F *phiHist = new TH1F("phiHist", "phiHist", 180, -90, 90);

   TH1F *multiplicityHist = new TH1F("multiplicityHist", "multiplicityHist", 10, -0.5, 9.5);

   TH1F *xSourceHist = new TH1F("xSourceHist", "xSourceHist", 1100, -100, 1000);
   TH1F *zSourceHist = new TH1F("zSourceHist", "zSourceHist", 100, -200, 200);
   TH2F *sourceHist = new TH2F("sourceHist", "sourceHist", 250, 0, 500, 250, -300, 200);

   // Initialize the FairRun, and load the output of the recontruct.C macro.
   FairRunAna *run = new FairRunAna();

   TFile *file = new TFile("/home/aurio/research/NUMEN/NUMEN_output/reconstructed_005.root", "READ");
   TTree *tree = (TTree *)file->Get("cbmsim");
   Int_t nEvents = tree->GetEntries();
   std::cout << " Number of events : " << nEvents << std::endl;

   TTreeReader Reader("cbmsim", file);
   TTreeReaderValue<TClonesArray> eventHArray(Reader, "AtHitClusterEventH");
   TTreeReaderValue<TClonesArray> patternArray(Reader, "AtPatternEvent");

   Double_t x0Prev{}, z0Prev{};
   Double_t xPrimePrev{}, zPrimePrev{};
   Bool_t firstTrack{true};

   // Iterate over events.
   std::cout << "Event:" << std::endl;
   for (Int_t i = 0; i < nEvents; i++) {
      if (i % 10000 == 0)
         std::cout << "\x1b[1A" << "\x1b[2K" << "Event: " << i << std::endl;

      Reader.Next();
      AtEvent *event = (AtEvent *)eventHArray->At(0);
      AtPatternEvent *patternEvent = (AtPatternEvent *)patternArray->At(0);

      if (event && patternEvent) {

         // Access the tracks in each event.
         auto &tracks = patternEvent->GetTrackCand();

         for (auto &track : tracks) {

            // Access the fitted pattern line in each track.
            const AtPatterns::AtPatternLine *patternLine = dynamic_cast<const AtPatterns::AtPatternLine *>(track.GetPattern());

            // Position in the entrance plane of the TPC. The z component should be 0 by construction.
            XYZPoint entrancePoint = patternLine->GetPointAt(0);
            Double_t x = entrancePoint.X();
            Double_t y = entrancePoint.Y();
            Double_t z = entrancePoint.Z();

            // Direction vector of the track. The z component should be 1 or -1 by construction, and the x and
            // y components the tangents of the entrance angles.
            XYZVector directionVector = patternLine->GetDirection();
            Double_t theta = TMath::ATan2(directionVector.X(), directionVector.Z()) * 180. / TMath::Pi();
            Double_t phi = TMath::ATan2(directionVector.Y(), directionVector.Z()) * 180. / TMath::Pi();

            // Fill the histograms.
            //if (y < 35 || 65 < y)
               //continue;

            xHist->Fill(x);
            yHist->Fill(y);
            zHist->Fill(z);

            thetaHist->Fill(theta);
            phiHist->Fill(phi);

            // Compute another (x,z) point for the source.
            if (!firstTrack) {
               Double_t zSource = (entrancePoint.X() - x0Prev + xPrimePrev / zPrimePrev * z0Prev - directionVector.X() / directionVector.Z() * entrancePoint.Z()) / (xPrimePrev / zPrimePrev - directionVector.X() / directionVector.Z());
               Double_t xSource = x0Prev + xPrimePrev * (zSource - z0Prev) / zPrimePrev;

               xSourceHist->Fill(xSource);
               zSourceHist->Fill(zSource);
               sourceHist->Fill(xSource, zSource);
            }

            firstTrack = false;

            x0Prev = entrancePoint.X();
            z0Prev = entrancePoint.Z();
            xPrimePrev = directionVector.X();
            zPrimePrev = directionVector.Z();
         }
         multiplicityHist->Fill(tracks.size());
      }
   }

   // Close the file.
   file->Close();

   // Draw the histograms.
   TCanvas *cX = new TCanvas();
   xHist->Draw();

   TCanvas *cY = new TCanvas();
   yHist->Draw();

   TCanvas *cZ = new TCanvas();
   zHist->Draw();

   TCanvas *cTheta = new TCanvas();
   thetaHist->Draw();

   TCanvas *cPhi = new TCanvas();
   phiHist->Draw();

   TCanvas *cMultiplicity = new TCanvas();
   multiplicityHist->Draw();

   TCanvas *cXSource = new TCanvas();
   xSourceHist->Draw();

   TCanvas *cZSource = new TCanvas();
   zSourceHist->Draw();

   TCanvas *cSource = new TCanvas();
   sourceHist->Draw("zcol");
}
