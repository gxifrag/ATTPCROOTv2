#include <stdio.h>

void tracking()
{

   FairRunAna *run = new FairRunAna();

   std::ofstream hitsFile;
   hitsFile.open("hits.txt");

   TFile *file = new TFile("/home/aurio/research/NUMEN/NUMEN_output/test.root", "READ");
   TTree *tree = (TTree *)file->Get("cbmsim");
   Int_t nEvents = tree->GetEntries();
   std::cout << " Number of events : " << nEvents << std::endl;

   TTreeReader Reader("cbmsim", file);
   TTreeReaderValue<TClonesArray> eventHArray(Reader, "AtEventH");
   TTreeReaderValue<TClonesArray> patternArray(Reader, "AtPatternEvent");

   TCanvas *c = new TCanvas("c", "c", 1600, 1000);
   c->Divide(2, 1);
   //c->SetSupportGL(true);


   TH2F *dummyHist = new TH2F("dummyHist", "", 1, 0, 110, 1, 0, 600);
   dummyHist->SetStats(kFALSE);
   dummyHist->SetMinimum(0);
   dummyHist->SetMaximum(1100);
   dummyHist->GetXaxis()->SetTitle("z");
   dummyHist->GetYaxis()->SetTitle("x");
   dummyHist->GetZaxis()->SetTitle("time");

   TH2F *chargeHist = new TH2F("chargeHist", "", 5, 0, 110, 300, 0, 600);
   chargeHist->SetStats(kFALSE);
   chargeHist->GetXaxis()->SetTitle("z");
   chargeHist->GetYaxis()->SetTitle("x");
   chargeHist->SetLineColorAlpha(0, 0);


   c->cd(2);

   c->Modified();
   c->Update();

   EColor colorList[6] = {kRed, kBlue, kGreen, kMagenta, kBlack, kPink};

   double maxX{}, maxZ{}, maxT{};
   double minX{1000}, minZ{1000}, minT{1000};
   for (Int_t i = 0; i < nEvents; i++) {

      Reader.Next();

      AtEvent *event = (AtEvent *)eventHArray->At(0);
      AtPatternEvent *patternEvent = (AtPatternEvent *)patternArray->At(0);

      chargeHist->Reset("ICESM");

      if (event && patternEvent) {

         auto &hitArray = event->GetHits();
         auto &tracks = patternEvent->GetTrackCand();
         std::cout << " Number of hits : " << hitArray.size() << std::endl;
         std::cout << " Number of tracks : " << tracks.size() << std::endl;

         /*for(auto &hit : hitArray) {
             auto pos = hit->GetPosition();
             auto charge = hit->GetCharge();
             auto time = hit->GetTimeStamp();
             hitsFile  << pos.X() << " " << pos.Y() << " " << pos.Z() << "  " << time << "  " << charge << std::endl;
         }*/

         c->Divide(2, 1);
         c->cd(1);
         dummyHist->Draw("SURF");

         c->cd(2);
         chargeHist->Draw("ZCOL");

         std::vector<TGraph2D*> scatters{};
         int j = 0;
         for (auto &track : tracks) {
            TGraph2D *scatterPlot3D = new TGraph2D();
            scatters.push_back(scatterPlot3D);
            scatterPlot3D->SetMarkerStyle(20);
            scatterPlot3D->SetMarkerColor(colorList[j]);

            TGraph *scatterPlot = new TGraph();
            scatterPlot->SetMarkerStyle(20);
            scatterPlot->SetMarkerColor(colorList[j]);

            auto &points = track.GetHitArray();
            std::cout << "Track ID: " << track.GetTrackID() << std::endl;

            int k = 0;
            for (auto &point : points) {
               auto pos = point->GetPosition();
               auto charge = point->GetCharge();
               auto time = point->GetTimeStamp();
               std::cout << pos.X() << " " << pos.Y() << " " << pos.Z() << "  " << time << "  " << charge << std::endl;

               scatterPlot3D->SetPoint(k, pos.Z(), pos.X(), time);
               scatterPlot->SetPoint(k, pos.Z(), pos.X());
               chargeHist->Fill(pos.Z(), pos.X(), charge);

               if(pos.X() < minX) minX = pos.X();
               if(pos.Z() < minZ) minZ = pos.Z();
               if(time    < minT) minT = time;
               if(pos.X() > maxX) maxX = pos.X();
               if(pos.Z() > maxZ) maxZ = pos.Z();
               if(time    > maxT) maxT = time;

               k++;
            }

            c->cd(1);
            scatterPlot3D->Draw("same p");

            c->cd(2);
            scatterPlot->Draw("same p");

            c->Modified();
            j++;
         }





         gPad->Modified();
         gPad->Update();

         gSystem->ProcessEvents();

         std::cout << "Continue?";
         gPad->WaitPrimitive();
         //std::cin.get();
         /*for(auto scatter : scatters) {
            gPad->GetListOfPrimitives()->Remove(scatter);
            gPad->Modified();
         }*/

         c->Clear();

         gPad->Modified();
         gPad->Update();
      }
   }
   std::cout << "X: "    << minX << " to " << maxX << "\n" <<
                "time: " << minT << " to " << maxT << "\n" <<
                "Z: "    << minZ << " to " << maxZ << "\n" << std::endl;

   hitsFile.close();
   file->Close();
}
