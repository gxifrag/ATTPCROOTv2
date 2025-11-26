 //----------------------------------- Plots for ground state analysis --------------------------------

   Double_t lowerBound = finalParams[1] - 2 * finalParams[2];
   Double_t upperBound = finalParams[1] + 2 * finalParams[2];

   Double_t integral = gaus1->Integral(lowerBound, upperBound);
   std::cout << "Integral del pico gaussiano (±2σ): " << integral << std::endl;
   //hGS_ExEnergy->Fit(gaus1, "RQ");


   TCanvas* cAngGS = new TCanvas("cAngGS", "Ground State Angular Distribution", 1000, 600);
   cAngGS->cd();
   gPad->SetLogy();  // ← Esto sí activa escala logarítmica en Y

   // Histograma experimental
   hGS_AngularDistr->Sumw2();
   hGS_AngularDistr->GetXaxis()->SetTitle("Theta_cm (deg)");
   hGS_AngularDistr->GetYaxis()->SetTitle("Counts");
   hGS_AngularDistr->SetLineColor(kBlue+2);
   hGS_AngularDistr->SetLineWidth(2);
   hGS_AngularDistr->Draw("E1");

   // Curva teórica
   auto *g = new TGraphErrors("/home/georgina/twofnr/21.groundState", "%lg %lg"); // Si solo hay X e Y
   g->SetLineColor(kRed);
   g->SetLineWidth(2);
   g->SetMarkerStyle(20);
   g->SetMarkerColor(kBlack);
   g->Draw("SAME"); // "P" para puntos, "L SAME" si prefieres línea

   // Leyenda opcional
   auto *legend = new TLegend(0.6, 0.7, 0.88, 0.88);
   legend->AddEntry(hGS_AngularDistr, "Datos experimentales", "l");
   legend->AddEntry(g, "twofnr", "lp");
   legend->Draw();

