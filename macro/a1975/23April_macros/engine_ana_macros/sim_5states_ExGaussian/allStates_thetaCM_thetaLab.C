#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// ── Función auxiliar de cinemática ────────────────────────────────────────────
Double_t omega(Double_t x, Double_t y, Double_t z)
{
   return sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}

std::tuple<double, double>
kine_2b(Double_t m1, Double_t m2, Double_t m3, Double_t m4, Double_t K_proj, Double_t thetalab, Double_t K_eject)
{
   double Et1 = K_proj + m1;
   double Et2 = m2;
   double Et3 = K_eject + m3;
   double Et4 = Et1 + Et2 - Et3;
   double m4_ex, Ex, theta_cm;
   double s, t, u;
   s = pow(m1, 2) + pow(m2, 2) + 2 * m2 * Et1;
   u = pow(m2, 2) + pow(m3, 2) - 2 * m2 * Et3;
   m4_ex = sqrt((cos(thetalab) * omega(s, pow(m1, 2), pow(m2, 2)) * omega(u, pow(m2, 2), pow(m3, 2)) -
                 (s - pow(m1, 2) - pow(m2, 2)) * (pow(m2, 2) + pow(m3, 2) - u)) /
                   (2 * pow(m2, 2)) +
                s + u - pow(m2, 2));
   Ex = m4_ex - m4;
   t = pow(m2, 2) + pow(m4_ex, 2) - 2 * m2 * Et4;
   theta_cm = TMath::Pi() - acos((pow(s, 2) + s * (2 * t - pow(m1, 2) - pow(m2, 2) - pow(m3, 2) - pow(m4_ex, 2)) +
                                  (pow(m1, 2) - pow(m2, 2)) * (pow(m3, 2) - pow(m4_ex, 2))) /
                                 (omega(s, pow(m1, 2), pow(m2, 2)) * omega(s, pow(m3, 2), pow(m4_ex, 2))));
   theta_cm *= TMath::RadToDeg();
   return std::make_tuple(Ex, theta_cm);
}

// ── Macro principal ───────────────────────────────────────────────────────────
void allStates_thetaCM_thetaLab()
{
   gStyle->SetCanvasPreferGL(kTRUE);
   gStyle->SetOptStat(0);

   // ── Masas (MeV/c²) ────────────────────────────────────────────────────────
   Double_t m_p = 938.272076;
   Double_t m_d = 1875.612931;
   Double_t m_C15 = 13979.218707;
   Double_t m_C16 = 14914.533798;

   Double_t Ebeam_buff = 11.5 * 16; // MeV total en buffer gas
   Double_t m_b = m_d;
   Double_t m_B = m_C15;
   int Z_ej = 1;
   Double_t m_ej = m_d;

   double densityH2 = 3.553e-5; // g/cm³
   AtTools::AtELossCATIMA elossH2(densityH2);
   elossH2.SetMaterial(catima::Material(1, 1));
   elossH2.SetProjectile(16, 6, 16.0147);

   // ── Definición de los estados ─────────────────────────────────────────────
   struct StateInfo {
      TString label;  // nombre legible
      TString folder; // subcarpeta dentro de engine_16C_pd_sim
      int color;
   };

   std::vector<StateInfo> states = {
      {"GS (0 keV)", "engine_GS", kBlack},
      {"740 keV", "engine_state740keV", kRed},
      {"3103 keV", "engine_state3103keV", kBlue},
      {"4780 keV", "engine_state4780keV", kGreen + 2},
      {"6841 keV", "engine_state6841keV", kMagenta + 1},
   };

   // ── Cargar líneas cinemáticas ─────────────────────────────────────────────
   std::vector<std::string> kineFiles = {"C16_pd_C15_gs_Ebeam11_5.txt", "C16_pd_C15_740keV_Ebeam11_5.txt",
                                         "C16_pd_C15_3103keV_Ebeam11_5.txt", "C16_pd_C15_4780keV_Ebeam11_5.txt",
                                         "C16_pd_C15_6841keV_Ebeam11_5.txt"};
   std::vector<std::string> kineLabels = {"GS (0 keV)", "740 keV", "3103 keV", "4780 keV", "6841 keV"};
   std::vector<int> kineColors = {kBlack, kRed, kBlue, kGreen + 2, kMagenta + 1};
   std::vector<TGraph *> kineGraphs;

   for (size_t i = 0; i < kineFiles.size(); i++) {
      TString fpath =
         Form("/home/georgina/fair_install/ATTPCROOTv2/macro/Kinematics/Decay_kinematics/%s", kineFiles[i].c_str());
      std::ifstream kin(fpath.Data());
      if (kin.fail()) {
         std::cout << "Warning: no kinematic file for " << kineLabels[i] << std::endl;
         kineGraphs.push_back(nullptr);
         continue;
      }
      std::vector<double> tLabRec, eLabRec;
      double tCMS, tLR, eLR, tLS, eLS;
      while (kin >> tCMS >> tLR >> eLR >> tLS >> eLS) {
         tLabRec.push_back(tLR); // theta_lab en grados
         eLabRec.push_back(eLR); // KE del deuterón en MeV
      }
      TGraph *g = new TGraph(tLabRec.size(), tLabRec.data(), eLabRec.data());
      g->SetLineColor(kineColors[i]);
      g->SetLineWidth(2);
      g->SetTitle(kineLabels[i].c_str());
      kineGraphs.push_back(g);
   }

   const TString baseDir = "/home/georgina/my_sim/engine_16C_pd_sim/";

   // ── Histogramas acumulados (todos los estados juntos) ─────────────────────
   auto *h_thetaCM_all = new TH2F("h_thetaCM_all", "All states: #theta_{CM} vs E_{x};#theta_{CM} (deg);E_{x} (MeV)",
                                  100, 0, 180, 200, -3, 10);
   auto *h_thetaLab_all = new TH2F("h_thetaLab_all", "All states: #theta_{lab} vs E_{x};#theta_{lab} (deg);E_{x} (MeV)",
                                   100, 0, 50, 200, -3, 10);

   // ── Histogramas acumulados (todos los estados juntos con cuts) ─────────────────────
   auto *h_thetaCM_all_cuts = new TH2F(
      "h_thetaCM_all_cuts", "All states: #theta_{CM} vs E_{x};#theta_{CM} (deg);E_{x} (MeV)", 100, 0, 180, 200, -3, 10);
   auto *h_thetaLab_all_cuts =
      new TH2F("h_thetaLab_all_cuts", "All states: #theta_{lab} vs E_{x};#theta_{lab} (deg);E_{x} (MeV)", 100, 0, 50,
               200, -3, 10);

   // ── Histogramas acumulados (todos los estados juntos) ─────────────────────
   auto *h_KEthetaCM_all =
      new TH2F("h_KEthetaCM_all", "All states: #theta_{CM} vs KE;#theta_{CM} (deg);KE (MeV)", 100, 0, 180, 200, 0, 80);
   auto *h_KEthetaLab_all = new TH2F("h_KEthetaLab_all", "All states: #theta_{lab} vs KE;#theta_{lab} (deg);KE (MeV)",
                                     100, 0, 50, 200, 0, 80);

   // ── Histogramas acumulados (todos los estados juntos con cuts) ─────────────────────
   auto *h_KEthetaCM_all_cuts = new TH2F(
      "h_KEthetaCM_all_cuts", "All states: #theta_{CM} vs KE;#theta_{CM} (deg);KE (MeV)", 100, 0, 180, 200, 0, 25);
   auto *h_KEthetaLab_all_cuts = new TH2F(
      "h_KEthetaLab_all_cuts", "All states: #theta_{lab} vs KE;#theta_{lab} (deg);KE (MeV)", 100, 0, 50, 200, 0, 25);

   // ── Histogramas acumulados KE vs Ex (todos los estados juntos) ─────────────────────

   auto *h_KE_Ex = new TH2F("h_KE_Ex", "All states: KE vs E_{x};KE (MeV);E_{x} (MeV)", 200, 0, 80, 200, -3, 10);
   auto *h_KE_Ex_cuts =
      new TH2F("h_KE_Ex_cuts", "All states: KE vs E_{x};KE (MeV);E_{x} (MeV)", 200, 0, 80, 200, -3, 10);

   // ── Histogramas excitation energy  ─────────────────────
   auto *h_Ex = new TH1F("h_Ex_GS", "GS: E_{x};E_{x} (MeV);Counts", 300, -0.5, 8);

   // Histogramas por estado (para canvas individual)
   std::vector<TH2F *> v_thetaCM, v_thetaLab;

   for (size_t s = 0; s < states.size(); s++) {
      TString hname_cm = TString::Format("h_thetaCM_%zu", s);
      TString hname_lab = TString::Format("h_thetaLab_%zu", s);
      v_thetaCM.push_back(
         new TH2F(hname_cm, states[s].label + ";#theta_{CM} (deg);E_{x} (MeV)", 100, 0, 180, 200, -3, 10));
      v_thetaLab.push_back(
         new TH2F(hname_lab, states[s].label + ";#theta_{lab} (deg);E_{x} (MeV)", 100, 0, 50, 200, -3, 10));
   }

   // ── Loop sobre estados ────────────────────────────────────────────────────
   for (size_t s = 0; s < states.size(); s++) {

      TString stateDir = baseDir + states[s].folder + "/";
      std::cout << "\n=== Estado: " << states[s].label << "  (" << stateDir << ")" << std::endl;

      TChain *chain = new TChain("parquettree");
      int filesAdded = 0;
      for (int i = 0; i <= 19; i++) {
         TString fname = TString::Format("run_%04d_2H.root", i);
         TString fpath = stateDir + fname;
         if (chain->Add(fpath, -1) > 0)
            filesAdded++;
      }
      std::cout << "  Ficheros añadidos: " << filesAdded << "  |  Entradas totales: " << chain->GetEntries()
                << std::endl;

      if (chain->GetEntries() == 0) {
         std::cout << "  [WARNING] TChain vacío, saltando estado." << std::endl;
         delete chain;
         continue;
      }

      // Ramas
      Double_t theta{}, phi{}, Brho{}, redchi{}, zPos{}, ke{};
      Double_t vx_pos{}, vy_pos{};
      chain->SetBranchAddress("polar", &theta);
      chain->SetBranchAddress("azimuthal", &phi);
      chain->SetBranchAddress("brho", &Brho);
      chain->SetBranchAddress("redchisq", &redchi);
      chain->SetBranchAddress("vertex_z", &zPos);
      chain->SetBranchAddress("ke", &ke);
      chain->SetBranchAddress("vertex_x", &vx_pos);
      chain->SetBranchAddress("vertex_y", &vy_pos);

      Long64_t nEntries = chain->GetEntries();
      for (Long64_t i = 0; i < nEntries; i++) {
         chain->GetEntry(i);

         Double_t p_ej = Brho * Z_ej * 2.99792458 / 10 * 1000;
         Double_t E_ej = TMath::Sqrt(p_ej * p_ej + m_ej * m_ej) - m_ej;

         double dist3D = TMath::Sqrt(vx_pos * vx_pos + vy_pos * vy_pos + zPos * zPos) * 100.0; // cm
         Double_t Ebeam_at_z = elossH2.GetEnergy(Ebeam_buff, dist3D * 10.0);                   // mm

         auto [ex, theta_cm] = kine_2b(m_C16, m_p, m_b, m_B, Ebeam_at_z, theta, E_ej);

         if (TMath::IsNaN(ex) || TMath::IsNaN(theta_cm) || TMath::IsNaN(E_ej) || TMath::IsNaN(Ebeam_at_z))
            continue;

         double theta_lab_deg = theta * TMath::RadToDeg();

         // Corte de calidad básico (igual que en la macro original)

         v_thetaCM[s]->Fill(theta_cm, ex);
         v_thetaLab[s]->Fill(theta_lab_deg, ex);

         h_thetaCM_all->Fill(theta_cm, ex);
         h_thetaLab_all->Fill(theta_lab_deg, ex);
         h_KEthetaCM_all->Fill(theta_cm, E_ej);
         h_KEthetaLab_all->Fill(theta_lab_deg, E_ej);
         h_KE_Ex->Fill(E_ej, ex);

         if (zPos * 100 > 2.0 && zPos * 100 < 60.0 && E_ej > 5.0 && E_ej < 20.0) {
            h_thetaCM_all_cuts->Fill(theta_cm, ex);
            h_thetaLab_all_cuts->Fill(theta_lab_deg, ex);
            h_KEthetaCM_all_cuts->Fill(theta_cm, E_ej);
            h_KEthetaLab_all_cuts->Fill(theta_lab_deg, E_ej);
            h_KE_Ex_cuts->Fill(E_ej, ex);

            h_Ex->Fill(ex);
         }
      }

      std::cout << "  Listo." << std::endl;
      delete chain;
   }

   // ── Canvas: un pad por estado, thetaCM ───────────────────────────────────
   TCanvas *c_cm = new TCanvas("c_thetaCM_states", "#theta_{CM} vs E_{x} por estado", 1200, 800);
   c_cm->Divide(3, 2); // 5 estados + 1 para el acumulado
   for (size_t s = 0; s < states.size(); s++) {
      c_cm->cd(s + 1);
      gPad->SetLogz(1);
      v_thetaCM[s]->SetTitle(states[s].label + " — #theta_{CM} vs E_{x}");
      v_thetaCM[s]->Draw("COLZ");
   }
   c_cm->cd(6);
   gPad->SetLogz(1);
   h_thetaCM_all->SetTitle("All states — #theta_{CM} vs E_{x}");
   h_thetaCM_all->Draw("COLZ");
   c_cm->Update();
   c_cm->SaveAs("ExThetaCM_allStates_simulation_Exgaussian.png");

   // ── Canvas: un pad por estado, thetaLab ──────────────────────────────────
   TCanvas *c_lab = new TCanvas("c_thetaLab_states", "#theta_{lab} vs E_{x} por estado", 1200, 800);
   c_lab->Divide(3, 2);
   for (size_t s = 0; s < states.size(); s++) {
      c_lab->cd(s + 1);
      gPad->SetLogz(1);
      v_thetaLab[s]->SetTitle(states[s].label + " — #theta_{lab} vs E_{x}");
      v_thetaLab[s]->Draw("COLZ");
   }
   c_lab->cd(6);
   gPad->SetLogz(1);
   h_thetaLab_all->SetTitle("All states — #theta_{lab} vs E_{x}");
   h_thetaLab_all->Draw("COLZ");
   c_lab->Update();
   c_lab->SaveAs("ExThetaLab_allStates_simulation_Exgaussian.png");

   // ── Canvas acumulado lado a lado ──────────────────────────────────────────
   TCanvas *c_all = new TCanvas("c_all", "All states combined", 1200, 800);
   c_all->Divide(2, 2);
   c_all->cd(1);
   gPad->SetLogz(1);
   h_thetaCM_all->Draw("COLZ");
   c_all->cd(2);
   gPad->SetLogz(1);
   h_thetaLab_all->Draw("COLZ");

   c_all->cd(3);
   gPad->SetLogz(1);
   h_thetaCM_all_cuts->Draw("COLZ");
   c_all->cd(4);
   gPad->SetLogz(1);
   h_thetaLab_all_cuts->Draw("COLZ");
   c_all->Update();
   c_all->SaveAs("ExThetaCM_allStates_combined_simulation_Exgaussian.png");

   TCanvas *c_all_KE = new TCanvas("c_all_KE", "All states combined - KE", 1200, 800);
   c_all_KE->Divide(2, 2);
   c_all_KE->cd(1);
   gPad->SetLogz(1);
   h_KEthetaCM_all->Draw("COLZ");

   c_all_KE->cd(2);
   gPad->SetLogz(1);
   h_KEthetaLab_all->Draw("COLZ");
   for (size_t i = 0; i < kineGraphs.size(); i++) {
      if (kineGraphs[i])
         kineGraphs[i]->Draw("L SAME");
   }
   auto *legKine = new TLegend(0.45, 0.65, 0.90, 0.90);
   for (size_t i = 0; i < kineGraphs.size(); i++) {
      if (kineGraphs[i])
         legKine->AddEntry(kineGraphs[i], kineLabels[i].c_str(), "l");
   }
   legKine->Draw();

   c_all_KE->cd(3);
   gPad->SetLogz(1);
   h_KEthetaCM_all_cuts->Draw("COLZ");

   c_all_KE->cd(4);
   gPad->SetLogz(1);
   h_KEthetaLab_all_cuts->Draw("COLZ");

   for (size_t i = 0; i < kineGraphs.size(); i++) {
      if (kineGraphs[i])
         kineGraphs[i]->Draw("L SAME");
   }

   c_all_KE->Update();
   c_all_KE->SaveAs("KEThetaCM_allStates_combined_simulation_Exgaussian.png");

   TCanvas *c_KE_Ex = new TCanvas("c_KE_Ex", "All states: KE vs E_{x}", 1200, 800);
   c_KE_Ex->cd();
   gPad->SetLogz(1);
   h_KE_Ex->Draw("COLZ");
   c_KE_Ex->Update();
   c_KE_Ex->SaveAs("KE_vs_Ex_allStates_simulation_Exgaussian.png");

   std::vector<TGraph *> dummies;
   for (const auto &state : states) {
      TGraph *g = new TGraph(1);
      g->SetLineColor(state.color);
      g->SetMarkerColor(state.color);
      g->SetMarkerStyle(20);
      g->SetMarkerSize(0.8);
      dummies.push_back(g);
   }

   TLegend *legX2 = new TLegend(0.60, 0.60, 0.95, 0.90);
   legX2->SetBorderSize(0);
   legX2->SetFillStyle(0);
   legX2->SetTextSize(0.04);
   legX2->AddEntry(dummies[4], "6841 keV", "p");
   legX2->AddEntry(dummies[3], "4780 keV", "p");
   legX2->AddEntry(dummies[2], "3103 keV", "p");
   legX2->AddEntry(dummies[1], "740 keV", "p");
   legX2->AddEntry(dummies[0], "GS (0 keV)", "p");
   legX2->Draw("same");
   c_KE_Ex->Update();
   c_KE_Ex->SaveAs("KE_vs_Ex_allStates_simulation_Exgaussian.png");

   // Fit GS
   TF1 *gaus_GS = new TF1("gaus_GS", "gaus", -0.3, 0.3);
   gaus_GS->SetLineColor(kRed);
   h_Ex->Fit("gaus_GS", "R");
   Double_t mean_GS = gaus_GS->GetParameter(1);
   Double_t sigma_GS = gaus_GS->GetParameter(2);
   cout << "GS:  mean = " << mean_GS << " MeV,  sigma = " << sigma_GS << " MeV,  FWHM = " << 2.355 * sigma_GS << " MeV"
        << endl;

   // Fit 2o estado (~0.74 MeV)
   TF1 *gaus_S2 = new TF1("gaus_S2", "gaus", 0.5, 1.1);
   gaus_S2->SetLineColor(kBlue);
   h_Ex->Fit("gaus_S2", "R+");
   Double_t mean_S2 = gaus_S2->GetParameter(1);
   Double_t sigma_S2 = gaus_S2->GetParameter(2);
   cout << "2nd: mean = " << mean_S2 << " MeV,  sigma = " << sigma_S2 << " MeV,  FWHM = " << 2.355 * sigma_S2 << " MeV"
        << endl;

   // Guardar resultados
   Double_t centroids[2] = {mean_GS, mean_S2};
   Double_t sigmas[2] = {sigma_GS, sigma_S2};
   Double_t err_centroids[2] = {gaus_GS->GetParError(1), gaus_S2->GetParError(1)};
   Double_t err_sigmas[2] = {gaus_GS->GetParError(2), gaus_S2->GetParError(2)};

   TGraphErrors *gr_sigma = new TGraphErrors(2, centroids, sigmas, err_centroids, err_sigmas);
   gr_sigma->SetMarkerStyle(20);
   gr_sigma->SetMarkerColor(kBlack);
   gr_sigma->SetLineColor(kBlack);
   gr_sigma->GetXaxis()->SetTitle("E_{x} (MeV)");
   gr_sigma->GetYaxis()->SetTitle("#sigma (MeV)");
   gr_sigma->SetTitle("Resolution vs E_{x}");

   Double_t centroids_exp[2] = {0.0362266, 0.816186};
   Double_t sigmas_exp[2] = {0.21964, 0.168605};
   Double_t err_centroids_exp[2] = {0.00782033, 0.00557057};
   Double_t err_sigmas_exp[2] = {0.00713929, 0.00407833};

   TGraphErrors *gr_sigma_exp = new TGraphErrors(2, centroids_exp, sigmas_exp, err_centroids_exp, err_sigmas_exp);
   gr_sigma_exp->SetMarkerStyle(21);
   gr_sigma_exp->SetMarkerColor(kRed);
   gr_sigma_exp->SetLineColor(kRed);

   // Canvas combinado
   TCanvas *c_combined = new TCanvas("c_combined", "Excitation Energy & Resolution", 1400, 600);
   c_combined->Divide(2, 1);

   c_combined->cd(1);
   h_Ex->Draw();
   gaus_GS->Draw("same");
   gaus_S2->Draw("same");

   c_combined->cd(2);
   gr_sigma->Draw("AP");
   gr_sigma->GetYaxis()->SetRangeUser(0.0, 0.30);
   gr_sigma->GetXaxis()->SetRangeUser(-0.1, 1.1);

   TF1 *lin = new TF1("lin", "pol1", -0.1, 1.0);
   lin->SetLineColor(kBlue);
   gr_sigma->Fit("lin", "R");
   cout << "Sigma(Ex) = " << lin->GetParameter(0) << " + " << lin->GetParameter(1) << " * Ex" << endl;

   TLatex *tex = new TLatex();
   tex->SetNDC();
   tex->SetTextSize(0.04);
   tex->SetTextColor(kBlue);
   tex->DrawLatex(0.18, 0.80, Form("#sigma = %.3f + %.3f E_{x}", lin->GetParameter(0), lin->GetParameter(1)));

   gr_sigma_exp->Draw("P same");

   TLegend *leg = new TLegend(0.55, 0.7, 0.88, 0.88);
   leg->AddEntry(gr_sigma, "Simulation", "p");
   leg->AddEntry(gr_sigma_exp, "Experiment", "p");
   leg->Draw();

   c_combined->Update();
   c_combined->SaveAs("ExGaussianFitting_and_ComparisonExperimentalSigmas.png");
}
