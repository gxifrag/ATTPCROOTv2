#include <fstream>
#include <iostream>

double Ebin_max = 9.;
double Ebin_min = -2.;
int NumberBins  = 100;

Double_t omega(Double_t x, Double_t y, Double_t z)
{
   return sqrt(x*x + y*y + z*z - 2*x*y - 2*y*z - 2*x*z);
}

std::tuple<double, double>
kine_2b(Double_t m1, Double_t m2, Double_t m3, Double_t m4,
        Double_t K_proj, Double_t thetalab, Double_t K_eject)
{
   double Et1 = K_proj + m1;
   double Et2 = m2;
   double Et3 = K_eject + m3;
   double Et4 = Et1 + Et2 - Et3;
   double m4_ex, Ex, theta_cm;
   double s, t, u;

   s = pow(m1,2) + pow(m2,2) + 2*m2*Et1;
   u = pow(m2,2) + pow(m3,2) - 2*m2*Et3;

   m4_ex = sqrt((cos(thetalab) * omega(s, pow(m1,2), pow(m2,2)) * omega(u, pow(m2,2), pow(m3,2)) -
                 (s - pow(m1,2) - pow(m2,2)) * (pow(m2,2) + pow(m3,2) - u)) /
                  (2*pow(m2,2)) + s + u - pow(m2,2));
   Ex = m4_ex - m4;

   t = pow(m2,2) + pow(m4_ex,2) - 2*m2*Et4;

   theta_cm = TMath::Pi() - acos(
      (pow(s,2) + s*(2*t - pow(m1,2) - pow(m2,2) - pow(m3,2) - pow(m4_ex,2)) +
       (pow(m1,2) - pow(m2,2)) * (pow(m3,2) - pow(m4_ex,2))) /
      (omega(s, pow(m1,2), pow(m2,2)) * omega(s, pow(m3,2), pow(m4_ex,2))));

   theta_cm *= TMath::RadToDeg();
   return std::make_tuple(Ex, theta_cm);
}

//-------------------------------main function---------------------------------------
void calculate_efficiency_C16_pd_cleanComparison()
{
   bool guardar_en_pdf = false;
   gROOT->ProcessLine(".X /home/georgina/fair_install/ATTPCROOTv2/macro/a1975/myStyle.C");
   gStyle->SetTitleAlign(23);
   gStyle->SetTitleX(0.5);
   gROOT->SetBatch(guardar_en_pdf ? kTRUE : kFALSE);

   // ---------------------------------------------------------------
   // Histogramas: con tilt (suffijo _tilt) y sin tilt (suffijo _notilt)
   // ---------------------------------------------------------------
   auto *hexCorr2_tilt   = new TH1F("hexCorr2_tilt",   "C16(p,d) con tilt",   NumberBins, Ebin_min, Ebin_max);
   auto *hexCorr2_notilt = new TH1F("hexCorr2_notilt", "C16(p,d) sin tilt",   NumberBins, Ebin_min, Ebin_max);

   auto *hDat_thetaCM_tilt   = new TH1F("hDat_thetaCM_tilt",   "Data #theta_{CM} (tilt)",   180, 0, 180);
   auto *hDat_thetaCM_notilt = new TH1F("hDat_thetaCM_notilt", "Data #theta_{CM} (no tilt)", 180, 0, 180);

   // Simulación: un solo hSim_thetaCM (no cambia, no tiene tilt)
   auto *hSim_thetaCM = new TH1F("hSim_thetaCM", "Sim #theta_{CM}", 180, 0, 180);
   auto *h_simEx      = new TH1F("h_simEx",       "simEx",          NumberBins, Ebin_min, Ebin_max);

   // --- Masas (MeV/c²) ---
   Double_t m_p   = 938.272076;
   Double_t m_d   = 1875.612931;
   Double_t m_C15 = 13979.218707;
   Double_t m_C16 = 14914.533798;

   cout << "Q-value = " << m_C16 + m_p - m_d - m_C15 << " MeV" << endl;

   Double_t Ebeam_buff = 11.5 * 16;
   Double_t m_ej = m_d;
   int Z_ej = 1;

   // --- Energy loss ---
   double densityH2 = 3.3084e-5;
   AtTools::AtELossCATIMA elossH2(densityH2);
   elossH2.SetMaterial(catima::Material(1, 1));
   elossH2.SetProjectile(16, 6, 16.0147);

   // ---------------------------------------------------------------
   // Loop sobre datos reconstruidos
   // ---------------------------------------------------------------
   std::vector<TString> filenames;
   std::set<int> excluded = {111, 121, 148, 149};
   TChain *chain = new TChain("parquettree");

   for (int i = 104; i <= 186; i++) {
      if (excluded.count(i)) continue;
      char name[64];
      std::snprintf(name, sizeof(name), "run_%04d_2H.root", i);
      filenames.push_back(name);
      chain->Add(("/home/georgina/engine_ExUniform_pd/engine_ExUniform_pd/InterpSolver/InterpSolverRoot/" + std::string(name)).c_str());
   }
   cout << "Total entries across all the reconstruction files: " << chain->GetEntries() << endl;

   for (auto filename : filenames) {
      TFile *runFile = new TFile(
         "/home/georgina/engine_ExUniform_pd/engine_ExUniform_pd/InterpSolver/InterpSolverRoot/" + filename, "R");
      TTree *Tphysics = (TTree*)runFile->Get("parquettree");

      Double_t theta{}, Brho{}, redchi{}, zPos{}, ke{};
      Double_t vx_pos{}, vy_pos{};
      Tphysics->SetBranchAddress("polar",     &theta);
      Tphysics->SetBranchAddress("brho",      &Brho);
      Tphysics->SetBranchAddress("redchisq",  &redchi);
      Tphysics->SetBranchAddress("vertex_z",  &zPos);
      Tphysics->SetBranchAddress("ke",        &ke);
      Tphysics->SetBranchAddress("vertex_x",  &vx_pos);
      Tphysics->SetBranchAddress("vertex_y",  &vy_pos);

   for (int i = 0; i < Tphysics->GetEntries(); i++) {
      Tphysics->GetEntry(i);

      Double_t p_ej = Brho * Z_ej * 2.99792458 / 10 * 1000;
      Double_t E_ej = TMath::Sqrt(p_ej*p_ej + m_ej*m_ej) - m_ej;

      double dist3D = TMath::Sqrt(vx_pos*vx_pos + vy_pos*vy_pos + zPos*zPos) * 100.0;
      Double_t Ebeam_at_z = elossH2.GetEnergy(Ebeam_buff, dist3D);

      // --- CON tilt ---
      double kethe = 13.;
      double theta_tilt = theta - (2.0*TMath::Pi()/4000) * (E_ej - kethe);
      auto [ex_tilt, theta_cm_tilt] =
         kine_2b(m_C16, m_p, m_d, m_C15, Ebeam_at_z, theta_tilt, E_ej);

      // --- SIN tilt (theta crudo, mismo Ebeam_at_z) ---
      auto [ex_notilt, theta_cm_notilt] =
         kine_2b(m_C16, m_p, m_d, m_C15, Ebeam_at_z, theta, E_ej);

       if (zPos*100 > 2.0 && zPos*100 < 60.0 && E_ej < 14.0) {
     
         hexCorr2_tilt->Fill(ex_tilt);
         hDat_thetaCM_tilt->Fill(theta_cm_tilt);

         hexCorr2_notilt->Fill(ex_notilt);
         hDat_thetaCM_notilt->Fill(theta_cm_notilt);
      
      }
   } // events
   runFile->Close();
   } // files
   // ---------------------------------------------------------------
   // Loop sobre simulación sin detector (sin tilt, igual que antes)
   // ---------------------------------------------------------------
   TFile *fSim = new TFile(
      "/home/georgina/fair_install/ATTPCROOTv2/macro/a1975/23April_macros/efficiencies/rawSim/output_16Cpd.root",
      "READ");
   TTree *tSim = (TTree*)fSim->Get("kinematics");

   Long64_t sim_Z, sim_A;
   Double_t sim_px, sim_py, sim_pz, sim_energy;
   Double_t sim_vx{}, sim_vy{}, sim_vz{};

   tSim->SetBranchAddress("Z",        &sim_Z);
   tSim->SetBranchAddress("A",        &sim_A);
   tSim->SetBranchAddress("px",       &sim_px);
   tSim->SetBranchAddress("py",       &sim_py);
   tSim->SetBranchAddress("pz",       &sim_pz);
   tSim->SetBranchAddress("energy",   &sim_energy);
   tSim->SetBranchAddress("vertex_x", &sim_vx);
   tSim->SetBranchAddress("vertex_y", &sim_vy);
   tSim->SetBranchAddress("vertex_z", &sim_vz);

   cout << "Total simulation entries: " << tSim->GetEntries() << endl;

   for (Long64_t i = 0; i < tSim->GetEntries(); i++) {
      tSim->GetEntry(i);
      if (!(sim_Z == 1 && sim_A == 2)) continue;

      Double_t sim_p  = TMath::Sqrt(sim_px*sim_px + sim_py*sim_py + sim_pz*sim_pz);
      Double_t sim_ke = sim_energy - m_d;

      double sim_dist3D = TMath::Sqrt(sim_vx*sim_vx + sim_vy*sim_vy + sim_vz*sim_vz) * 100.0;
      Double_t sim_Ebeam_at_z = elossH2.GetEnergy(Ebeam_buff, sim_dist3D);

      Double_t sim_theta_lab = TMath::ACos(sim_pz / sim_p);
      auto [sim_ex, sim_theta_cm] =
         kine_2b(m_C16, m_p, m_d, m_C15, sim_Ebeam_at_z, sim_theta_lab, sim_ke);

      hSim_thetaCM->Fill(sim_theta_cm);
      h_simEx->Fill(sim_ex);
         
   }

   // ---------------------------------------------------------------
   // Canvas 1: Espectros de Ex comparados
   // ---------------------------------------------------------------
   TCanvas *c_ExEner = new TCanvas("c_ExEner", "Excitation Energy: tilt vs no tilt", 1200, 800);
   c_ExEner->cd();

   h_simEx->SetLineColor(kGray+2);
   h_simEx->SetLineStyle(2);
   h_simEx->GetXaxis()->SetTitle("Excitation Energy (MeV)");
   h_simEx->GetYaxis()->SetTitle("Counts");
   h_simEx->GetYaxis()->SetMaxDigits(3);
   gPad->SetLogy();
   h_simEx->Draw("HIST"); 

   hexCorr2_tilt->GetXaxis()->SetTitle("Excitation Energy (MeV)");
   hexCorr2_tilt->GetYaxis()->SetTitle("Counts");
   hexCorr2_tilt->SetLineColor(kBlue);
   hexCorr2_tilt->Draw("HIST same");

   hexCorr2_notilt->SetLineColor(kRed);
   hexCorr2_notilt->Draw("HIST same");

   auto leg_ex = new TLegend(0.65, 0.45, 0.9, 0.7);
   leg_ex->AddEntry(hexCorr2_tilt,   "Data con tilt",  "l");
   leg_ex->AddEntry(hexCorr2_notilt, "Data sin tilt",  "l");
   leg_ex->AddEntry(h_simEx,         "Sim (raw)",      "l");
   leg_ex->Draw();
   c_ExEner->Update();

   // ---------------------------------------------------------------
   // Canvas 2: theta_CM datos vs simulación (tilt vs no tilt)
   // ---------------------------------------------------------------
   TCanvas *c_thetaCM = new TCanvas("c_thetaCM", "#theta_{CM}: tilt vs no tilt", 1200, 800);
   c_thetaCM->cd();
   gPad->SetTopMargin(0.15);
   hSim_thetaCM->SetLineColor(kBlack);
   hSim_thetaCM->SetLineStyle(2);
   hSim_thetaCM->GetXaxis()->SetTitle("#theta_{CM} (#circ)");
   hSim_thetaCM->GetYaxis()->SetTitle("Counts");
   hSim_thetaCM->GetYaxis()->SetMaxDigits(3);
   hSim_thetaCM->Draw("HIST");
   
   hDat_thetaCM_tilt->SetTitle("Con tilt");
   hDat_thetaCM_tilt->SetLineColor(kBlue);
   hDat_thetaCM_tilt->Draw("HIST same");
   hDat_thetaCM_notilt->SetLineColor(kRed);
   hDat_thetaCM_notilt->Draw("HIST same");

   auto leg1 = new TLegend(0.65, 0.75, 0.9, 0.9);
   leg1->AddEntry(hSim_thetaCM,        "Sim", "l");
   leg1->AddEntry(hDat_thetaCM_tilt,   "Data with tilt", "l");
   leg1->AddEntry(hDat_thetaCM_notilt, "Data without tilt", "l");
   leg1->Draw();
   c_thetaCM->Update();

   // ---------------------------------------------------------------
   // Canvas 3: Eficiencia con tilt vs sin tilt superpuestas
   // ---------------------------------------------------------------
   TEfficiency *pEff_tilt   = new TEfficiency(*hDat_thetaCM_tilt,   *hSim_thetaCM);
   TEfficiency *pEff_notilt = new TEfficiency(*hDat_thetaCM_notilt, *hSim_thetaCM);

   pEff_tilt->SetTitle("Efficiency vs #theta_{CM};#theta_{CM} (#circ);#epsilon(#theta_{CM})");
   pEff_tilt->SetMarkerStyle(20);
   pEff_tilt->SetMarkerColor(kBlue);
   pEff_tilt->SetLineColor(kBlue);

   pEff_notilt->SetMarkerStyle(24);   // círculo vacío para distinguirlo
   pEff_notilt->SetMarkerColor(kRed);
   pEff_notilt->SetLineColor(kRed);

   TCanvas *c_eff = new TCanvas("c_eff", "Efficiency: tilt vs no tilt", 900, 600);
   c_eff->cd();
   
   // 1. Dibujas el primer objeto
   pEff_tilt->Draw("AP");
   
   // 2. IMPORTANTE: Forzar el renderizado para que ROOT cree el gráfico interno
   /*c_eff->Update(); 
   
   // 3. Ahora sí modificamos el rango del eje X de 0 a 60 grados
   if (pEff_tilt->GetPaintedGraph()) {
      pEff_tilt->GetPaintedGraph()->GetXaxis()->SetRangeUser(0, 60);
   }
*/
   // 4. Dibujas el segundo objeto encima
   pEff_notilt->Draw("P same");

   auto leg_eff = new TLegend(0.65, 0.75, 0.9, 0.9);
   leg_eff->AddEntry(pEff_tilt,   "Con tilt", "lp");
   leg_eff->AddEntry(pEff_notilt, "Sin tilt", "lp");
   leg_eff->Draw();
   
   // 5. Volvemos a actualizar para aplicar todos los cambios visuales
   c_eff->Update();

   // ---------------------------------------------------------------
   // Guardar ambas eficiencias en ROOT file
   // ---------------------------------------------------------------
   TFile *fEff = new TFile("efficiency.root", "RECREATE");
   pEff_tilt->Write("pEff_tilt");
   pEff_notilt->Write("pEff_notilt");
   fEff->Close();
}
