#include <fstream>
#include <iostream>

double Ebin_max = 9.;
double Ebin_min = -1.;
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
void calculate_efficiency_C16_pd_clean()
{
   bool guardar_en_pdf = false;
   gROOT->ProcessLine(".X /home/georgina/fair_install/ATTPCROOTv2/macro/a1975/myStyle.C");
   gStyle->SetTitleAlign(23);
   gStyle->SetTitleX(0.5);
   gROOT->SetBatch(guardar_en_pdf ? kTRUE : kFALSE);

   // --- Histogramas que sí se usan ---
   auto *hexCorr2     = new TH1F("hexCorr2",     "C16(p,d)",       NumberBins, Ebin_min, Ebin_max);
   auto *hDat_lab_corr = new TH1F("hDat_lab_corr", "Data #theta_{CM}", 45, 0, 45);
   auto *hSim_lab_corr = new TH1F("hSim_lab_corr", "Sim #theta_{CM}",  180, 0, 180);
   auto *h_simEx      = new TH1F("h_simEx",       "simEx",           NumberBins, Ebin_min, Ebin_max);
   auto *ExCorrvsZpos = new TH2F("ExCorrvsZpos", "ExCorrvsZpos", 130, -3.0, 10.0, 110, -5, 105);

   auto * h_Eff_lab_corr = new TH1F("h_Eff_lab_corr", "Efficiency vs #theta_{CM}", 180, 0, 180);

      
   // --- Histogramas 2D para eficiencia en (Ex, theta_lab) ---
   int nBinsEx  = NumberBins; // mismo número de bins que en hexCorr2

auto *hDat_2D = new TH2F("hDat_2D", "Data; Ex (MeV); #theta_{lab} (#circ)",
                          nBinsEx, Ebin_min, Ebin_max,
                          45, 0, 45);  // ajusta según tu rango real
auto *hSim_2D = new TH2F("hSim_2D", "Sim; Ex (MeV); #theta_{lab} (#circ)",
                          nBinsEx, Ebin_min, Ebin_max,
                          45, 0, 45);

   // --- Masas (MeV/c²) ---
   Double_t m_p   = 938.272076;
   Double_t m_d   = 1875.612931;
   Double_t m_C15 = 13979.218707;
   Double_t m_C16 = 14914.533798;

   cout << "Q-value = " << m_C16 + m_p - m_d - m_C15 << " MeV" << endl;

   Double_t Ebeam_buff = 11.5 * 16; // MeV
   Double_t m_b  = m_d;
   Double_t m_B  = m_C15;
   Double_t m_ej = m_d;
   int Z_ej = 1;

   // --- Energy loss ---
   double densityH2 = 3.553e-5; // g/cm³
   AtTools::AtELossCATIMA elossH2(densityH2);
   elossH2.SetMaterial(catima::Material(1, 1));
   elossH2.SetProjectile(16, 6, 16.0147);

   // --- Datos reconstruidos ---
   std::vector<TString> filenames;
   //std::set<int> excluded = {111, 121, 148, 149};
   TChain *chain = new TChain("parquettree");
   for (int i = 104; i <= 186; i++) {
   //for (int i = 0; i <= 19; i++) {
      //if (excluded.count(i)) continue;
      char name[64];
     std::snprintf(name, sizeof(name), "run_%04d_2H.root", i);
      filenames.push_back(name);
      chain->Add(("/home/georgina/engine_ExUniform_pd/engine_ExUniform_pd/InterpSolver/InterpSolverRoot/" + std::string(name)).c_str());
      //chain->Add(("/home/georgina/my_sim/engine_Ex_GS_C16_pd/InterpSolver/interpSolverRoot/" + std::string(name)).c_str());

   }
   cout << "Total entries (reconstruction): " << chain->GetEntries() << endl;

   Double_t theta{}, phi{}, Brho{}, redchi{}, zPos{}, ke{};
   Double_t vx_pos{}, vy_pos{};
   chain->SetBranchAddress("polar",     &theta);
   chain->SetBranchAddress("azimuthal", &phi);
   chain->SetBranchAddress("brho",      &Brho);
   chain->SetBranchAddress("redchisq",  &redchi);
   chain->SetBranchAddress("vertex_z",  &zPos);
   chain->SetBranchAddress("ke",        &ke);
   chain->SetBranchAddress("vertex_x",  &vx_pos); //m
   chain->SetBranchAddress("vertex_y",  &vy_pos);

   for (Long64_t i = 0; i < chain->GetEntries(); i++) {
      chain->GetEntry(i);

      Double_t p_ej = Brho * Z_ej * 2.99792458 / 10 * 1000;
      Double_t E_ej = TMath::Sqrt(p_ej*p_ej + m_ej*m_ej) - m_ej;

      double dist3D = TMath::Sqrt(vx_pos*vx_pos + vy_pos*vy_pos + zPos*zPos) * 100.0;
      Double_t Ebeam_at_z = elossH2.GetEnergy(Ebeam_buff, dist3D*10.0); // Convertir dist3D a mm para la corrección de energía

      double kethe = 13.;
      double theta_lab_corr_tilt = theta; //theta-(2.0*TMath::Pi()/4000) * (E_ej - kethe);

      auto [ex_corr_tilt, theta_cm_corr_tilt] =
         kine_2b(m_C16, m_p, m_b, m_B, Ebeam_at_z, theta_lab_corr_tilt, E_ej);

      if (zPos*100 > 2.0 && zPos*100 < 60.0 && E_ej < 14.0) {
         hexCorr2->Fill(ex_corr_tilt);
         hDat_lab_corr->Fill(theta_lab_corr_tilt);
         ExCorrvsZpos->Fill(ex_corr_tilt, zPos*100);

         double theta_lab_deg = theta_lab_corr_tilt * TMath::RadToDeg();
         hDat_2D->Fill(ex_corr_tilt, theta_lab_deg);  
      }
   }
   // --- Simulación sin detector ---
   TFile *fSim = new TFile(
      "/home/georgina/fair_install/ATTPCROOTv2/macro/a1975/23April_macros/efficiencies/rawSim/output_16Cpd.root",
      //"/home/georgina/my_sim/output_generateKin_gs.root",
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

      double sim_dist3D = TMath::Sqrt(sim_vx*sim_vx + sim_vy*sim_vy + sim_vz*sim_vz) * 100.0; //cm
      Double_t sim_Ebeam_at_z = elossH2.GetEnergy(Ebeam_buff, sim_dist3D*10.0); // Convertir sim_dist3D a mm para la corrección de energía

      Double_t sim_theta_lab = TMath::ACos(sim_pz / sim_p);
      auto [sim_ex, sim_theta_cm] =
         kine_2b(m_C16, m_p, m_d, m_C15, sim_Ebeam_at_z, sim_theta_lab, sim_ke);

      hSim_lab_corr->Fill(sim_theta_cm);
      h_simEx->Fill(sim_ex);
      ExCorrvsZpos->Fill(sim_ex, sim_vz*100.0);
      //cout << "sim_ex = " << sim_ex << "  sim_theta_cm = " << sim_theta_cm << "  sim_vz = " << sim_vz*100.0 << endl;

     double sim_theta_lab_deg = sim_theta_lab * TMath::RadToDeg();
     hSim_2D->Fill(sim_ex, sim_theta_lab_deg);
   }

   Double_t scale = (hSim_2D->Integral() > 0)
                 ? hDat_2D->Integral() / hSim_2D->Integral()
                 : 1.0;
   hSim_2D->Scale(scale);

   TH2F *hEff_2D = (TH2F*)hDat_2D->Clone("hEff_2D");
   hEff_2D->Divide(hSim_2D);

   // --- Canvas 1: Espectro de energía de excitación ---
   TCanvas *c_ExEner = new TCanvas("ExEner", "Excited Energy spectra", 1200, 800);
   c_ExEner->cd();
   h_simEx->GetXaxis()->SetTitle("Excitation Energy (MeV)");
   h_simEx->GetYaxis()->SetTitle("Counts");
   h_simEx->SetLineColor(kRed);
   h_simEx->Draw("HIST");
   hexCorr2->Draw("HIST same");

   auto leg1 = new TLegend(0.6, 0.7, 0.9, 0.9);
   leg1->AddEntry(hexCorr2, "Data (corrected)", "l");
   leg1->AddEntry(h_simEx,  "Sim (raw)", "l");
   leg1->Draw();
   c_ExEner->Update();

   // --- Canvas 2: Check theta_CM ---
   TCanvas *c_eff_check = new TCanvas("c_eff_check", "Efficiency check", 1200, 600);
   c_eff_check->Divide(2, 1);

   c_eff_check->cd(1);
   hSim_lab_corr->SetLineColor(kRed);
   hSim_lab_corr->GetXaxis()->SetTitle("#theta_{CM} (#circ)");
   hSim_lab_corr->GetYaxis()->SetTitle("Counts");
   hSim_lab_corr->GetYaxis()->SetMaxDigits(3);
   hSim_lab_corr->SetTitle("Simulacion sin detector");
   gPad->SetTopMargin(0.15);
   hSim_lab_corr->Draw("HIST");

   c_eff_check->cd(2);
   hDat_lab_corr->SetLineColor(kBlue);
   hDat_lab_corr->GetXaxis()->SetTitle("#theta_{CM} (#circ)");
   hDat_lab_corr->GetYaxis()->SetTitle("Counts");
   hDat_lab_corr->GetYaxis()->SetMaxDigits(3);
   hDat_lab_corr->SetTitle("Datos (con cortes)");
   gPad->SetTopMargin(0.15);
   hDat_lab_corr->Draw("HIST");
   c_eff_check->Update();

   // --- Canvas 3: Eficiencia ---
   TCanvas *c_eff = new TCanvas("c_eff", "Efficiency vs #theta_{CM}", 800, 600);
   c_eff->cd();
   TEfficiency *pEff = new TEfficiency(*hDat_lab_corr, *hSim_lab_corr);
   pEff->SetMarkerStyle(20);
   pEff->SetMarkerColor(kBlue);
   pEff->SetLineColor(kBlue);
   pEff->Draw("AP");
   c_eff->Update();

   // --- Canvas: Eficiencia 2D (Ex vs theta_lab) ---
   TCanvas *c_eff2D = new TCanvas("c_eff2D", "Efficiency vs Ex and #theta_{lab}", 900, 700);
   c_eff2D->cd();

   // Dividir bin a bin: eff = datos / sim (con protección contra división por cero)

   hEff_2D->Divide(hSim_2D);
   hEff_2D->SetTitle("Efficiency; Ex (MeV); #theta_{lab} (#circ); Efficiency");
   hEff_2D->SetMaximum(1.0);
   hEff_2D->SetMinimum(0.0);
   hEff_2D->Draw("COLZ");
   c_eff2D->Update();
  

   // --- Guardar eficiencia ---
   TFile *fEff = new TFile("efficiency.root", "RECREATE");
   pEff->Write("pEff");
   fEff->Close();
}