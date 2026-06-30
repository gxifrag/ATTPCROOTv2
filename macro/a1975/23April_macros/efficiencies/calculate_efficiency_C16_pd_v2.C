#include <fstream>
#include <iostream>

double Ebin_max = 9.;
double Ebin_min = -1.;
int NumberBins = 100;

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

//-------------------------------main function---------------------------------------
void calculate_efficiency_C16_pd_v2()
{
   bool guardar_en_pdf = false;
   gROOT->ProcessLine(".X /home/georgina/fair_install/ATTPCROOTv2/macro/a1975/myStyle.C");
   gStyle->SetTitleAlign(23);
   gStyle->SetTitleX(0.5);
   gROOT->SetBatch(guardar_en_pdf ? kTRUE : kFALSE);

   // --- Binning común para θ_CM [grados] ---
   // FIX: ambos histogramas usan θ_CM en grados, mismo rango y bins
   const int nBinsCM = 80;
   const double tCM_min = 10.;
   const double tCM_max = 90.;

   // Histogramas de θ_CM: reconstruidos (numerador) y generados (denominador)
   // FIX: nombres claros para evitar confusión lab/CM
   auto *hReco_tCM = new TH1F("hReco_tCM", "Reco #theta_{CM};#theta_{CM} (#circ);Counts", nBinsCM, tCM_min, tCM_max);
   auto *hGen_tCM = new TH1F("hGen_tCM", "Gen #theta_{CM};#theta_{CM} (#circ);Counts", nBinsCM, tCM_min, tCM_max);

   // Histogramas auxiliares de diagnóstico
   auto *hexCorr2 = new TH1F("hexCorr2", "C16(p,d) Data Ex;Ex (MeV);Counts", NumberBins, Ebin_min, Ebin_max);
   auto *h_simEx = new TH1F("h_simEx", "Sim Ex;Ex (MeV);Counts", NumberBins, Ebin_min, Ebin_max);
   auto *ExCorrvsZpos = new TH2F("ExCorrvsZpos", "Ex vs Zpos;Ex (MeV);z (cm)", 130, -3.0, 10.0, 110, -5, 105);

   // --- Masas (MeV/c²) ---
   Double_t m_p = 938.272076;
   Double_t m_d = 1875.612931;
   Double_t m_C15 = 13979.218707;
   Double_t m_C16 = 14914.533798;

   cout << "Q-value = " << m_C16 + m_p - m_d - m_C15 << " MeV" << endl;

   Double_t Ebeam_buff = 11.5 * 16; // MeV
   Double_t m_b = m_d;
   Double_t m_B = m_C15;
   Double_t m_ej = m_d;
   int Z_ej = 1;

   // --- Energy loss ---
   double densityH2 = 3.553e-5; // g/cm³
   AtTools::AtELossCATIMA elossH2(densityH2);
   elossH2.SetMaterial(catima::Material(1, 1));
   elossH2.SetProjectile(16, 6, 16.0147);

   // =========================================================
   // --- DATOS RECONSTRUIDOS (numerador de eficiencia) ---
   // =========================================================
   TChain *chain = new TChain("parquettree");
   for (int i = 104; i <= 186; i++) {
      char name[64];
      std::snprintf(name, sizeof(name), "run_%04d_2H.root", i);
      chain->Add(
         ("/home/georgina/engine_ExUniform_pd/engine_ExUniform_pd/InterpSolver/InterpSolverRoot/" + std::string(name))
            .c_str());
   }
   cout << "Total entries (reconstruction): " << chain->GetEntries() << endl;

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

   for (Long64_t i = 0; i < chain->GetEntries(); i++) {
      chain->GetEntry(i);

      Double_t p_ej = Brho * Z_ej * 2.99792458 / 10 * 1000;
      Double_t E_ej = TMath::Sqrt(p_ej * p_ej + m_ej * m_ej) - m_ej;

      double dist3D = TMath::Sqrt(vx_pos * vx_pos + vy_pos * vy_pos + zPos * zPos) * 100.0;
      Double_t Ebeam_at_z = elossH2.GetEnergy(Ebeam_buff, dist3D * 10.0);

      double kethe = 13.;
      double theta_lab_corr_tilt = theta - (2.0 * TMath::Pi() / 4000) * (E_ej - kethe);

      // FIX: kine_2b devuelve (Ex, theta_CM en grados) — usamos theta_cm_corr_tilt
      auto [ex_corr_tilt, theta_cm_corr_tilt] = kine_2b(m_C16, m_p, m_b, m_B, Ebeam_at_z, theta_lab_corr_tilt, E_ej);

      if (zPos * 100 > 2.0 && zPos * 100 < 60.0 && E_ej > 5.0 && E_ej < 15.0) {
         hexCorr2->Fill(ex_corr_tilt);
         ExCorrvsZpos->Fill(ex_corr_tilt, zPos * 100);

         // FIX: rellenar con θ_CM en grados (no θ_lab en radianes)
         hReco_tCM->Fill(theta_cm_corr_tilt);
      }
   }

   // =========================================================
   // --- SIMULACIÓN GENERADA (denominador de eficiencia) ---
   // =========================================================
   TFile *fSim = new TFile(
      "/home/georgina/fair_install/ATTPCROOTv2/macro/a1975/23April_macros/efficiencies/rawSim/output_16Cpd.root",
      "READ");
   TTree *tSim = (TTree *)fSim->Get("kinematics");

   Long64_t sim_Z, sim_A;
   Double_t sim_px, sim_py, sim_pz, sim_energy;
   Double_t sim_vx{}, sim_vy{}, sim_vz{};

   tSim->SetBranchAddress("Z", &sim_Z);
   tSim->SetBranchAddress("A", &sim_A);
   tSim->SetBranchAddress("px", &sim_px);
   tSim->SetBranchAddress("py", &sim_py);
   tSim->SetBranchAddress("pz", &sim_pz);
   tSim->SetBranchAddress("energy", &sim_energy);
   tSim->SetBranchAddress("vertex_x", &sim_vx);
   tSim->SetBranchAddress("vertex_y", &sim_vy);
   tSim->SetBranchAddress("vertex_z", &sim_vz);

   // Contamos total de eventos generados (deuterones) para normalización
   Long64_t N_gen_total = 0;
   cout << "Total simulation entries: " << tSim->GetEntries() << endl;

   for (Long64_t i = 0; i < tSim->GetEntries(); i++) {
      tSim->GetEntry(i);
      if (!(sim_Z == 1 && sim_A == 2))
         continue;
      N_gen_total++;

      Double_t sim_p = TMath::Sqrt(sim_px * sim_px + sim_py * sim_py + sim_pz * sim_pz);
      Double_t sim_ke = sim_energy - m_d;

      double sim_dist3D = TMath::Sqrt(sim_vx * sim_vx + sim_vy * sim_vy + sim_vz * sim_vz) * 100.0;
      Double_t sim_Ebeam_at_z = elossH2.GetEnergy(Ebeam_buff, sim_dist3D * 10.0);

      Double_t sim_theta_lab = TMath::ACos(sim_pz / sim_p);
      // FIX: kine_2b devuelve theta_CM en grados directamente
      auto [sim_ex, sim_theta_cm] = kine_2b(m_C16, m_p, m_d, m_C15, sim_Ebeam_at_z, sim_theta_lab, sim_ke);

      // Rellenamos θ_CM del generado (denominador)
      if (sim_vz * 100 > 2.0 && sim_vz * 100 < 60.0 && sim_ke > 5.0 && sim_ke < 20.0) {
         hGen_tCM->Fill(sim_theta_cm);
         h_simEx->Fill(sim_ex);
      }
   }
   cout << "N_gen_total (deuterons) = " << N_gen_total << endl;

   // =========================================================
   // --- CÁLCULO DE EFICIENCIA ---
   // =========================================================
   // ε(θ_CM) = N_reco(θ_CM) / N_gen(θ_CM)
   // El denominador ya tiene el jacobiano implícito (PolarUniform → no plano en θ_CM)
   // Al dividir bin a bin se cancela.
   //
   // FIX: NO escalamos hGen_tCM a hReco_tCM antes de dividir.
   // La eficiencia absoluta requiere que hGen_tCM esté en unidades de
   // "eventos generados por bin", que ya lo está.
   //
   // Sin embargo, como los datos reales tienen N_ev != N_sim, necesitamos
   // escalar el denominador por N_reco_sim / N_gen para obtener ε ∈ [0,1].
   // Lo hacemos con un TH1 clonado normalizado:

   // =========================================================
   // --- CÁLCULO DE EFICIENCIA CON TEfficiency ---
   // =========================================================

   if (!TEfficiency::CheckConsistency(*hReco_tCM, *hGen_tCM)) {
      std::cerr << "Histograms not consistent for TEfficiency" << std::endl;
      return;
   }

   TEfficiency *pEff = new TEfficiency(*hReco_tCM, *hGen_tCM);
   pEff->SetStatisticOption(TEfficiency::kFCP); // Clopper-Pearson (recommended)
   pEff->SetTitle("Efficiency vs #theta_{CM};#theta_{CM} (deg);#varepsilon");

   // Guardar en archivo
   TFile *fEff = new TFile("efficiency_tCM.root", "RECREATE");
   pEff->Write("efficiency_tCM");
   hReco_tCM->Write("hReco_tCM");
   hGen_tCM->Write("hGen_tCM");
   fEff->Close();

   // =========================================================
   // --- PLOTS ---
   // =========================================================

   // Canvas 1: Espectro Ex
   TCanvas *c_ExEner = new TCanvas("ExEner", "Excited Energy spectra", 1200, 800);
   h_simEx->SetLineColor(kRed);
   h_simEx->Draw("HIST");
   hexCorr2->Draw("HIST same");
   auto leg1 = new TLegend(0.6, 0.7, 0.9, 0.9);
   leg1->AddEntry(hexCorr2, "Data (corrected)", "l");
   leg1->AddEntry(h_simEx, "Sim (raw)", "l");
   leg1->Draw();
   c_ExEner->Update();

   // Canvas 2: θ_CM generado vs reconstruido
   TCanvas *c_eff_check = new TCanvas("c_eff_check", "Efficiency check", 1400, 600);
   c_eff_check->Divide(3, 1);

   c_eff_check->cd(1);
   hGen_tCM->SetLineColor(kRed);
   hGen_tCM->GetYaxis()->SetMaxDigits(3);
   hGen_tCM->SetTitle("Generado (sin detector)");
   gPad->SetTopMargin(0.15);
   hGen_tCM->Draw("HIST");

   c_eff_check->cd(2);
   hReco_tCM->SetLineColor(kBlue);
   hReco_tCM->GetYaxis()->SetMaxDigits(3);
   hReco_tCM->SetTitle("Reconstruido (SPYRAL, con cortes)");
   gPad->SetTopMargin(0.15);
   hReco_tCM->Draw("HIST");

   // Canvas 3: Eficiencia
   c_eff_check->cd(3);
   pEff->SetTitle("Efficiency vs #theta_{CM};#theta_{CM} (deg);#varepsilon, E_ej cut: 5-15 MeV");
   pEff->Draw("AP"); // A = axis, P = points

   c_eff_check->Update();
   c_eff_check->SaveAs("efficiency_thetaCM_KE_lt_15.png");
}