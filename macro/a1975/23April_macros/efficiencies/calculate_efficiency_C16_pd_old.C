#include <fstream>
#include <iostream>

double Ebin_max = 9.;
double Ebin_min = -2.;
int NumberBins = 100; // 140
int NumberBinsAux = 200;

Double_t omega(Double_t x, Double_t y, Double_t z)
{
   return sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}

std::tuple<double, double>
kine_2b(Double_t m1, Double_t m2, Double_t m3, Double_t m4, Double_t K_proj, Double_t thetalab, Double_t K_eject)
{
   // in this definition: m1(projectile); m2(target); m3(ejectile); and m4(recoil);
   double Et1 = K_proj + m1;
   double Et2 = m2;
   double Et3 = K_eject + m3;
   double Et4 = Et1 + Et2 - Et3;
   double m4_ex, Ex, theta_cm;
   double s, t, u; //---Mandelstam variables

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

   theta_cm = theta_cm * TMath::RadToDeg();

   /*cout << "m1 = " << m1 << " MeV/c^2" << "m2" << m2 << " MeV/c^2" << "m3" << m3 << " MeV/c^2" << "m4" << m4 << "
   MeV/c^2" << endl; cout << "K_proj = " << K_proj << " MeV" << "thetalab = " << thetalab << " radianes" << "K_eject = "
   << K_eject << " MeV" << endl;*/
   return std::make_tuple(Ex, theta_cm);
}

void GetEnergy(Double_t M, Double_t IZ, Double_t BRO, Double_t &E);

//-------------------------------main function---------------------------------------
void calculate_efficiency_C16_pd()
{
   bool guardar_en_pdf = false; // ← cambia a false si quieres solo verlos en pantalla
   gROOT->ProcessLine(".X /home/georgina/fair_install/ATTPCROOTv2/macro/a1975/myStyle.C");
   gStyle->SetTitleAlign(23); // Centrado horizontal (2) y alineado arriba (3)
   gStyle->SetTitleX(0.5);    // Posición X en el centro (coordenadas de 0 a 1)

   // Activar modo batch si estás guardando en PDF
   if (guardar_en_pdf) {
      gROOT->SetBatch(kTRUE); // ← esto evita que se abran ventanas
   } else {
      gROOT->SetBatch(kFALSE); // ← esto permite ver los canvas en pantalla
   }

   TH2F *Ang_Ener_Corr = new TH2F("Ang_Ener_Corr", "Ang_Ener_Corr", 720, 10, 45, 1000, 0, 60.0);
   TH2F *Ang_Ener_Cal = new TH2F("Ang_Ener_Cal", "Ang_Ener_Cal_tilt", 720, 10, 45, 1000, 0, 60.0);

   TH2F *Ebeam_test = new TH2F("Ebeam_test", "Ebeam_test", 1000, -2, 10, 60, 0, 300);

   TH2F *ELossvsBrho = new TH2F("ELossvsBrho", "ELossvsBrho", 4000, 0, 25000, 1000, 0, 4);
   TH2F *dedxvsBrho = new TH2F("dedxvsBrho", "dedxvsBrho", 4000, 0, 10000, 1000, 0, 4);
   TH2F *hVxVy = new TH2F("hVxVy", "hVxVy", 1000, 0, 4, 1000, 0, 4);
   auto *QvsEb = new TH2F("QvsEb", "QvsEb", 1000, -5, 15, 100, 0, 300);
   auto *QvsZpos = new TH2F("QvsZpos", "QvsZpos", 1000, -10, 50, 200, -100, 100);

   auto *hexCorr = new TH1F("hexCorr", "", NumberBins, Ebin_min, Ebin_max);
   auto *hexCorr2 = new TH1F("hexCorr2", "C16(p,d)", NumberBins, Ebin_min, Ebin_max);

   auto *AngDistr = new TH1F("Ang_Distr", "Ang_Distr", 90, 0, 90);
   auto *AngDistrCM = new TH1F("Ang_Distr_CM", "Ang_Distr_CM", 180, 0, 180);
   auto *ExvsZpos = new TH2F("ExvsZpos", "ExvsZpos", 1000, -5, 10, 200, -5, 80);
   // auto *ExvsTrackLength = new TH2F("ExvsTrackLength", "ExvsTrackLength", 1000, -5, 15, 200, -20, 150);
   auto *ExCorrvsZpos = new TH2F("ExCorrvsZpos", "ExCorrvsZpos", 1000, -5, 10, 200, -5, 80);
   auto *KineticEnergy = new TH1F("KineticEnergy", "KineticEnergy", 100, 0, 75);

   TH1F *h_PS_1n = new TH1F("h_PS_1n", "h_PS_1n", NumberBins, Ebin_min, Ebin_max);

   /*auto *hredchi2 = new TH1F("redchi2", "redchi2", 1000, 0, 0.0001);
   auto *hbredchi2 = new TH1F("bredchi2", "bredchi2", 1000, 0, 5);*/

   auto *hexvstheta_CM = new TH2F("hexVStheta_CM", "hexVStheta_CM", NumberBins, Ebin_min, Ebin_max, 180, 0, 180);
   auto *hexvstheta_lab = new TH2F("hexVStheta_lab", "hexVStheta_lab", NumberBins, Ebin_min, Ebin_max, 50, 0, 50);
   //--------------------------------------------------------------------------------------

   auto *hSim_thetaCM = new TH1F("hSim_thetaCM", "Sim #theta_{CM}", 180, 0, 180); // bins de 1 grados
   auto *hDat_thetaCM = new TH1F("hDat_thetaCM", "Data #theta_{CM}", 180, 0, 180);
   auto *h_simEx = new TH1F("h_simEx", "simEx", NumberBins, Ebin_min, Ebin_max);

   //--------------------------------------------------------------------------------------
   Double_t nc_tot[200];
   Double_t nc_PS_1n[200];
   Double_t x[200];
   Int_t nbins;

   // Some useful transformation constants.
   Double_t u_to_MeV = 931.49410242; // MeV/u (CODATA 2018)
   Double_t Brho_to_p = 1.602176634E-19;

   // Some masses that may be useful for the experiment
   Double_t m_p = 938.272076 / 1.0;     // masa nuclear del protón
   Double_t m_d = 1875.612931 / 1.0;    // masa nuclear del deuterón
   Double_t m_C15 = 13979.218707 / 1.0; // masa nuclear del 15C
   Double_t m_C16 = 14914.533798 / 1.0; // masa nuclear del 16C

   // Correct nuclear masses (MeV/c^2)
   // = atomic mass (u) * 931.494 - Z * 0.511 (electron mass)

   // Verify Q-value
   double Q = m_C16 + m_p - m_d - m_C15;
   cout << "Q-value = " << Q << " MeV " << endl;

   // Beam and target parameters.
   // Double_t Ebeam_buff = 186.39; //11.5 * 16; // MeV, energía del haz en el buffer gas
   Double_t Ebeam_buff = 11.5 * 16; // MeV, energía del haz en el buffer gas
   Double_t m_b = m_d;
   Double_t m_B = m_C15;

   // Ejectile parameters:deuterium
   int A_ej = 2;
   int Z_ej = 1;
   Double_t m_ej = m_d;

   // -----------------------------KINEMATICS FOR DIFFERENT EXCITATION ENERGIES

   std::vector<TString> filenames;
   double densityH2 = 3.3084e-5; // g/cm³
   AtTools::AtELossCATIMA elossH2(densityH2);
   double mass{16.0147};                        // Mass of C16 in u
   elossH2.SetMaterial(catima::Material(1, 1)); // Set material to H2
   elossH2.SetProjectile(16, 6, mass);          // Set projectile to proton

   std::set<int> excluded = {111, 121, 148, 149};
   TChain *chain = new TChain("parquettree");

   for (int i = 104; i <= 186; i++) {
      if (excluded.count(i))
         continue;
      char name[64];
      std::snprintf(name, sizeof(name), "run_%04d_2H.root", i);
      filenames.push_back(name);
      chain->Add(
         ("/home/georgina/engine_ExUniform_pd/engine_ExUniform_pd/InterpSolver/InterpSolverRoot/" + std::string(name))
            .c_str());
   }

   cout << "Total entries across all the reconstruction files: " << chain->GetEntries() << endl;

   // Loop over the chain entries of the reconstructed simulation
   for (auto filename : filenames) {
      TFile *runFile = new TFile(
         "/home/georgina/engine_ExUniform_pd/engine_ExUniform_pd/InterpSolver/InterpSolverRoot/" + filename, "R");
      TTree *Tphysics = (TTree *)runFile->Get("parquettree");

      Double_t theta{};
      Double_t phi{};
      Double_t Brho{};
      Double_t redchi{};
      Double_t zPos{}; // zPos is in meters, I will convert it to cm when filling the histograms
      Double_t ke{};

      Tphysics->SetBranchAddress("polar", &theta);   // rad
      Tphysics->SetBranchAddress("azimuthal", &phi); // rad
      Tphysics->SetBranchAddress("brho", &Brho);     // T·m
      Tphysics->SetBranchAddress("redchisq", &redchi);
      Tphysics->SetBranchAddress("vertex_z",
                                 &zPos);     // zPos is in meters, I will convert it to cm when filling the histograms
      Tphysics->SetBranchAddress("ke", &ke); // ke is in MeV

      Double_t vx_pos{},
         vy_pos{}; // vx_pos and vy_pos are in meters, I will convert them to cm when filling the histograms
      Tphysics->SetBranchAddress("vertex_x", &vx_pos);
      Tphysics->SetBranchAddress("vertex_y", &vy_pos);

      for (int i = 0; i < Tphysics->GetEntries(); i++) {

         Tphysics->GetEntry(i);

         Double_t p_ej = Brho * Z_ej * 2.99792458 / 10 * 1000;
         Double_t E_ej = (TMath::Sqrt(p_ej * p_ej + m_ej * m_ej) - m_ej);

         double dist3D = TMath::Sqrt(vx_pos * vx_pos + vy_pos * vy_pos + zPos * zPos) * 100.0; // mm
         // Double_t Ebeam_at_z = elossH2.GetEnergy(Ebeam_buff, zPos * 100); //
         Double_t Ebeam_at_z =
            elossH2.GetEnergy(Ebeam_buff, dist3D); // Usar la distancia 3D para la corrección de energía, MeV/mm

         double theta_lab_corr = theta; // Por ahora sin corrección de ángulo, solo para probar la implementación de la
                                        // corrección de energía en el cálculo de Ex y theta_cm
         auto [ex_energy_corr, theta_cm_corr] = kine_2b(m_C16, m_p, m_b, m_B, Ebeam_at_z, theta_lab_corr, E_ej);

         // Corrección cinemática
         double kethe = 13.;
         double theta_lab_corr_tilt =
            (theta -
             (2.0 * TMath::Pi() / 4000) * (E_ej - kethe)); // theta: rad; theta_lab_corr: rad; E_ej-kethe: MeV 29.5

         auto [ex_energy_corr_tilt, theta_cm_corr_tilt] = kine_2b(m_C16, m_p, m_b, m_B, Ebeam_at_z, theta_lab_corr_tilt,
                                                                  E_ej); // energies: MeV, angles: radians

         // Fill uncorrected histogram
         ExvsZpos->Fill(ex_energy_corr, zPos * 100.0); // MeV, cm
         KineticEnergy->Fill(E_ej); // Usar energía calibrada para el histograma de energía cinética

         Ang_Ener_Corr->Fill(theta_lab_corr * TMath::RadToDeg(),
                             E_ej); // theta lab!! -> I still have to implement the correction of catima?
         Ang_Ener_Cal->Fill(theta_lab_corr_tilt * TMath::RadToDeg(), E_ej); // calibrado

         double theta_deg = theta * TMath::RadToDeg();

         // Fill corrected histogram
         if (zPos * 100 > 2.0 && zPos * 100 < 60.0 && E_ej < 14.0) { // cm y MeV (zPos is in meters, E_ej is in MeV)
            ExCorrvsZpos->Fill(ex_energy_corr, zPos * 100.0);
            hexCorr->Fill(ex_energy_corr_tilt);
            hexCorr2->Fill(ex_energy_corr_tilt);
            hDat_thetaCM->Fill(theta_cm_corr_tilt);

            // Histograms
            Double_t vx = TMath::Sin(theta) * TMath::Sqrt(ke);
            Double_t vy = TMath::Cos(theta) * TMath::Sqrt(ke);

            hVxVy->Fill(vx, vy);

            AngDistr->Fill(theta * TMath::RadToDeg()); // degrees
            AngDistrCM->Fill(theta_cm_corr);           // degrees
            hexvstheta_CM->Fill(ex_energy_corr,
                                theta_cm_corr); // theta_cm_corr is already in degrees, ex_energy_corr is in MeV
            hexvstheta_lab->Fill(
               ex_energy_corr,
               theta_lab_corr *
                  TMath::RadToDeg()); // theta_lab_corr_tilt is in radians, convert to degrees for the histogram
            // ExvsTrackLength->Fill(ex_energy_corr, arclength);
         }

         // tEvents->Fill();
      } // events
   }    // Files

   //---------------- SIMULACION SIN DETECTOR ----------------//
   TFile *fSim = new TFile(
      "/home/georgina/fair_install/ATTPCROOTv2/macro/a1975/23April_macros/efficiencies/rawSim/output_16Cpd.root",
      "READ");
   TTree *tSim = (TTree *)fSim->Get("kinematics");

   Long64_t sim_Z, sim_A, sim_event;
   Double_t sim_px, sim_py, sim_pz, sim_energy;
   Double_t sim_vx{}, sim_vy{}, sim_vz{};

   tSim->SetBranchAddress("event", &sim_event);
   tSim->SetBranchAddress("Z", &sim_Z);
   tSim->SetBranchAddress("A", &sim_A);
   tSim->SetBranchAddress("px", &sim_px);
   tSim->SetBranchAddress("py", &sim_py);
   tSim->SetBranchAddress("pz", &sim_pz);
   tSim->SetBranchAddress("energy", &sim_energy);
   tSim->SetBranchAddress("vertex_x", &sim_vx);
   tSim->SetBranchAddress("vertex_y", &sim_vy);
   tSim->SetBranchAddress("vertex_z", &sim_vz);

   cout << "Total simulation entries: " << tSim->GetEntries() << endl;
   for (int i = 0; i < tSim->GetEntries(); i++) {
      tSim->GetEntry(i);

      // Solo el deuterón ejectil (índice 2: Z==1, A==2)
      // Excluimos el target (Z==1, A==1)
      if (!(sim_Z == 1 && sim_A == 2))
         continue;

      // Calcular energía cinética del ejectil
      Double_t sim_p = TMath::Sqrt(sim_px * sim_px + sim_py * sim_py + sim_pz * sim_pz);
      Double_t sim_ke = sim_energy - m_d; // energía cinética = energía total - masa

      double sim_dist3D = TMath::Sqrt(sim_vx * sim_vx + sim_vy * sim_vy + sim_vz * sim_vz) * 100.0; // mm
      Double_t sim_Ebeam_at_z = elossH2.GetEnergy(Ebeam_buff, sim_dist3D);

      // Calcular theta_lab
      Double_t sim_theta_lab = TMath::ACos(sim_pz / sim_p); // radianes

      // Calcular theta_cm con kine_2b (misma función que usas con datos)
      auto [sim_ex, sim_theta_cm] = kine_2b(m_C16, m_p, m_d, m_C15, sim_Ebeam_at_z, sim_theta_lab, sim_ke);

      hSim_thetaCM->Fill(sim_theta_cm);
      h_simEx->Fill(sim_ex);
   }
   //--------------------------------------------------------------------------------------------
   TCanvas *c_ExEner = new TCanvas("ExEner", "Corrected Excited Energy spectra", 1200, 800);
   c_ExEner->cd();
   h_simEx->GetXaxis()->SetTitle("Excitation Energy (MeV)");
   h_simEx->GetYaxis()->SetTitle("Counts");
   h_simEx->SetLineColor(kRed);
   h_simEx->Draw("hist");
   hexCorr2->Draw("E1 same"); // E1 para mostrar errores

   auto leg1 = new TLegend(0.6, 0.7, 0.9, 0.9);
   leg1->AddEntry(hexCorr2, "Data (corrected)", "l");
   leg1->AddEntry(h_simEx, "Sim (raw)", "l");
   leg1->Draw();
   c_ExEner->Update();

   TCanvas *c_combined = new TCanvas("c_combined", "Combined Ang Energy", 1200, 800);
   c_combined->Divide(2, 1); // 2 columns, 2 rows

   // --- PAD 1: Ang_Ener_Corr (uncorrected kinematics) ---
   c_combined->cd(1);
   Ang_Ener_Corr->Draw("col");
   Ang_Ener_Corr->GetXaxis()->SetTitle("#theta_{lab} (deg)");
   Ang_Ener_Corr->GetYaxis()->SetTitle("Kinetic Energy (MeV)");

   /*for (size_t i = 0; i < graphs.size(); i++) {
      if(graphs[i] == nullptr || graphs[i]->GetN() == 0) continue;
      graphs[i]->Draw("L SAME");
   }
   auto leg0 = new TLegend(0.45, 0.7, 0.9, 0.9);
   for (size_t i = 0; i < graphs.size(); i++) {
      if(graphs[i] == nullptr) continue;
      leg0->AddEntry(graphs[i], labels[i].c_str(), "l");
   }
   leg0->Draw();*/

   // --- PAD 2: Ang_Ener_Cal (corrected angle kinematics) ---
   c_combined->cd(2);
   Ang_Ener_Cal->Draw("col");
   Ang_Ener_Cal->GetXaxis()->SetTitle("#theta_{lab} (deg)");
   Ang_Ener_Cal->GetYaxis()->SetTitle("Kinetic Energy (MeV)");

   /*for (size_t i = 0; i < graphs.size(); i++) {
      if(graphs[i] == nullptr || graphs[i]->GetN() == 0) continue;
      graphs[i]->Draw("L SAME");
   }
   auto leg2 = new TLegend(0.45, 0.7, 0.9, 0.9);
   for (size_t i = 0; i < graphs.size(); i++) {
      if(graphs[i] == nullptr) continue;
      leg2->AddEntry(graphs[i], labels[i].c_str(), "l");
   }
   leg2->Draw();*/

   c_combined->Update();

   TCanvas *c_AngDistr = new TCanvas("AngDistr", "Angular Distribution", 1200, 600);
   c_AngDistr->cd();
   c_AngDistr->Divide(2, 1);
   c_AngDistr->cd(1);
   c_AngDistr->SetTitle("Angular distribution in lab frame");
   AngDistr->GetXaxis()->SetTitle("#theta_{lab} (deg)");
   AngDistr->GetYaxis()->SetTitle("Counts");
   AngDistr->Sumw2();
   gPad->SetTopMargin(0.15);
   AngDistr->GetYaxis()->SetMaxDigits(3);
   AngDistr->Draw("E1");

   c_AngDistr->cd(2);
   c_AngDistr->SetTitle("Angular distribution in center-of-mass frame");
   AngDistrCM->GetXaxis()->SetTitle("#theta_{CM} (deg)");
   AngDistrCM->GetYaxis()->SetTitle("#frac{d#sigma}{d#Omega} (a.u.)");
   AngDistrCM->Sumw2();
   gPad->SetTopMargin(0.15);
   AngDistrCM->GetYaxis()->SetMaxDigits(3);
   AngDistrCM->Draw("E1");

   TCanvas *kin = new TCanvas("kin", "kin", 1200, 800);
   KineticEnergy->Sumw2();
   kin->cd();
   KineticEnergy->Draw("E1");

   //---------------- EFICIENCIA ----------------//
   TGraphErrors *gEfficiency = new TGraphErrors();
   gEfficiency->SetTitle("Efficiency; #theta_{CM} (#circ); #epsilon(#theta_{CM})");

   int nSim_total = tSim->GetEntries() / 4; // 4 particulas por evento, raw simulation sin detector, entonces el número
                                            // total de eventos es el número de entradas dividido por 4

   /*cout << "Total simulated events (raw, no detector, 4 tree events, 1 full event): " << nSim_total << endl;
   cout << "hDat_thetaCM bins: " << hDat_thetaCM->GetNbinsX() << endl;
   cout << "hSim_thetaCM bins: " << hSim_thetaCM->GetNbinsX() << endl;

   if (hDat_thetaCM->GetNbinsX() != hSim_thetaCM->GetNbinsX()) {
      cout << "Error: Number of bins in data and simulation histograms do not match!" << endl;
      return;
   }

   for (int b = 1; b <= hDat_thetaCM->GetNbinsX(); b++) {
      double theta_center = hDat_thetaCM->GetBinCenter(b);
      double N_det = hDat_thetaCM->GetBinContent(b);//post detector
      double N_gen = hSim_thetaCM->GetBinContent(b); ///pre detector (simulación sin detector)

      if (N_gen <= 0) continue;
      if (N_det > N_gen) {
       cout << "Warning: bin " << b << " (theta_cm = " << theta_center
            << " deg), N_det=" << N_det << " > N_gen=" << N_gen << endl;
       continue; // ok si ya sabes que vas a cortar theta > 170°
   }

      double eff = N_det / N_gen;
      double eff_err = TMath::Sqrt(eff * (1 - eff) / N_gen); // error sobre N original

      int point = gEfficiency->GetN();
      gEfficiency->SetPoint(point, theta_center, eff);
      gEfficiency->SetPointError(point, 0, eff_err);
   }*/

   TCanvas *c_eff_check = new TCanvas("c_eff_check", "Efficiency check", 1200, 600);
   c_eff_check->Divide(2, 1);

   c_eff_check->cd(1);
   hSim_thetaCM->SetLineColor(kRed);
   hSim_thetaCM->GetXaxis()->SetTitle("#theta_{CM} (#circ)");
   hSim_thetaCM->GetYaxis()->SetTitle("Counts");
   hSim_thetaCM->GetYaxis()->SetMaxDigits(3);
   hSim_thetaCM->SetTitle("Simulacion sin detector");
   gPad->SetTopMargin(0.15);
   hSim_thetaCM->Draw("HIST");

   c_eff_check->cd(2);
   hDat_thetaCM->SetLineColor(kBlue);
   hDat_thetaCM->GetXaxis()->SetTitle("#theta_{CM} (#circ)");
   hDat_thetaCM->GetYaxis()->SetTitle("Counts");
   hDat_thetaCM->GetYaxis()->SetMaxDigits(3);
   hDat_thetaCM->SetTitle("Datos (con cortes)");
   gPad->SetTopMargin(0.15);
   hDat_thetaCM->Draw("E1");

   TCanvas *c_eff = new TCanvas("c_eff", "Efficiency vs #theta_{CM}", 800, 600);
   c_eff->cd();
   /*gEfficiency->SetMarkerStyle(20);
   gEfficiency->SetMarkerColor(kBlue);
   gEfficiency->SetLineColor(kBlue);
   gEfficiency->Draw("APE");
   gEfficiency->GetYaxis()->SetRangeUser(0, 0.3);
   gEfficiency->GetXaxis()->SetRangeUser(0, 60);
   gEfficiency->GetXaxis()->SetTitle("#theta_{CM} (#circ)");
   gEfficiency->GetYaxis()->SetTitle("Efficiency");*/

   // Suponiendo que hDat_thetaCM son tus datos y hSim_thetaCM es tu simulación
   TEfficiency *pEff = new TEfficiency(*hDat_thetaCM, *hSim_thetaCM);

   pEff->SetMarkerStyle(20);
   pEff->SetMarkerColor(kBlue);
   pEff->SetLineColor(kBlue);

   // "AP" -> A: Dibuja los ejes alrededor, P: Dibuja los puntos con sus errores binomiales corregidos
   pEff->Draw("AP");

   // Guardar en ROOT file para usar en angular_distribution.C
   TFile *fEff = new TFile("efficiency.root", "RECREATE");
   pEff->Write("pEff");
   fEff->Close();
}

void GetEnergy(Double_t M, Double_t IZ, Double_t BRO, Double_t &E)
{

   // Energy per nucleon
   Float_t AM = 931.5;
   Float_t X = BRO / 0.1439 * IZ / M;
   X = pow(X, 2);
   X = 2. * AM * X;
   X = X + pow(AM, 2);
   E = TMath::Sqrt(X) - AM;
}
