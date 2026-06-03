#include <fstream>
#include <iostream>
#include <set>

double Ebin_max = 9.;
double Ebin_min = -1.;
int NumberBins = 100; // 140

// ── Función auxiliar de cinemática (igual que en la macro principal) ──────────
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

// ── Función principal ─────────────────────────────────────────────────────────
void simulationCorrelations_Exuniform_check()
{
   // ── Masas (MeV/c²) — mismos valores que en la macro principal ────────────
   Double_t m_p = 938.272076;
   Double_t m_d = 1875.612931;
   Double_t m_C15 = 13979.218707;
   Double_t m_C16 = 14914.533798;

   double Q = m_C16 + m_p - m_d - m_C15;
   cout << "Q-value = " << Q << " MeV " << endl;

   // Beam and target parameters.
   Double_t Ebeam_buff = 11.5 * 16; // MeV, energía del haz en el buffer gas
   Double_t m_b = m_d;
   Double_t m_B = m_C15;

   // Ejectile parameters:deuterium
   int A_ej = 2;
   int Z_ej = 1;
   Double_t m_ej = m_d;

   // Parámetros del fit (se rellenarán después del primer loop)
   double parA = 0.0;
   double parB = 0.0;

   // ── Curvas teóricas: GS para 3 energías ──────────────────────────────────
   struct KineFile {
      TString label;
      TString path;
      int color;
   };

   std::vector<KineFile> kineFiles = {
      //{"11.0 AMeV",
      //"/home/georgina/fair_install/ATTPCROOTv2/macro/Kinematics/Decay_kinematics/16C_pd_15C_GS_11_26May.txt", kBlue},
      {"GS kine line (184.0 MeV)",
       "/home/georgina/fair_install/ATTPCROOTv2/macro/Kinematics/Decay_kinematics/16C_pd_15C_GS_11_5_26May.txt", kRed},
      //{"Maximum Eloss GS (171.874 MeV)",
      //"/home/georgina/fair_install/ATTPCROOTv2/macro/Kinematics/Decay_kinematics/16C_pd_15C_GS_171_874MeV.txt",
      // kOrange+7},
      //{"12.0 AMeV",
      //"/home/georgina/fair_install/ATTPCROOTv2/macro/Kinematics/Decay_kinematics/16C_pd_15C_GS_12_26May.txt",
      // kViolet+2},
      //{"11.65 AMeV",
      //"/home/georgina/fair_install/ATTPCROOTv2/macro/Kinematics/Decay_kinematics/16C_pd_15C_GS_11_65_26May.txt",
      // kGreen+2}
   };

   std::vector<TGraph *> gKine;

   for (auto &kf : kineFiles) {
      std::ifstream fin(kf.path.Data());
      if (!fin.is_open()) {
         std::cout << "[ERROR] No se pudo abrir: " << kf.path << std::endl;
         gKine.push_back(nullptr);
         continue;
      }
      // Saltar cabecera
      std::string header;
      std::getline(fin, header);

      std::vector<Double_t> theta, ener;
      Double_t tCMS, tLabRec, eLabRec, tLabSca, eLabSca;
      while (fin >> tCMS >> tLabRec >> eLabRec >> tLabSca >> eLabSca) {
         theta.push_back(tLabRec);
         ener.push_back(eLabRec);
      }

      if (theta.empty()) {
         std::cout << "[WARNING] Fichero vacío o formato incorrecto: " << kf.path << std::endl;
         gKine.push_back(nullptr);
         continue;
      }

      TGraph *g = new TGraph(theta.size(), theta.data(), ener.data());
      g->SetLineColor(kf.color);
      g->SetLineWidth(4);
      g->SetTitle(kf.label);
      gKine.push_back(g);
      std::cout << " -> " << kf.label << " cargada (" << theta.size() << " puntos)." << std::endl;
   }

   //....................
   // Estructura para guardar la verdad simulada
   struct TruthData {
      double K_real;
      double theta_real; // en grados
   };

   std::map<Long64_t, TruthData> truthMap;
   // ==========================================
   // 1. LEER LA SIMULACIÓN Y LLENAR EL MAPA
   // ==========================================
   std::cout << "Leyendo archivo de simulacion (.root)..." << std::endl;
   TFile *file_sim = TFile::Open("output_16Cpd_ExUniform-2_9MeV.root");

   // NECESITAMOS EL NOMBRE DEL ÁRBOL AQUÍ:
   TTree *tree_sim = (TTree *)file_sim->Get("kinematics");

   // Variables para leer las ramas de la simulación
   Long64_t sim_eventID, sim_Z, sim_A;
   Double_t sim_px, sim_py, sim_pz;

   // NECESITAMOS CONECTAR LAS RAMAS AQUÍ:
   tree_sim->SetBranchAddress("event", &sim_eventID);
   tree_sim->SetBranchAddress("Z", &sim_Z);
   tree_sim->SetBranchAddress("A", &sim_A);
   tree_sim->SetBranchAddress("px", &sim_px);
   tree_sim->SetBranchAddress("py", &sim_py);
   tree_sim->SetBranchAddress("pz", &sim_pz);

   for (Long64_t i = 0; i < tree_sim->GetEntries(); i++) {
      tree_sim->GetEntry(i);

      // Filtrar al deuterón
      if (sim_Z == 1 && sim_A == 2) {
         double p = std::sqrt(sim_px * sim_px + sim_py * sim_py + sim_pz * sim_pz);

         TruthData data;
         data.K_real = std::sqrt(p * p + m_d * m_d) - m_d;
         data.theta_real = std::acos(sim_pz / p) * TMath::RadToDeg();

         truthMap[sim_eventID] = data;
      }
   }
   file_sim->Close();
   std::cout << "Se cargaron " << truthMap.size() << " deuterones desde la simulacion." << std::endl;
   //....................................
   // Histogramas de Residuos (Correlaciones)
   auto *h_dK_vs_theta =
      new TH2F("h_dK_vs_theta", "#Delta K vs #theta_{rec};#theta_{rec} (deg);K_{rec} - K_{real} (MeV)", 100, 5, 70, 100,
               -20, 20);
   auto *h_dK_vs_K =
      new TH2F("h_dK_vs_K", "#Delta K vs K_{rec};K_{rec} (MeV);K_{rec} - K_{real} (MeV)", 100, 5, 70, 100, -20, 20);
   auto *h_dTheta_vs_K =
      new TH2F("h_dTheta_vs_K", "#Delta#theta vs K_{rec};K_{rec} (MeV);#theta_{rec} - #theta_{real} (deg)", 100, 5, 70,
               100, -20, 20);
   auto *h_dTheta_vs_theta = new TH2F(
      "h_dTheta_vs_theta", "#Delta#theta vs #theta_{rec};#theta_{rec} (deg);#theta_{rec} - #theta_{real} (deg)", 100, 5,
      70, 100, -20, 20);

   auto *h_dK_vs_theta_real =
      new TH2F("h_dK_vs_theta_real", "#Delta K vs #theta_{real};#theta_{real} (deg);K_{rec} - K_{real} (MeV)", 100, 5,
               70, 100, -20, 20);
   auto *h_dK_vs_K_real = new TH2F("h_dK_vs_K_real", "#Delta K vs K_{real};K_{real} (MeV);K_{rec} - K_{real} (MeV)",
                                   100, 5, 70, 100, -20, 20);
   auto *h_dTheta_vs_K_real =
      new TH2F("h_dTheta_vs_K_real", "#Delta#theta vs K_{real};K_{real} (MeV);#theta_{rec} - #theta_{real} (deg)", 100,
               5, 70, 100, -20, 20);
   auto *h_dTheta_vs_theta_real = new TH2F(
      "h_dTheta_vs_theta_real", "#Delta#theta vs #theta_{real};#theta_{real} (deg);#theta_{rec} - #theta_{real} (deg)",
      100, 5, 70, 100, -20, 20);

   auto *h_ke_vs_Eej = new TH2F("h_ke_vs_Eej", "ke vs E_ej;E_ej (MeV);ke (MeV)", 200, 0, 70, 200, 0, 70);
   auto *h_dKE = new TH1F("h_dKE", "E_ej - ke;E_ej - ke (MeV);Counts", 200, -10, 10);

   // ── Histograma 2D: Ex vs Z ───────────────────────────────────────────
   auto *ExvsZpos = new TH2F("ExvsZpos", "", 100, -6, 10, 200, -5, 100);
   auto *hex = new TH1F("hex", "C16(p,d)", NumberBins, -7.0, 10.0);
   TH2F *Ang_Ener = new TH2F("Ang_Ener", "Ang_Ener", 400, 10, 50, 100, 0, 50.0);

   auto *ExvsZpos_cuts = new TH2F("ExvsZpos_cuts", "", 100, -6, 10, 200, -5, 100);
   auto *hex_cuts = new TH1F("hex_cuts", "C16(p,d)", 200, -6.0, 10.0);
   TH2F *Ang_Ener_cuts = new TH2F("Ang_Ener_cuts", "Ang_Ener_cuts", 400, 10, 50, 100, 0, 50.0);

   // ── ELoss CATIMA (igual que en la macro principal) ────────────────────────
   double densityH2 = 3.553e-5; // g/cm³
   AtTools::AtELossCATIMA elossH2(densityH2);
   elossH2.SetMaterial(catima::Material(1, 1));
   elossH2.SetProjectile(16, 6, 16.0147);

   // ── Datos reconstruidos ──────────────────────────────────────────────────
   std::vector<TString> filenames;
   TChain *chain = new TChain("parquettree"); // parquettree
   for (int i = 104; i <= 186; i++) {
      char name[64];
      std::snprintf(name, sizeof(name), "run_%04d_2H.root", i);
      filenames.push_back(name);
      chain->Add(("/home/georgina/my_sim/engine_ExUniform_pd/InterpSolver_pd_root/" + std::string(name)).c_str());
   }
   std::cout << "Total entries (reconstruction): " << chain->GetEntries() << std::endl;

   Long64_t rec_eventID = 0;
   chain->SetBranchAddress("event", &rec_eventID); // Asegúrate de que "event" es el nombre correcto en tu árbol
   Double_t theta{}, phi{}, Brho{}, redchi{}, zPos{}, ke{};
   Double_t vx_pos{}, vy_pos{};
   chain->SetBranchAddress("polar", &theta);
   chain->SetBranchAddress("azimuthal", &phi);
   chain->SetBranchAddress("brho", &Brho);
   chain->SetBranchAddress("redchisq", &redchi);
   chain->SetBranchAddress("vertex_z", &zPos);
   chain->SetBranchAddress("ke", &ke);
   chain->SetBranchAddress("vertex_x", &vx_pos); // m
   chain->SetBranchAddress("vertex_y", &vy_pos);

   for (Long64_t i = 0; i < chain->GetEntries(); i++) {
      chain->GetEntry(i);

      Double_t p_ej = Brho * Z_ej * 2.99792458 / 10 * 1000;
      Double_t E_ej = TMath::Sqrt(p_ej * p_ej + m_ej * m_ej) - m_ej;

      double dist3D = TMath::Sqrt(vx_pos * vx_pos + vy_pos * vy_pos + zPos * zPos) * 100.0;
      Double_t Ebeam_at_z =
         elossH2.GetEnergy(Ebeam_buff, dist3D * 10.0); // Convertir dist3D a mm para la corrección de energía

      double theta_lab_corr = theta;
      auto [ex_corr, theta_cm_corr] = kine_2b(m_C16, m_p, m_b, m_B, Ebeam_at_z, theta_lab_corr, E_ej);

      // --- FILTRO DE PROTECCIÓN (El reemplazo del cout) ---
      // Si el cálculo cinemático da un error matemático, ignoramos este evento
      if (TMath::IsNaN(ex_corr) || TMath::IsNaN(E_ej) || TMath::IsNaN(Ebeam_at_z)) {
         continue;
      }

      // --- EL MATCH MÁGICO DE RESIDUOS ---
      // Buscamos si el evento de la TPC existe en nuestro mapa de verdad simulada
      if (truthMap.find(rec_eventID) != truthMap.end()) {

         // 1. Extraemos la verdad
         double K_real = truthMap[rec_eventID].K_real;
         double theta_real_deg = truthMap[rec_eventID].theta_real;

         // 2. Definimos lo reconstruido
         double K_rec = E_ej;
         double theta_rec_deg = theta_lab_corr * TMath::RadToDeg();

         // 3. Calculamos los residuos (Reconstruido - Real)
         double dK = K_rec - K_real;
         double dTheta = theta_rec_deg - theta_real_deg;

         // 4. Llenamos los histogramas
         // (Aplicamos el mismo corte Z y Ex que usas para quedarte con el GS limpio)
         // if (zPos*100 > 2.0 && zPos*100 < 60.0 && E_ej < 20.0) {

         h_dK_vs_theta->Fill(theta_rec_deg, dK);
         h_dK_vs_K->Fill(K_rec, dK);
         h_dTheta_vs_K->Fill(K_rec, dTheta);
         h_dTheta_vs_theta->Fill(theta_rec_deg, dTheta);

         h_dK_vs_theta_real->Fill(theta_real_deg, dK);
         h_dK_vs_K_real->Fill(K_real, dK);
         h_dTheta_vs_K_real->Fill(K_real, dTheta);
         h_dTheta_vs_theta_real->Fill(theta_real_deg, dTheta);

         h_ke_vs_Eej->Fill(E_ej, ke);
         h_dKE->Fill(E_ej - ke);

         //}
      }

      hex->Fill(ex_corr);
      ExvsZpos->Fill(ex_corr, zPos * 100);

      if (ex_corr > 0) {
         Ang_Ener->Fill(theta_lab_corr * TMath::RadToDeg(), E_ej); // calibrado
      }
      // Ang_Ener->Fill(theta_lab_corr* TMath::RadToDeg(), E_ej);

      if (zPos * 100 > 2.0 && zPos * 100 < 60.0 &&
          E_ej < 20.0) { // cm y MeV, y un corte en Ex para quedarnos solo con la GS
         hex_cuts->Fill(ex_corr);
         ExvsZpos_cuts->Fill(ex_corr, zPos * 100);
         // Ang_Ener_cuts->Fill(theta_lab_corr* TMath::RadToDeg(), E_ej);
      }
   }

   // 1. Extraemos el perfil y ajustamos el POLINOMIO (pol3)
   TProfile *prof_dK_K = h_dK_vs_K->ProfileX("prof_dK_K");
   TF1 *myFit = new TF1("myFit", "pol3", 5, 60);
   prof_dK_K->Fit(myFit, "R");

   // Ya no necesitamos extraer parámetros manualmente, Eval() hace el trabajo.

   // 2. Bucle para aplicar la corrección
   for (Long64_t i = 0; i < chain->GetEntries(); i++) {
      chain->GetEntry(i);

      Double_t p_ej = Brho * Z_ej * 2.99792458 / 10 * 1000;
      Double_t E_ej = TMath::Sqrt(p_ej * p_ej + m_ej * m_ej) - m_ej;
      double theta_rec_deg = theta * TMath::RadToDeg();

      // 3. Calculamos la corrección usando el polinomio
      Double_t correction = myFit->Eval(E_ej);

      // Declaramos E_ej_corr UNA SOLA VEZ
      Double_t E_ej_corr = E_ej - correction;

      // 4. Llenamos los histogramas de la simulación
      if (truthMap.find(rec_eventID) != truthMap.end()) {
         double K_real = truthMap[rec_eventID].K_real;

         double dK_corr_res = E_ej_corr - K_real;

         h_dK_vs_theta_corr->Fill(theta_rec_deg, dK_corr_res);
         h_dK_vs_K_corr->Fill(E_ej_corr, dK_corr_res);

         // Llenamos el de cinemática solo con eventos simulados (Opción 2 de antes)
         Ang_Ener_cuts->Fill(theta_rec_deg, E_ej_corr);
      }
   }

   // ── Canvas y dibujo ───────────────────────────────────────────────────────
   /*TCanvas *c = new TCanvas("c_ZvsEx_GS", "Z vs Ex: Ground State", 900, 700);
   c->Divide(2,1);
   c->cd(1);
   ExvsZpos->GetXaxis()->SetTitle("Ex (MeV)");
   ExvsZpos->GetYaxis()->SetTitle("Z (cm)");
   ExvsZpos->SetTitle("Ex vs Z");
   ExvsZpos->Draw("COLZ");
   c->cd(2);
    hex->GetXaxis()->SetTitle("E_{x} (MeV)");
    hex->GetYaxis()->SetTitle("Counts");
    hex->SetTitle("Excitation Energy");
    hex->Draw("same");

   c->Update();

   TCanvas *c_cuts = new TCanvas("c_ZvsEx_GS_cuts", "Z vs Ex with cuts: Ground State", 900, 700);
   c_cuts->Divide(2,1);
   c_cuts->cd(1);
   ExvsZpos_cuts->GetXaxis()->SetTitle("Ex (MeV)");
   ExvsZpos_cuts->GetYaxis()->SetTitle("Z (cm)");
   ExvsZpos_cuts->SetTitle("Ex vs Z");
   ExvsZpos_cuts->Draw("COLZ");
   c_cuts->cd(2);
   hex_cuts->GetXaxis()->SetTitle("E_{x} (MeV)");
   hex_cuts->GetYaxis()->SetTitle("Counts");
   hex_cuts->SetTitle("Excitation Energy");
   hex_cuts->Draw("same");
   c_cuts->Update();*/

   // ── Canvas para los gráficos de residuos ─────────────────────────────────
   TCanvas *c_correlaciones = new TCanvas("c_corr", "Analisis de Residuos", 1200, 800);
   c_correlaciones->Divide(2, 2);

   c_correlaciones->cd(1);
   gPad->SetLogz();
   h_dK_vs_theta->Draw("colz");

   c_correlaciones->cd(2);
   gPad->SetLogz();
   h_dK_vs_K->Draw("colz");

   c_correlaciones->cd(3);
   gPad->SetLogz();
   h_dTheta_vs_K->Draw("colz");

   c_correlaciones->cd(4);
   gPad->SetLogz();
   h_dTheta_vs_theta->Draw("colz");

   c_correlaciones->Update();

   // ── Canvas para los gráficos de residuos ─────────────────────────────────
   /*TCanvas *c_correlaciones2 = new TCanvas("c_corr2", "Analisis de Residuos Real", 1200, 800);
   c_correlaciones2->Divide(2, 2);

   c_correlaciones2->cd(1);
   gPad->SetLogz();
   h_dK_vs_theta_real->Draw("colz");

   c_correlaciones2->cd(2);
   gPad->SetLogz();
   h_dK_vs_K_real->Draw("colz");

   c_correlaciones2->cd(3);
   gPad->SetLogz();
   h_dTheta_vs_K_real->Draw("colz");

   c_correlaciones2->cd(4);
   gPad->SetLogz();
   h_dTheta_vs_theta_real->Draw("colz");

   c_correlaciones2->Update();*/

   /*TCanvas *c_ke = new TCanvas("c_ke", "ke vs E_ej", 1200, 500);
   c_ke->Divide(2, 1);

   c_ke->cd(1);
   gPad->SetLogz();
   h_ke_vs_Eej->Draw("colz");
   // Línea diagonal y=x para referencia
   TLine *diag = new TLine(0, 0, 70, 70);
   diag->SetLineColor(kRed);
   diag->SetLineWidth(2);
   diag->Draw("same");

   c_ke->cd(2);
   h_dKE->Draw();
   // Línea en 0 para referencia
   TLine *zero = new TLine(0, 0, 0, h_dKE->GetMaximum());
   zero->SetLineColor(kRed);
   zero->SetLineWidth(2);
   zero->Draw("same");

   c_ke->Update();*/

   TCanvas *c_perfil = new TCanvas("c_perfil", "Perfil del sesgo", 1200, 500);
   c_perfil->Divide(2, 1);

   c_perfil->cd(1);
   TProfile *prof_dK_theta = h_dK_vs_theta->ProfileX("prof_dK_theta");
   prof_dK_theta->SetTitle("#Delta K vs #theta_{rec} - Perfil;#theta_{rec} (deg);<K_{rec} - K_{real}> (MeV)");
   prof_dK_theta->SetMarkerStyle(20);
   prof_dK_theta->SetMarkerColor(kRed);
   prof_dK_theta->Draw();

   c_perfil->cd(2);
   // Ya NO usamos "TProfile *" ni "TF1 *" porque ya existen
   prof_dK_K->SetTitle("#Delta K vs K_{rec} - Perfil;K_{rec} (MeV);<K_{rec} - K_{real}> (MeV)");
   prof_dK_K->SetMarkerStyle(20);
   prof_dK_K->SetMarkerColor(kBlue);
   prof_dK_K->Draw();

   myFit->SetLineColor(kRed);
   myFit->Draw("same");

   c_perfil->Update();

   TCanvas *c_corr_test = new TCanvas("c_corr_test", "Residuos corregidos", 1200, 500);
   c_corr_test->Divide(2, 1);

   c_corr_test->cd(1);
   gPad->SetLogz();
   h_dK_vs_theta_corr->Draw("colz");

   c_corr_test->cd(2);
   gPad->SetLogz();
   h_dK_vs_K_corr->Draw("colz");

   c_corr_test->Update();

   TCanvas *c_ang_ener_cuts = new TCanvas("c_ang_ener_cuts", "Angle vs Energy with cuts", 900, 700);
   c_ang_ener_cuts->Divide(2, 1);
   c_ang_ener_cuts->cd(1);
   Ang_Ener->GetXaxis()->SetTitle("#theta_{lab} (deg)");
   Ang_Ener->GetYaxis()->SetTitle("KE (MeV)");
   Ang_Ener->SetTitle("Theta vs KE");
   Ang_Ener->Draw("colz");
   TLegend *leg2 = new TLegend(0.15, 0.65, 0.40, 0.85);
   leg2->SetBorderSize(0);
   leg2->SetFillStyle(0);
   leg2->SetTextSize(0.02);

   for (auto *g : gKine) {
      if (g) {
         g->Draw("L SAME");
         leg2->AddEntry(g, g->GetTitle(), "l");
      }
   }
   leg2->Draw("same");

   c_ang_ener_cuts->cd(2);
   Ang_Ener_cuts->GetXaxis()->SetTitle("#theta_{lab} (deg)");
   Ang_Ener_cuts->GetYaxis()->SetTitle("KE (MeV)");
   Ang_Ener_cuts->SetTitle("Theta vs KE");
   Ang_Ener_cuts->Draw("colz");

   TLegend *leg3 = new TLegend(0.15, 0.65, 0.40, 0.85);
   leg3->SetBorderSize(0);
   leg3->SetFillStyle(0);
   leg3->SetTextSize(0.02);

   for (auto *g : gKine) {
      if (g) {
         g->Draw("L SAME");
         leg3->AddEntry(g, g->GetTitle(), "l");
      }
   }
   leg3->Draw("same");
   c_ang_ener_cuts->Update();

   // ── Canvas para ver los perfiles corregidos centrados en 0 ──────────────────
   TCanvas *c_perfil_corr = new TCanvas("c_perfil_corr", "Perfiles de sesgo corregidos", 1200, 500);
   c_perfil_corr->Divide(2, 1);

   // --- 1. Perfil de Delta K vs Theta (Corregido) ---
   c_perfil_corr->cd(1);
   TProfile *prof_dK_theta_corr = h_dK_vs_theta_corr->ProfileX("prof_dK_theta_corr");
   prof_dK_theta_corr->SetTitle(
      "Sesgo Corregido: #Delta K vs #theta_{rec};#theta_{rec} (deg);<K_{rec,corr} - K_{real}> (MeV)");
   prof_dK_theta_corr->SetMarkerStyle(20);
   prof_dK_theta_corr->SetMarkerColor(kGreen + 2); // Verde para indicar éxito/corregido

   // Forzamos el eje Y para que el 0 quede exactamente en el medio
   prof_dK_theta_corr->SetMinimum(-5);
   prof_dK_theta_corr->SetMaximum(5);
   prof_dK_theta_corr->Draw();

   // Línea de referencia en 0
   TLine *linea_cero_theta = new TLine(5, 0, 70, 0); // Ajusta el 5-70 según tus ejes de Theta
   linea_cero_theta->SetLineStyle(2);
   linea_cero_theta->SetLineColor(kBlack);
   linea_cero_theta->Draw("same");

   // --- 2. Perfil de Delta K vs K_rec (Corregido) ---
   c_perfil_corr->cd(2);
   TProfile *prof_dK_K_corr = h_dK_vs_K_corr->ProfileX("prof_dK_K_corr");
   prof_dK_K_corr->SetTitle("Sesgo Corregido: #Delta K vs K_{rec};K_{rec,corr} (MeV);<K_{rec,corr} - K_{real}> (MeV)");
   prof_dK_K_corr->SetMarkerStyle(20);
   prof_dK_K_corr->SetMarkerColor(kGreen + 2);

   // Forzamos el eje Y
   prof_dK_K_corr->SetMinimum(-5);
   prof_dK_K_corr->SetMaximum(5);
   prof_dK_K_corr->Draw();

   // Línea de referencia en 0
   TLine *linea_cero_K = new TLine(5, 0, 70, 0);
   linea_cero_K->SetLineStyle(2);
   linea_cero_K->SetLineColor(kBlack);
   linea_cero_K->Draw("same");

   c_perfil_corr->Update();
}