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
void groundStateSim_check()
{

   gStyle->SetCanvasPreferGL(kTRUE);
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
      {"Maximum Eloss GS (171.874 MeV)",
       "/home/georgina/fair_install/ATTPCROOTv2/macro/Kinematics/Decay_kinematics/16C_pd_15C_GS_171_874MeV.txt",
       kMagenta + 2},
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

   // Estructura para guardar la verdad simulada
   struct TruthData {
      double K_real;
      double theta_real; // en grados
   };

   std::map<Long64_t, TruthData> truthMap;
   // ==========================================
   // 1. LEER LA SIMULACIÓN Y LLENAR EL MAPA
   // ==========================================
   // std::cout << "Leyendo archivo de simulacion (.root)..." << std::endl;
   TFile *file_sim = TFile::Open("output_16Cpd_GS_184MeV.root");

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
   // std::cout << "Se cargaron " << truthMap.size() << " deuterones desde la simulacion." << std::endl;
   //....................................
   //  Histogramas de Residuos (Correlaciones)
   auto *h_dK_vs_theta = new TH2F(
      "h_dK_vs_theta", "#Delta K vs #theta_{rec};#theta_{rec} (deg);K_{rec} - K_{real} (MeV)", 100, 5, 40, 100, -2, 2);
   auto *h_dK_vs_K =
      new TH2F("h_dK_vs_K", "#Delta K vs K_{rec};K_{rec} (MeV);K_{rec} - K_{real} (MeV)", 100, 5, 25, 100, -2, 2);
   auto *h_dTheta_vs_K =
      new TH2F("h_dTheta_vs_K", "#Delta#theta vs K_{rec};K_{rec} (MeV);#theta_{rec} - #theta_{real} (deg)", 100, 5, 25,
               100, -2, 2);
   auto *h_dTheta_vs_theta = new TH2F(
      "h_dTheta_vs_theta", "#Delta#theta vs #theta_{rec};#theta_{rec} (deg);#theta_{rec} - #theta_{real} (deg)", 100, 5,
      40, 100, -2, 2);

   // ── Histograma 2D: Ex vs Z ───────────────────────────────────────────
   auto *ExvsZpos = new TH2F("ExvsZpos", "", 100, -5, 5, 200, -5, 80);
   auto *hex = new TH1F("hex", "C16(p,d)", NumberBins, -4.0, 4.0);
   TH2F *Ang_Ener = new TH2F("Ang_Ener", "Ang_Ener", 400, 10, 40, 1000, 0, 40.0);

   auto *ExvsZpos_cuts = new TH2F("ExvsZpos_cuts", "", 100, -5, 5, 200, -5, 100);
   auto *hex_thetaCM = new TH2F("hex_thetaCM", "C16(p,d)", 100, 0, 80, 200, -1, 1);
   auto *hex_thetaLab = new TH2F("hex_thetaLab", "C16(p,d)", 100, 0, 60, 200, -1, 1);
   auto *hex_cuts = new TH1F("hex_cuts", "C16(p,d)", 200, -0.5, 0.5);
   TH2F *Ang_Ener_cuts = new TH2F("Ang_Ener_cuts", "Ang_Ener_cuts", 400, 10, 40, 1000, 0, 40.0);

   // ── ELoss CATIMA (igual que en la macro principal) ────────────────────────
   double densityH2 = 3.553e-5; // g/cm³
   AtTools::AtELossCATIMA elossH2(densityH2);
   elossH2.SetMaterial(catima::Material(1, 1));
   elossH2.SetProjectile(16, 6, 16.0147);

   // ── Datos reconstruidos ──────────────────────────────────────────────────
   std::vector<TString> filenames;
   TChain *chain = new TChain("parquettree"); // parquettree
   for (int i = 0; i <= 19; i++) {
      char name[64];
      std::snprintf(name, sizeof(name), "run_%04d_2H.root", i);
      filenames.push_back(name);
      chain->Add(
         ("/home/georgina/my_sim/engine_Ex_GS_C16_pd/InterpSolver/interpSolverRoot/" + std::string(name)).c_str());
      // chain->Add(("/home/georgina/my_sim/engine_Ex_GS_196MeV/InterpSolver/InterpSolverRoot/" +
      // std::string(name)).c_str());
      // chain->Add(("/home/georgina/my_sim/engine_Ex_GS_186_39MeV/InterpSolver/InterpSolverRoot/" +
      // std::string(name)).c_str()); //tree kinematics: error en convertir
   }
   // std::cout << "Total entries (reconstruction): " << chain->GetEntries() << std::endl;

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

      auto [ex, theta_cm] = kine_2b(m_C16, m_p, m_b, m_B, Ebeam_at_z, theta, E_ej);

      // --- FILTRO DE PROTECCIÓN (El reemplazo del cout) ---
      // Si el cálculo cinemático da un error matemático, ignoramos este evento
      if (TMath::IsNaN(ex) || TMath::IsNaN(E_ej) || TMath::IsNaN(Ebeam_at_z)) {
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
         double theta_rec_deg = theta * TMath::RadToDeg();

         // 3. Calculamos los residuos (Reconstruido - Real)
         double dK = K_rec - K_real;
         double dTheta = theta_rec_deg - theta_real_deg;

         // 4. Llenamos los histogramas
         // (Aplicamos el mismo corte Z y Ex que usas para quedarte con el GS limpio)
         if (zPos * 100 > 2.0 && zPos * 100 < 60.0 && E_ej < 20.0) {
            h_dK_vs_theta->Fill(theta_rec_deg, dK);
            h_dK_vs_K->Fill(K_rec, dK);
            h_dTheta_vs_K->Fill(K_rec, dTheta);
            h_dTheta_vs_theta->Fill(theta_rec_deg, dTheta);
         }
      }
      // ------------------------------------

      hex->Fill(ex);
      ExvsZpos->Fill(ex, zPos * 100);

      Ang_Ener->Fill(theta * TMath::RadToDeg(), E_ej); // calibrado

      if (zPos * 100 > 2.0 && zPos * 100 < 60.0 &&
          E_ej < 20.0) { // cm y MeV, y un corte en Ex para quedarnos solo con la GS
         hex_cuts->Fill(ex);
         ExvsZpos_cuts->Fill(ex, zPos * 100);
         Ang_Ener_cuts->Fill(theta * TMath::RadToDeg(), E_ej);
         hex_thetaCM->Fill(theta_cm, ex);
         hex_thetaLab->Fill(theta * TMath::RadToDeg(), ex);
      }
   }

   // ── Canvas y dibujo ───────────────────────────────────────────────────────
   TCanvas *c = new TCanvas("c_ZvsEx_GS", "Z vs Ex: Ground State", 900, 700);
   c->Divide(2, 1);
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
   c->SaveAs("ExvsZ_GSengine.png");

   TCanvas *c_cuts = new TCanvas("c_ZvsEx_GS_cuts", "Z vs Ex with cuts: Ground State", 900, 700);
   c_cuts->Divide(2, 1);
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
   c_cuts->Update();
   c_cuts->SaveAs("ExvsZ_GS_cutsengine.png");

   TCanvas *c_hex_fit = new TCanvas("c_hex_fit", "Fit Excitation Energy", 900, 700);
   c_hex_fit->cd();
   hex_cuts->GetXaxis()->SetTitle("E_{x} (MeV)");
   hex_cuts->GetYaxis()->SetTitle("Counts");
   hex_cuts->SetTitle("Excitation Energy with cuts");
   hex_cuts->Draw();

   TF1 *f_gaus = new TF1("f_gaus", "gaus", -2.0, 2.0);
   f_gaus->SetParameters(hex_cuts->GetMaximum(), 0, 0.3); // Estima inicial: altura, media, sigma
   hex_cuts->Fit(f_gaus, "R");

   c_hex_fit->Update();
   c_hex_fit->SaveAs("hex_cuts_fit_GSengine.png");

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
   c_ang_ener_cuts->SaveAs("Theta_vs_KE_GS_cutsengine.png");

   // ── Canvas para los gráficos de residuos ─────────────────────────────────
   TCanvas *c_correlaciones = new TCanvas("c_corr", "Analisis de Residuos con cuts", 1200, 800);
   c_correlaciones->Divide(2, 2);

   c_correlaciones->cd(1);
   h_dK_vs_theta->Draw("colz");

   c_correlaciones->cd(2);
   h_dK_vs_K->Draw("colz");

   c_correlaciones->cd(3);
   h_dTheta_vs_K->Draw("colz");

   c_correlaciones->cd(4);
   h_dTheta_vs_theta->Draw("colz");

   c_correlaciones->Update();
   c_correlaciones->SaveAs("correlaciones_residuos_GSengine.png");

   TCanvas *angulos = new TCanvas("angulos", "Angulos", 900, 700);
   angulos->Divide(2, 1);
   angulos->cd(1);
   hex_thetaCM->GetXaxis()->SetTitle("#theta_{CM} (deg)");
   hex_thetaCM->GetYaxis()->SetTitle("E_{x} (MeV)");
   hex_thetaCM->SetTitle("Theta CM vs Ex");
   hex_thetaCM->Draw("colz");

   angulos->cd(2);
   hex_thetaLab->GetXaxis()->SetTitle("#theta_{lab} (deg)");
   hex_thetaLab->GetYaxis()->SetTitle("E_{x} (MeV)");
   hex_thetaLab->SetTitle("Theta Lab vs Ex");
   hex_thetaLab->Draw("colz");
   angulos->Update();
   angulos->SaveAs("angulos_GSengine.png");
}