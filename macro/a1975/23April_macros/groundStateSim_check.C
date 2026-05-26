#include <fstream>
#include <iostream>
#include <set>


double Ebin_max = 9. ;
double Ebin_min = -1. ;
int NumberBins = 100; //140

// ── Función auxiliar de cinemática (igual que en la macro principal) ──────────
Double_t omega(Double_t x, Double_t y, Double_t z)
{
   return sqrt(x*x + y*y + z*z - 2*x*y - 2*y*z - 2*x*z);
}

std::tuple<double,double>
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

   m4_ex = sqrt((cos(thetalab) * omega(s,pow(m1,2),pow(m2,2)) * omega(u,pow(m2,2),pow(m3,2)) -
                 (s-pow(m1,2)-pow(m2,2)) * (pow(m2,2)+pow(m3,2)-u)) /
                (2*pow(m2,2)) + s + u - pow(m2,2));
   Ex = m4_ex - m4;

   t = pow(m2,2) + pow(m4_ex,2) - 2*m2*Et4;

   theta_cm = TMath::Pi() - acos((pow(s,2) + s*(2*t - pow(m1,2) - pow(m2,2) - pow(m3,2) - pow(m4_ex,2)) +
                                  (pow(m1,2)-pow(m2,2))*(pow(m3,2)-pow(m4_ex,2))) /
                                 (omega(s,pow(m1,2),pow(m2,2)) * omega(s,pow(m3,2),pow(m4_ex,2))));
   theta_cm *= TMath::RadToDeg();

   return std::make_tuple(Ex, theta_cm);
}

// ── Función principal ─────────────────────────────────────────────────────────
void groundStateSim_check()
{
   // ── Masas (MeV/c²) — mismos valores que en la macro principal ────────────
   Double_t m_p   = 938.272076;
   Double_t m_d   = 1875.612931;
   Double_t m_C15 = 13979.218707;
   Double_t m_C16 = 14914.533798;

   double Q = m_C16 + m_p - m_d - m_C15;
   cout << "Q-value = " << Q << " MeV " << endl;
   
   // Beam and target parameters.
   Double_t Ebeam_buff = 12.25 *16; //11.5 * 16; // MeV, energía del haz en el buffer gas
   //Double_t Ebeam_buff = 11.5 * 16; // MeV, energía del haz en el buffer gas
   Double_t m_b = m_d;
   Double_t m_B = m_C15;

   // Ejectile parameters:deuterium
   int A_ej = 2;
   int Z_ej = 1;
   Double_t m_ej = m_d;

   // ── Curva teórica del Ground State (Un solo archivo) ─────────────────────
   /*TGraph *gGS = nullptr; // Aquí guardaremos nuestra gráfica

   TString fileKine = "/home/georgina/fair_install/ATTPCROOTv2/macro/Kinematics/Decay_kinematics/16C_pd_15C_GS_11_5_26May.txt";
   std::ifstream kineStr(fileKine.Data());

   if (!kineStr.is_open()) {
      std::cout << " [ERROR] No se pudo abrir el archivo de cinemática!" << std::endl;
   } else {
      std::vector<Double_t> ThetaLabRec, EnerLabRec;
      Double_t tCMS, tLabRec, eLabRec, tLabSca, eLabSca; 
      
      // Saltar la primera línea si tiene letras (nombres de columnas)
      std::string header;
      std::getline(kineStr, header);

      // Leer los datos
      while (kineStr >> tCMS >> tLabRec >> eLabRec >> tLabSca >> eLabSca) {
         ThetaLabRec.push_back(tLabRec);
         EnerLabRec.push_back(eLabRec);
      }

      // Si leímos datos correctamente, creamos la gráfica
      if(!ThetaLabRec.empty()) {
         gGS = new TGraph(ThetaLabRec.size(), ThetaLabRec.data(), EnerLabRec.data());
         gGS->SetLineColor(2); // 2 = kRed
         gGS->SetLineWidth(2);
         gGS->SetTitle("Ground State");
         std::cout << " -> Curva teórica cargada con éxito (" << ThetaLabRec.size() << " puntos)." << std::endl;
      } else {
         std::cout << " [WARNING] El archivo de cinemática está vacío o tiene un formato incorrecto." << std::endl;
      }
   }*/


   // ── Curvas teóricas: GS para 3 energías ──────────────────────────────────
struct KineFile {
   TString label;
   TString path;
   int     color;
};

std::vector<KineFile> kineFiles = {
   {"11.0 AMeV", "/home/georgina/fair_install/ATTPCROOTv2/macro/Kinematics/Decay_kinematics/16C_pd_15C_GS_11_26May.txt",   kBlue},
   {"11.5 AMeV", "/home/georgina/fair_install/ATTPCROOTv2/macro/Kinematics/Decay_kinematics/16C_pd_15C_GS_11_5_26May.txt", kRed},
   {"12.0 AMeV", "/home/georgina/fair_install/ATTPCROOTv2/macro/Kinematics/Decay_kinematics/16C_pd_15C_GS_12_26May.txt", kViolet+2},
   {"11.65 AMeV", "/home/georgina/fair_install/ATTPCROOTv2/macro/Kinematics/Decay_kinematics/16C_pd_15C_GS_11_65_26May.txt", kGreen+2},
};

std::vector<TGraph*> gKine;

for (auto& kf : kineFiles) {
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

   TGraph* g = new TGraph(theta.size(), theta.data(), ener.data());
   g->SetLineColor(kf.color);
   g->SetLineWidth(4);
   g->SetTitle(kf.label);
   gKine.push_back(g);
   std::cout << " -> " << kf.label << " cargada (" << theta.size() << " puntos)." << std::endl;
}

   // ── Histograma 2D: Ex vs Z ───────────────────────────────────────────
   auto *ExvsZpos = new TH2F("ExvsZpos", "", 100, -6, 6, 200, -5, 100);
   auto *hex = new TH1F("hex", "C16(p,d)", NumberBins, -7.0, 7.0);
   TH2F *Ang_Ener = new TH2F("Ang_Ener", "Ang_Ener", 400, 10, 40, 1000, 0, 40.0); 
   
   auto *ExvsZpos_cuts = new TH2F("ExvsZpos_cuts", "", 100, -6, 6, 200, -5, 100);
   auto *hex_cuts = new TH1F("hex_cuts", "C16(p,d)", 200, -7.0, 7.0);
   TH2F *Ang_Ener_cuts = new TH2F("Ang_Ener_cuts", "Ang_Ener_cuts", 400, 10, 40, 1000, 0, 40.0);

   // ── ELoss CATIMA (igual que en la macro principal) ────────────────────────
   double densityH2 = 3.553e-5; // g/cm³
   AtTools::AtELossCATIMA elossH2(densityH2);
   elossH2.SetMaterial(catima::Material(1, 1));
   elossH2.SetProjectile(16, 6, 16.0147);

   // ── Datos reconstruidos ──────────────────────────────────────────────────
   std::vector<TString> filenames;
   TChain *chain = new TChain("parquettree"); //parquettree
   for (int i = 0; i <= 19; i++) {
      char name[64];
      std::snprintf(name, sizeof(name), "run_%04d_2H.root", i);
      filenames.push_back(name);
      //chain->Add(("/home/georgina/my_sim/engine_Ex_GS_C16_pd/InterpSolver/interpSolverRoot/" + std::string(name)).c_str());
      chain->Add(("/home/georgina/my_sim/engine_Ex_GS_196MeV/InterpSolver/InterpSolverRoot/" + std::string(name)).c_str());

      //chain->Add(("/home/georgina/my_sim/engine_Ex_GS_186_39MeV/InterpSolver/InterpSolverRoot/" + std::string(name)).c_str()); //tree kinematics: error en convertir
   }
   std::cout << "Total entries (reconstruction): " << chain->GetEntries() << std::endl;

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

      double theta_lab_corr = theta; 
      auto [ex_corr, theta_cm_corr] = kine_2b(m_C16, m_p, m_b, m_B, Ebeam_at_z, theta_lab_corr, E_ej);

         // --- FILTRO DE PROTECCIÓN (El reemplazo del cout) ---
         // Si el cálculo cinemático da un error matemático, ignoramos este evento
         if (TMath::IsNaN(ex_corr) || TMath::IsNaN(E_ej) || TMath::IsNaN(Ebeam_at_z)) {
            continue; 
         }
      
        hex->Fill(ex_corr);
        ExvsZpos->Fill(ex_corr, zPos*100);
        Ang_Ener->Fill(theta_lab_corr* TMath::RadToDeg(), E_ej);
         
        if (zPos*100 > 2.0 && zPos*100 < 60.0 && E_ej < 20.0) {
            hex_cuts->Fill(ex_corr);
            ExvsZpos_cuts->Fill(ex_corr, zPos*100);
            Ang_Ener_cuts->Fill(theta_lab_corr* TMath::RadToDeg(), E_ej);
        }
   }

   // ── Canvas y dibujo ───────────────────────────────────────────────────────
   TCanvas *c = new TCanvas("c_ZvsEx_GS", "Z vs Ex: Ground State", 900, 700);
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
   c_cuts->Update();

   TCanvas *c_ang_ener = new TCanvas("c_ang_ener", "Angle vs Energy elation", 900, 700);
   c_ang_ener->cd();
   Ang_Ener->GetXaxis()->SetTitle("Theta_{lab} (deg)");
   Ang_Ener->GetYaxis()->SetTitle("Kinetic energy} (MeV)");
   Ang_Ener->SetTitle("Theta vs KE");
   Ang_Ener->Draw("colz");
  

   // ── Dibujar sobre el canvas existente ────────────────────────────────────
   TLegend* leg = new TLegend(0.15, 0.65, 0.40, 0.85);
   leg->SetBorderSize(0);
   leg->SetFillStyle(0);

   for (auto* g : gKine) {
      if (g) {
         g->Draw("L SAME");
         leg->AddEntry(g, g->GetTitle(), "l");
      }
   }
   leg->Draw();
    c_ang_ener->Update();

   TCanvas *c_ang_ener_cuts = new TCanvas("c_ang_ener_cuts", "Angle vs Energy with cuts", 900, 700);
   c_ang_ener_cuts->cd();
   Ang_Ener_cuts->GetXaxis()->SetTitle("Theta_{lab} (deg)");
   Ang_Ener_cuts->GetYaxis()->SetTitle("Kinetic energy} (MeV)");
   Ang_Ener_cuts->SetTitle("Theta vs KE");
   Ang_Ener_cuts->Draw("colz");
   // 4. Actualizas el canvas
   c_ang_ener->Update();
}