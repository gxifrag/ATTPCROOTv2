#include <fstream>
#include <iostream>

double Ebin_max = 9.;
double Ebin_min = -1.;
int NumberBins = 100; // 140
int NumberBinsAux = 200;

class SpectralModel {
public:
   TGraph *graphPS; // Phase-space TGraph

   SpectralModel(TGraph *g) : graphPS(g) {}

   double operator()(double *x, double *p)
   {
      double val = 0;
      val += p[0] * TMath::Gaus(x[0], p[1], p[2], false);
      val += p[3] * TMath::Gaus(x[0], p[4], p[5], false);
      val += p[6] * TMath::BreitWigner(x[0], p[7], p[8]);
      val += p[9] * TMath::BreitWigner(x[0], p[10], p[11]);
      val += p[12] * TMath::BreitWigner(x[0], p[13], p[14]);

      // Phase space from TGraph
      double ps_val = graphPS->Eval(x[0]);
      val += p[15] * ps_val;

      return val;
   }
};

TGraph *histoToTgraph(TH1F *h)
{

   auto g = new TGraph();
   for (int i = 1; i <= h->GetNbinsX(); i++) {
      g->SetPoint(i - 1, h->GetBinCenter(i), h->GetBinContent(i));
   }

   g->SetName("graphPS_1n");
   return g;
}

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

double BW_func(double *x, double *p)
{
   return p[0] * TMath::BreitWigner(x[0], p[1], p[2]);
}

TF1 *BW_pefit(const char *name, double xmin, double xmax)
{
   TF1 *fBW = new TF1(name, BW_func, xmin, xmax, 3);

   fBW->SetParNames("A", "mean", "gamma");
   fBW->SetNpx(5000);
   return fBW;
}

//---------------------------main function---------------------------------------
void C16_pd_tiltingCorrection()
{
   bool guardar_en_pdf = false; // ← cambia a false si quieres solo verlos en pantalla
   // gROOT->ProcessLine(".X /home/georgina/fair_install/ATTPCROOTv2/macro/a1975/myStyle.C");

   // Activar modo batch si estás guardando en PDF
   if (guardar_en_pdf) {
      gROOT->SetBatch(kTRUE); // ← esto evita que se abran ventanas
   } else {
      gROOT->SetBatch(kFALSE); // ← esto permite ver los canvas en pantalla
   }

   TH2F *Ang_Ener_Corr = new TH2F("Ang_Ener_Corr", "Ang_Ener_Corr", 120, 10, 40, 120, 0, 60.0);
   TH2F *Ang_Ener_tiltCorr = new TH2F("Ang_Ener_tiltCorr", "Ang_Ener_tiltCorr", 120, 10, 40, 120, 0, 60.0);

   TH2F *ELossvsBrho = new TH2F("ELossvsBrho", "ELossvsBrho", 4000, 0, 25000, 1000, 0, 4);
   TH2F *dedxvsBrho = new TH2F("dedxvsBrho", "dedxvsBrho", 4000, 0, 10000, 1000, 0, 4);
   TH2F *hVxVy = new TH2F("hVxVy", "hVxVy", 1000, 0, 4, 1000, 0, 4);
   auto *QvsEb = new TH2F("QvsEb", "QvsEb", 1000, -5, 15, 100, 0, 300);
   auto *QvsZpos = new TH2F("QvsZpos", "QvsZpos", 1000, -10, 50, 200, -100, 100);

   auto *hex_noEff = new TH1F("hex_noEff", "Sin corregir por eff", NumberBins, Ebin_min, Ebin_max);
   auto *hex = new TH1F("hex", "C16(p,d)", NumberBins, Ebin_min, Ebin_max);
   auto *hexCorr = new TH1F("hexCorr", "", NumberBins, Ebin_min, Ebin_max);
   auto *hexCorr1 = new TH1F("hexCorr1", "C16(p,d)", NumberBins, Ebin_min, Ebin_max);
   auto *hexCorr2 = new TH1F("hexCorr2", "C16(p,d)", NumberBins, Ebin_min, Ebin_max);

   auto *AngDistr = new TH1F("Ang_Distr", "Ang_Distr", 128, 0, 120);
   auto *AngDistrCM = new TH1F("Ang_Distr_CM", "Ang_Distr_CM", NumberBins, 0, 180);
   auto *ExvsZpos = new TH2F("ExvsZpos", "ExvsZpos", NumberBins, -4.5, Ebin_max, 100, -5, 105);
   auto *ExCorrvsZpos = new TH2F("ExCorrvsZpos", "ExCorrvsZpos", NumberBins, Ebin_min, Ebin_max, 100, -5, 65);
   auto *KineticEnergy = new TH1F("KineticEnergy", "KineticEnergy", 100, 0, 75);

   TH1F *h_PS_1n = new TH1F("h_PS_1n", "h_PS_1n", NumberBins, Ebin_min, Ebin_max);
   // auto *ExvsTrackLength = new TH2F("ExvsTrackLength", "ExvsTrackLength", 1000, -5, 15, 200, -20, 150);
   /*auto *hredchi2 = new TH1F("redchi2", "redchi2", 1000, 0, 0.0001);
   auto *hbredchi2 = new TH1F("bredchi2", "bredchi2", 1000, 0, 5);*/

   auto *hexvstheta_CM = new TH2F("hexVStheta_CM", "hexVStheta_CM", NumberBins, Ebin_min, Ebin_max, 180, 0, 180);
   auto *hexvstheta_lab = new TH2F("hexVStheta_lab", "hexVStheta_lab", NumberBins, Ebin_min, Ebin_max, 50, 0, 50);

   auto *KEvsEx = new TH2F("KEvsEx", "KEvsEx", 20, 0, 60, NumberBins, -4.5, Ebin_max);
   auto *ThetavsEx =
      new TH2F("ThetavsEx", "Ex vs Theta", 50, 0.0, 1.0, NumberBins, -4.5, 1.5); // Asumiendo theta entre 0 y 1 radian
   //--------------------------------------------------------------------------------------

   Double_t nc_tot[200];
   Double_t nc_PS_1n[200];
   Double_t x[200];
   Int_t nbins;

   // Some useful transformation constants.
   Double_t u_to_MeV = 931.49410242; // MeV/u (CODATA 2018)
   Double_t Brho_to_p = 1.602176634E-19;

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
   Double_t Ebeam_buff = 11.5 * 16; // MeV, energía del haz en el buffer gas
   Double_t m_b = m_d;
   Double_t m_B = m_C15;

   // Ejectile parameters:deuterium
   int A_ej = 2;
   int Z_ej = 1;
   Double_t m_ej = m_d;

   hexCorr->Sumw2();
   hexCorr1->Sumw2();
   hexCorr2->Sumw2();

   // -----------------------------KINEMATICS FOR DIFFERENT EXCITATION ENERGIES

   std::vector<std::string> files = {"C16_pd_C15_gs_Ebeam11_5.txt", "C16_pd_C15_740keV_Ebeam11_5.txt",
                                     "C16_pd_C15_3103keV_Ebeam11_5.txt", "C16_pd_C15_4780keV_Ebeam11_5.txt",
                                     "C16_pd_C15_6841keV_Ebeam11_5.txt"};

   std::vector<std::string> labels = {"Ground State", "1st Excited State", "2nd Excited State", "3rd Excited State",
                                      "4th Excited State"};

   std::vector<int> colors = {kOrange + 7, kBlue, kGreen + 2, kMagenta, kRed + 2};
   std::vector<TGraph *> graphs;

   // Loop de carga — solo carga, no dibuja
   for (size_t i = 0; i < files.size(); i++) {
      TString fileKine =
         Form("/home/georgina/fair_install/ATTPCROOTv2/macro/Kinematics/Decay_kinematics/%s", files[i].c_str());
      std::ifstream kineStr(fileKine.Data());

      if (kineStr.fail()) {
         std::cout << " Warning : No Kinematics file found for " << labels[i] << "!" << std::endl;
         continue;
      }

      std::vector<Double_t> ThetaCMS, ThetaLabRec, EnerLabRec, ThetaLabSca, EnerLabSca;
      Double_t tCMS, tLabRec, eLabRec, tLabSca, eLabSca; // recoil:deuteron, scattered:carbon-15
      while (kineStr >> tCMS >> tLabRec >> eLabRec >> tLabSca >> eLabSca) {
         ThetaCMS.push_back(tCMS);
         ThetaLabRec.push_back(tLabRec);
         EnerLabRec.push_back(eLabRec);
         ThetaLabSca.push_back(tLabSca);
         EnerLabSca.push_back(eLabSca);
      }

      if (ThetaLabRec.empty()) {
         std::cout << " Warning : Empty kinematics file for " << labels[i] << "!" << std::endl;
         continue;
      }

      TGraph *g = new TGraph(ThetaLabRec.size(), ThetaLabRec.data(), EnerLabRec.data());
      g->SetLineColor(colors[i]);
      g->SetLineWidth(2);
      g->SetTitle(labels[i].c_str());
      graphs.push_back(g);
   }

   std::vector<TString> filenames;
   double densityH2 = 3.553e-5; // g/cm³
   AtTools::AtELossCATIMA elossH2(densityH2);
   double mass{16.0147};                        // Mass of C16 in u
   elossH2.SetMaterial(catima::Material(1, 1)); // Set material to H2
   elossH2.SetProjectile(16, 6, mass);          // Set projectile to proton

   std::set<int> excluded = {111, 121, 148, 149};

   for (int i = 104; i <= 189; i++) {
      if (excluded.count(i))
         continue;
      char name[64];
      std::snprintf(name, sizeof(name), "run_%04d_2H.root", i);
      filenames.push_back(name);
   }

   cout << "Ebeam at the end of the TPC = " << elossH2.GetEnergy(Ebeam_buff, 100.0 * 10.0) << " MeV" << endl;

   for (auto filename : filenames) {
      TString fullPath =
         "/home/georgina/C16_analysis/C16_H2/output_a1975_22April_tb510_MMG20/InterpSolver/InterpSolver_pd_root/" +
         filename;
      TFile *runFile = new TFile(fullPath, "READ");

      // 1. SAFETY CHECK: Did the file actually open?
      if (!runFile || runFile->IsZombie()) {
         std::cout << "[ERROR] File missing or corrupted: " << filename << std::endl;
         if (runFile)
            runFile->Close();
         continue; // Skip to the next file
      }

      TTree *Tphysics = (TTree *)runFile->Get("parquettree");

      // 2. SAFETY CHECK: Does 'parquettree' exist inside this file?
      if (!Tphysics) {
         std::cout << "[ERROR] Tree 'parquettree' not found in: " << filename << std::endl;
         runFile->Close();
         continue; // Skip to the next file
      }

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
         Double_t E_ej = (TMath::Sqrt(p_ej * p_ej + m_ej * m_ej) - m_ej); // p_ej: MeV/c, m_ej: MeV/c^2, E_ej: MeV

         auto [ex_energy, theta_cm] = kine_2b(m_C16, m_p, m_b, m_B, Ebeam_buff, theta, E_ej); // ke

         double dist3D = TMath::Sqrt(vx_pos * vx_pos + vy_pos * vy_pos + zPos * zPos) * 100.0; // cm
         // Double_t Ebeam_at_z = elossH2.GetEnergy(Ebeam_buff, zPos * 100); //
         Double_t Ebeam_at_z =
            elossH2.GetEnergy(Ebeam_buff, dist3D * 10.0); // Usar la distancia 3D para la corrección de energía, MeV/mm

         double theta_lab_corr = theta; // Por ahora sin corrección de ángulo, solo para probar la implementación de la
                                        // corrección de energía en el cálculo de Ex y theta_cm
         auto [ex_energy_corr, theta_cm_corr] = kine_2b(m_C16, m_p, m_b, m_B, Ebeam_at_z, theta_lab_corr, E_ej);

         KineticEnergy->Fill(E_ej); // Usar energía calibrada para el histograma de energía cinética

         // if ( ex_energy_corr_tiltCorr > 0.0) {
         Ang_Ener_Corr->Fill(theta_lab_corr * TMath::RadToDeg(),
                             E_ej); // theta lab!! -> I still have to implement the correction of catima?
         // }
         double theta_deg = theta * TMath::RadToDeg();

         // Fill corrected histogram
         KEvsEx->Fill(E_ej, ex_energy_corr); // MeV, MeV sin corrección tilt

         ExvsZpos->Fill(ex_energy_corr, zPos * 100.0); // MeV, cm
         hexCorr2->Fill(ex_energy_corr);               // MeV, sin corrección tilt
         ThetavsEx->Fill(theta, ex_energy_corr);       // Radianes, MeV sin corrección tilt

         if (zPos * 100 > 2.0 && zPos * 100 < 60.0 && E_ej < 14.0) { // cm y MeV (zPos is in meters, E_ej is in MeV)
            hex->Fill(ex_energy_corr); // Llenar el histograma con corrección de eficiencia
         }

         // Histograms

         Double_t vx = TMath::Sin(theta) * TMath::Sqrt(ke);
         Double_t vy = TMath::Cos(theta) * TMath::Sqrt(ke);

         hVxVy->Fill(vx, vy);

         AngDistr->Fill(theta * TMath::RadToDeg());
         AngDistrCM->Fill(theta_cm);
         hexvstheta_CM->Fill(ex_energy_corr, theta_cm);
         hexvstheta_lab->Fill(ex_energy_corr, theta_lab_corr * TMath::RadToDeg());
         // ExvsTrackLength->Fill(ex_energy_corr, dist3D);

         // tEvents->Fill();
      } // events
   }    // Files

   AngDistrCM->Divide(new TF1("sin", "sin(x * TMath::DegToRad())", 0, 180));

   //------------------------------------------------------------------------------------------------
   TCanvas *c_KEvsEx = new TCanvas("c_KEvsEx", "Excitation Energy vs Kinetic Energy", 1200, 800);
   KEvsEx->GetYaxis()->SetTitle("Excitation Energy (MeV)");
   KEvsEx->GetXaxis()->SetTitle("Kinetic Energy of Ejectile (MeV)");
   KEvsEx->GetYaxis()->SetRangeUser(-4.5, 1.5);

   TProfile *perfil_GS = KEvsEx->ProfileX("perfil_GS");

   // Hacer el ajuste lineal (Paso 4)
   perfil_GS->Fit("pol1", "Q", "", 5, 45);
   TF1 *ajuste = perfil_GS->GetFunction("pol1");

   // Guardar los parámetros para tu corrección cinemática
   double pendiente = ajuste->GetParameter(1);
   double interseccion = ajuste->GetParameter(0);

   std::cout << "======================================" << std::endl;
   std::cout << "PARAMETROS DEL AJUSTE (Ground State)" << std::endl;
   std::cout << "Intersección (p0): " << interseccion << " MeV" << std::endl;
   std::cout << "Pendiente (p1)   : " << pendiente << " (Cambio de Ex por MeV de Eej)" << std::endl;
   std::cout << "======================================" << std::endl;

   // Restaurar el eje Y para verlo completo
   KEvsEx->GetYaxis()->UnZoom();

   // Formato visual
   perfil_GS->SetMarkerStyle(20);
   perfil_GS->SetMarkerColor(kRed);
   perfil_GS->SetLineColor(kRed);

   // Dibujar todo en el Canvas
   TCanvas *c1 = new TCanvas("c1", "Cinematica", 800, 600);
   KEvsEx->Draw("colz");
   perfil_GS->Draw("same E1");

   //-----------------------------------------------------------------------------------
   // 1. Sacar el perfil y ajustar
   TProfile *perfil_Theta = ThetavsEx->ProfileX("perfil_Theta");
   perfil_Theta->Fit("pol1", "Q"); // Ajuste lineal
   TF1 *ajuste_theta = perfil_Theta->GetFunction("pol1");

   // 2. Obtener la pendiente (m2)
   double m2 = ajuste_theta->GetParameter(1); // Esto es dEx / dTheta (MeV/rad)

   // 3. Calcular TU constante
   double m1 = -0.0516547;             // La pendiente que ya sacaste antes de KEvsEx
   double mi_constante_tilt = m1 / m2; // (rad / MeV)

   std::cout << "Pendiente Ex vs Theta (m2): " << m2 << " MeV/rad" << std::endl;
   std::cout << "--> MI CONSTANTE CALCULADA: " << mi_constante_tilt << " rad/MeV" << std::endl;

   TCanvas *hexthetaCanvas_CMandlab = new TCanvas("hexvsthetaCanvas_CMandlab", "hexvstheta_CM and lab", 1200, 800);
   hexthetaCanvas_CMandlab->Divide(2, 1); // dos pads: arriba el 2D, abajo la proyección

   // Pad 1: el 2D completo
   hexthetaCanvas_CMandlab->cd(1);
   hexvstheta_CM->GetXaxis()->SetTitle("Excitation Energy (MeV)");
   hexvstheta_CM->GetYaxis()->SetTitle("#theta_{CM} (#circ)");
   hexvstheta_CM->Draw("colz");

   // Pad 2: el 2D completo
   hexthetaCanvas_CMandlab->cd(2);
   hexvstheta_lab->GetXaxis()->SetTitle("Excitation Energy (MeV)");
   hexvstheta_lab->GetYaxis()->SetTitle("#theta_{lab} (#circ)");
   hexvstheta_lab->Draw("colz");

   //---------------- Combined Canvas ----------------//
   TCanvas *c_combined = new TCanvas("c_combined", "Combined Ang Energy", 1200, 800);
   c_combined->Divide(2, 1); // 2 columns, 2 rows

   // --- PAD 1: Ang_Ener_Corr (uncorrected kinematics) ---
   c_combined->cd(1);
   Ang_Ener_Corr->Draw("col");
   Ang_Ener_Corr->GetXaxis()->SetTitle("#theta_{lab} (deg)");
   Ang_Ener_Corr->GetYaxis()->SetTitle("Kinetic Energy (MeV)");

   for (size_t i = 0; i < graphs.size(); i++) {
      if (graphs[i] == nullptr || graphs[i]->GetN() == 0)
         continue;
      graphs[i]->Draw("L SAME");
   }
   auto leg0 = new TLegend(0.45, 0.7, 0.9, 0.9);
   for (size_t i = 0; i < graphs.size(); i++) {
      if (graphs[i] == nullptr)
         continue;
      leg0->AddEntry(graphs[i], labels[i].c_str(), "l");
   }
   leg0->Draw();

   // --- PAD 2: Ang_Ener_tiltCorr (corrected angle kinematics) ---
   c_combined->cd(2);
   Ang_Ener_tiltCorr->Draw("col");
   Ang_Ener_tiltCorr->GetXaxis()->SetTitle("#theta_{lab} (deg)");
   Ang_Ener_tiltCorr->GetYaxis()->SetTitle("Kinetic Energy (MeV)");

   for (size_t i = 0; i < graphs.size(); i++) {
      if (graphs[i] == nullptr || graphs[i]->GetN() == 0)
         continue;
      graphs[i]->Draw("L SAME");
   }
   auto leg2 = new TLegend(0.45, 0.7, 0.9, 0.9);
   for (size_t i = 0; i < graphs.size(); i++) {
      if (graphs[i] == nullptr)
         continue;
      leg2->AddEntry(graphs[i], labels[i].c_str(), "l");
   }
   leg2->Draw();

   c_combined->Update();

   //---------------------combined 2
   TCanvas *c_combined2 = new TCanvas("c_combined2", "Combined Ang Energy and Ex spectra", 1200, 800);
   c_combined2->Divide(2, 1); // 1 column, 3 rows

   // --- PAD 1: Ang_Ener_Corr (cut E_ej < 14) ---
   c_combined2->cd(1);
   Ang_Ener_Corr->Draw("col");
   Ang_Ener_Corr->GetXaxis()->SetTitle("#theta_{lab} (deg)");
   Ang_Ener_Corr->GetYaxis()->SetTitle("Kinetic Energy (MeV)");

   TLine *line1 = new TLine(Ang_Ener_Corr->GetXaxis()->GetXmin(), 20, Ang_Ener_Corr->GetXaxis()->GetXmax(), 20);
   line1->SetLineColor(kRed);
   line1->SetLineStyle(2);
   line1->SetLineWidth(2);
   line1->Draw("SAME");
   for (size_t i = 0; i < graphs.size(); i++) {
      if (graphs[i] == nullptr || graphs[i]->GetN() == 0)
         continue;
      graphs[i]->Draw("L SAME");
   }

   TLine *line2 = new TLine(Ang_Ener_Corr->GetXaxis()->GetXmin(), 14, Ang_Ener_Corr->GetXaxis()->GetXmax(), 14);
   line2->SetLineColor(kRed);
   line2->SetLineStyle(2);
   line2->SetLineWidth(2);
   line2->Draw("SAME");

   auto leg12 = new TLegend(0.45, 0.7, 0.9, 0.9);
   for (size_t i = 0; i < graphs.size(); i++) {
      if (graphs[i] == nullptr)
         continue;
      leg12->AddEntry(graphs[i], labels[i].c_str(), "l");
   }
   leg12->Draw();

   // --- PAD 2: hexCorr (uncorrected excitation energy) ---
   c_combined2->cd(2);
   hexCorr2->GetXaxis()->SetTitle("Excitation Energy (MeV)");
   hexCorr2->GetYaxis()->SetTitle("Counts");
   hexCorr2->SetLineColor(kGreen + 2);
   hexCorr2->SetLineWidth(2);
   hexCorr2->Draw("HIST");
   hex->SetLineColor(kBlue);
   hex->Draw("same HIST");

   auto leg42 = new TLegend(0.45, 0.8, 0.9, 0.9);
   leg42->AddEntry(hexCorr2, "hexCorr2", "l");
   leg42->AddEntry(hex, "hex", "l");
   leg42->Draw();
   c_combined2->Update();

   TCanvas *c_ExvsZpos = new TCanvas("ExvsZpos", "Excited Energy spectra: comparison", 1200, 800);
   c_ExvsZpos->Divide(2, 1);
   c_ExvsZpos->cd(1);
   ExvsZpos->Draw("zcol");
   ExvsZpos->GetYaxis()->SetTitle("Z position (mm)");
   ExvsZpos->GetXaxis()->SetTitle("Excitation Energy (MeV)");
   c_ExvsZpos->cd(2);
   ExCorrvsZpos->GetYaxis()->SetTitle("Z position (mm)");
   ExCorrvsZpos->GetXaxis()->SetTitle("Excitation Energy (MeV)");
   ExCorrvsZpos->Draw("zcol");

   TCanvas *c_AngDistr = new TCanvas("AngDistr", "Angular Distribution", 1200, 600);
   c_AngDistr->cd();
   c_AngDistr->Divide(2, 1);
   c_AngDistr->cd(1);
   AngDistr->Sumw2();
   AngDistr->Draw("E1");
   c_AngDistr->cd(2);
   AngDistrCM->Sumw2();
   AngDistrCM->Draw("E1");
   AngDistrCM->GetXaxis()->SetTitle("Angle (deg)");
   AngDistrCM->GetYaxis()->SetTitle("#frac{d#sigma}{d#Omega} (a.u.)");

   /*TCanvas *kin = new TCanvas( "kin", "kin", 1200, 800);
   KineticEnergy->Sumw2();
   kin->cd();
   KineticEnergy->Draw("E1");
   TLine *line_ke = new TLine(14, 0, 14, 720);
   line_ke->SetLineColor(kRed);
   line_ke->SetLineStyle(2);
   line_ke->SetLineWidth(2);
   line_ke->Draw("same");*/

   /*TCanvas * c_hexCompare = new TCanvas("c_hexCompare", "Comparison hex vs hexCorr2", 1200, 800);
   c_hexCompare->cd();
   hex->SetLineColor(kGreen+2);
   hex->SetLineWidth(2);
   hex->Draw("HIST");
   hexCorr2->SetLineColor(kBlue);
   hexCorr2->SetLineWidth(2);
   hexCorr2->Draw("same HIST");

   auto leg52 = new TLegend(0.75, 0.8, 0.9, 0.9);
   leg52->AddEntry(hexCorr2, "hexCorr2", "l");
   leg52->AddEntry(hex, "hex", "l");
   leg52->Draw("same");*/
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
