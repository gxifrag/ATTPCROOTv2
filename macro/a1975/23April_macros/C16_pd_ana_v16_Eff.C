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
void C16_pd_ana_v16_Eff()
{
   bool guardar_en_pdf = false; // ← cambia a false si quieres solo verlos en pantalla
   // gROOT->ProcessLine(".X /home/georgina/fair_install/ATTPCROOTv2/macro/a1975/myStyle.C");

   // Activar modo batch si estás guardando en PDF
   if (guardar_en_pdf) {
      gROOT->SetBatch(kTRUE); // ← esto evita que se abran ventanas
   } else {
      gROOT->SetBatch(kFALSE); // ← esto permite ver los canvas en pantalla
   }
   TH2F *Ang_Ener_Corr = new TH2F("Ang_Ener_Corr", "Ang_Ener_Corr", 720, 10, 40, 1000, 0, 60.0);
   TH2F *Ang_Ener_tiltCorr = new TH2F("Ang_Ener_tiltCorr", "Ang_Ener_tiltCorr_tilt", 720, 10, 40, 1000, 0, 60.0);

   TH2F *Ebeam_test = new TH2F("Ebeam_test", "Ebeam_test", 1000, -2, 10, 60, 0, 300);

   TH2F *ELossvsBrho = new TH2F("ELossvsBrho", "ELossvsBrho", 4000, 0, 25000, 1000, 0, 4);
   TH2F *dedxvsBrho = new TH2F("dedxvsBrho", "dedxvsBrho", 4000, 0, 10000, 1000, 0, 4);
   TH2F *hVxVy = new TH2F("hVxVy", "hVxVy", 1000, 0, 4, 1000, 0, 4);
   auto *QvsEb = new TH2F("QvsEb", "QvsEb", 1000, -5, 15, 100, 0, 300);
   auto *QvsZpos = new TH2F("QvsZpos", "QvsZpos", 1000, -10, 50, 200, -100, 100);

   auto *hex_noEff = new TH1F("hex_noEff", "Sin corregir por eff", NumberBins, Ebin_min, Ebin_max);
   auto *hexCorr = new TH1F("hexCorr", "", NumberBins, Ebin_min, Ebin_max);
   auto *hexCorr1 = new TH1F("hexCorr1", "C16(p,d)", NumberBins, Ebin_min, Ebin_max);
   auto *hexCorr2 = new TH1F("hexCorr2", "C16(p,d)", NumberBins, Ebin_min, Ebin_max);

   auto *AngDistr = new TH1F("Ang_Distr", "Ang_Distr", 128, 0, 120);
   auto *AngDistrCM = new TH1F("Ang_Distr_CM", "Ang_Distr_CM", NumberBins, 0, 100);
   auto *ExvsZpos = new TH2F("ExvsZpos", "ExvsZpos", 1000, -5, 10, 200, -5, 80);
   // auto *ExvsTrackLength = new TH2F("ExvsTrackLength", "ExvsTrackLength", 1000, -5, 15, 200, -20, 150);
   auto *ExCorrvsZpos = new TH2F("ExCorrvsZpos", "ExCorrvsZpos", NumberBins, Ebin_min, Ebin_max, 200, -5, 80);
   auto *KineticEnergy = new TH1F("KineticEnergy", "KineticEnergy", 100, 0, 100);

   TH1F *h_PS_1n = new TH1F("h_PS_1n", "h_PS_1n", NumberBins, Ebin_min, Ebin_max);

   /*auto *hredchi2 = new TH1F("redchi2", "redchi2", 1000, 0, 0.0001);
   auto *hbredchi2 = new TH1F("bredchi2", "bredchi2", 1000, 0, 5);*/

   auto *hexvstheta_CM = new TH2F("hexVStheta_CM", "hexVStheta_CM", NumberBins, Ebin_min, Ebin_max, 80, 0, 80);
   auto *hexvstheta_lab = new TH2F("hexVStheta_lab", "hexVStheta_lab", NumberBins, Ebin_min, Ebin_max, 50, 0, 50);

   //--------------------------------------------------------------------------------------

   Double_t nc_tot[200];
   Double_t nc_PS_1n[200];
   Double_t x[200];
   Int_t nbins;

   // Some useful transformation constants.
   Double_t u_to_MeV = 931.49410242; // MeV/u (CODATA 2018)
   Double_t Brho_to_p = 1.602176634E-19;

   // Some masses that may be useful for the experiment.

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

   // Cargar eficiencia
   TFile *fEff = TFile::Open(
      "/home/georgina/fair_install/ATTPCROOTv2/macro/a1975/23April_macros/efficiencies/efficiency.root", "READ");
   if (!fEff || fEff->IsZombie()) {
      std::cerr << "ERROR: No se pudo abrir el archivo de eficiencia!" << std::endl;
      return;
   }

   TEfficiency *gEff = (TEfficiency *)fEff->Get("pEff"); // keep only this one
   if (!gEff) {
      std::cerr << "ERROR: No se encontró pEff!" << std::endl;
      return;
   }
   cout << "Total histogram entries: " << gEff->GetTotalHistogram()->GetEntries() << endl;
   cout << "Passed histogram entries: " << gEff->GetPassedHistogram()->GetEntries() << endl;

   hexCorr->Sumw2();
   hexCorr1->Sumw2();
   hexCorr2->Sumw2();

   // -----------------------------KINEMATICS FOR DIFFERENT EXCITATION ENERGIES
   std::vector<std::string> files = {"C16_pd_C15_gs_Ebeam11_5.txt", "C16_pd_C15_740keV_Ebeam11_5.txt",
                                     "C16_pd_C15_3103keV_Ebeam11_5.txt", "C16_pd_C15_4780keV_Ebeam11_5.txt",
                                     "C16_pd_C15_6841keV_Ebeam11_5.txt"};
   /*std::vector<std::string> labels = {
      "Ground State", "1st Excited State (740keV)", "2nd Excited State (3103keV)",
      "3rd Excited State (4780keV)", "4th Excited State(6841keV)"
   };*/

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

   // ELoss tables.
   // AtTools::AtELossTable *elossTableH2 = new AtTools::AtELossTable();
   // elossTableH2->LoadSrimTable("StoppingPower_SRIM_C16_H2.txt"); //SRIM no me va.
   // elossTableH2->LoadLiseTable("StoppingPower_C16_H2.txt", 2.0158,3.3084e-5);

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

   for (auto filename : filenames) {
      TFile *runFile = new TFile(
         "/home/georgina/C16_analysis/C16_H2/output_a1975_22April_tb510_MMG20/InterpSolver/InterpSolver_pd_root/" +
            filename,
         "R");
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

         auto [ex_energy, theta_cm] = kine_2b(m_C16, m_p, m_b, m_B, Ebeam_buff, theta, E_ej); // ke

         double dist3D = TMath::Sqrt(vx_pos * vx_pos + vy_pos * vy_pos + zPos * zPos) * 100.0; // cm
         // Double_t Ebeam_at_z = elossH2.GetEnergy(Ebeam_buff, zPos * 100); //
         Double_t Ebeam_at_z =
            elossH2.GetEnergy(Ebeam_buff, dist3D * 10.0); // Usar la distancia 3D para la corrección de energía, MeV/mm

         // Corrección cinemática
         double kethe = 13.;
         double theta_lab_corr_tilt =
            theta -
            (2.0 * TMath::Pi() / 4000) *
               (E_ej -
                kethe); // corrección lineal de ángulo, con pendiente de 2pi/4000 rad/MeV, y con un offset de 13 MeV
                        // (ke-the) para que la corrección sea cero en ke=13 MeV. Este valor de ke-the lo elegí porque
                        // es el valor de ke donde veo que la corrección de energía hace que ex_energy_corr_tilt sea
                        // igual a ex_energy sin corrección, es decir, donde la corrección de energía no tiene efecto en
                        // el cálculo de ex_energy, por lo que la corrección de ángulo también debería ser cero para que
                        // ex_energy_corr_tilt sea igual a ex_energy sin corrección.

         double theta_lab_corr = theta; // Por ahora sin corrección de ángulo, solo para probar la implementación de la
                                        // corrección de energía en el cálculo de Ex y theta_cm
         auto [ex_energy_corr, theta_cm_corr] = kine_2b(m_C16, m_p, m_b, m_B, Ebeam_at_z, theta_lab_corr, E_ej);

         auto [ex_energy_corr_tilt, theta_cm_corr_tilt] = kine_2b(m_C16, m_p, m_b, m_B, Ebeam_at_z, theta_lab_corr_tilt,
                                                                  E_ej); // energies: MeV, angles: radians

         // Convertir theta_cm a grados para evaluar la eficiencia
         double theta_cm_deg = theta_cm_corr_tilt; // theta_cm_corr_tilt ya está en grados según la función kine_2b
                                                   // Protección contra valores no numéricos o fuera de rango
         if (std::isnan(theta_cm_deg) || std::isinf(theta_cm_deg)) {
            continue;
         }

         // cout << "theta_cm_deg" << theta_cm_deg << "ex_ener_corr_tilt" << ex_energy_corr_tilt << endl;
         hex_noEff->Fill(ex_energy_corr_tilt); // Llenar el histograma sin corregir por eficiencia, sin cortes, "RAW"

         KineticEnergy->Fill(E_ej); // Usar energía calibrada para el histograma de energía cinética

         Ang_Ener_Corr->Fill(theta_lab_corr * TMath::RadToDeg(),
                             E_ej); // theta lab!! -> I still have to implement the correction of catima?
         Ang_Ener_tiltCorr->Fill(theta_lab_corr_tilt * TMath::RadToDeg(), E_ej); // calibrado

         double theta_deg = theta * TMath::RadToDeg();

         // Corte en región plana de eficiencia
         if (theta_cm_deg < 20 || theta_cm_deg > 160)
            continue;

         // Instead of: double eff = gEff->Eval(theta_cm_deg);
         int bin = gEff->GetTotalHistogram()->FindFixBin(theta_cm_deg);
         double eff = gEff->GetEfficiency(bin);

         if (eff <= 0)
            continue;

         // Fill corrected histogram

         if (zPos * 100 > 2.0 && zPos * 100 < 60.0 && E_ej < 14.0) { // cm y MeV (zPos is in meters, E_ej is in MeV)
            ExCorrvsZpos->Fill(ex_energy_corr_tilt, zPos * 100.0);   // MeV, cm
            hexCorr->Fill(ex_energy_corr_tilt, 1.0 / eff);  // Llenar el histograma con corrección de eficiencia
            hexCorr2->Fill(ex_energy_corr_tilt, 1.0 / eff); // Llenar el histograma con corrección de eficiencia
            /*hexCorr2->Fill(ex_energy_corr_tilt);
            hexCorr->Fill(ex_energy_corr_tilt);*/
         }

         if (zPos * 100 > 2.0 && zPos * 100 < 60.0 && E_ej < 20.0) { // cm y MeV
            hexCorr1->Fill(ex_energy_corr_tilt); // Llenar el histograma con corrección de eficiencia
         }
         // Histograms

         Double_t vx = TMath::Sin(theta) * TMath::Sqrt(ke);
         Double_t vy = TMath::Cos(theta) * TMath::Sqrt(ke);

         hVxVy->Fill(vx, vy);

         AngDistr->Fill(theta * TMath::RadToDeg());
         AngDistrCM->Fill(theta_cm_corr_tilt);
         hexvstheta_CM->Fill(
            ex_energy_corr_tilt,
            theta_cm_corr_tilt); // theta_cm_corr_tilt is already in degrees, ex_energy_corr_tilt is in MeV
         hexvstheta_lab->Fill(
            ex_energy_corr_tilt,
            theta_lab_corr_tilt *
               TMath::RadToDeg()); // theta_lab_corr_tilt is in radians, convert to degrees for the histogram
         // ExvsTrackLength->Fill(ex_energy_corr, arclength);

         // tEvents->Fill();
      } // events
   } // Files

   AngDistrCM->Divide(new TF1("sin", "sin(x * TMath::DegToRad())", 0, 180));

   TCanvas *c_ExvsZpos = new TCanvas("ExCorrvsZpos", "Excitation Energy vs z position and track length", 1200, 800);
   c_ExvsZpos->cd();
   ExCorrvsZpos->GetXaxis()->SetTitle("Excitation Energy (MeV)");
   ExCorrvsZpos->GetYaxis()->SetTitle("z (cm)");
   ExCorrvsZpos->Draw("zcol");

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

   //-----------------------------------------------------
   TCanvas *hexthetaCanvas = new TCanvas("hexvsthetaCanvas", "hexvstheta_CM", 1200, 800);
   hexthetaCanvas->Divide(4, 3); // dos pads: arriba el 2D, abajo la proyección

   // Pad 1: el 2D completo
   hexthetaCanvas->cd(1);
   hexvstheta_CM->GetXaxis()->SetTitle("Excitation Energy (MeV)");
   hexvstheta_CM->GetYaxis()->SetTitle("#theta_{CM} (#circ)");

   hexvstheta_CM->Draw("colz");

   // Pad 2: proyección en Ex para 20° < θ_cm < 30°
   hexthetaCanvas->cd(2);
   int bin_min = hexvstheta_CM->GetYaxis()->FindBin(20);
   int bin_max = hexvstheta_CM->GetYaxis()->FindBin(30);
   TH1D *hEx_slice = hexvstheta_CM->ProjectionX("hEx_20_30", bin_min, bin_max);
   hEx_slice->GetXaxis()->SetTitle("Excitation Energy (MeV)");
   hEx_slice->GetYaxis()->SetTitle("Counts");

   hEx_slice->SetTitle("20#circ < #theta_{CM} < 30#circ");
   hEx_slice->Draw("hist");

   // Pad 3: proyección en Ex para 30° < θ_cm < 40°
   hexthetaCanvas->cd(3);
   bin_min = hexvstheta_CM->GetYaxis()->FindBin(30);
   bin_max = hexvstheta_CM->GetYaxis()->FindBin(40);
   hEx_slice = hexvstheta_CM->ProjectionX("hEx_30_40", bin_min, bin_max);
   hEx_slice->GetXaxis()->SetTitle("Excitation Energy (MeV)");
   hEx_slice->GetYaxis()->SetTitle("Counts");

   hEx_slice->SetTitle("30#circ < #theta_{CM} < 40#circ");
   hEx_slice->Draw("hist");

   // Pad 4: proyección en Ex para 40° < θ_cm < 50°
   hexthetaCanvas->cd(4);
   bin_min = hexvstheta_CM->GetYaxis()->FindBin(40);
   bin_max = hexvstheta_CM->GetYaxis()->FindBin(50);
   hEx_slice = hexvstheta_CM->ProjectionX("hEx_40_50", bin_min, bin_max);
   hEx_slice->GetXaxis()->SetTitle("Excitation Energy (MeV)");
   hEx_slice->GetYaxis()->SetTitle("Counts");

   hEx_slice->SetTitle("40#circ < #theta_{CM} < 50#circ");
   hEx_slice->Draw("hist");

   // Pad 5: proyección en Ex para 50° < θ_cm < 60°
   hexthetaCanvas->cd(5);
   bin_min = hexvstheta_CM->GetYaxis()->FindBin(50);
   bin_max = hexvstheta_CM->GetYaxis()->FindBin(60);
   hEx_slice = hexvstheta_CM->ProjectionX("hEx_50_60", bin_min, bin_max);
   hEx_slice->GetXaxis()->SetTitle("Excitation Energy (MeV)");
   hEx_slice->GetYaxis()->SetTitle("Counts");

   hEx_slice->SetTitle("50#circ < #theta_{CM} < 60#circ");
   hEx_slice->Draw("hist");

   // Pad 6: proyección en Ex para 60° < θ_cm < 70°
   hexthetaCanvas->cd(6);
   bin_min = hexvstheta_CM->GetYaxis()->FindBin(60);
   bin_max = hexvstheta_CM->GetYaxis()->FindBin(70);
   hEx_slice = hexvstheta_CM->ProjectionX("hEx_60_70", bin_min, bin_max);
   hEx_slice->GetXaxis()->SetTitle("Excitation Energy (MeV)");
   hEx_slice->GetYaxis()->SetTitle("Counts");

   hEx_slice->SetTitle("60#circ < #theta_{CM} < 70#circ");
   hEx_slice->Draw("hist");

   // Pad 7: proyección en Ex para 70° < θ_cm < 80°
   hexthetaCanvas->cd(7);
   bin_min = hexvstheta_CM->GetYaxis()->FindBin(70);
   bin_max = hexvstheta_CM->GetYaxis()->FindBin(80);
   hEx_slice = hexvstheta_CM->ProjectionX("hEx_70_80", bin_min, bin_max);
   hEx_slice->GetXaxis()->SetTitle("Excitation Energy (MeV)");
   hEx_slice->GetYaxis()->SetTitle("Counts");

   hEx_slice->SetTitle("70#circ < #theta_{CM} < 80#circ");
   hEx_slice->Draw("hist");

   // Pad 8: proyección en Ex para 80° < θ_cm < 90°
   hexthetaCanvas->cd(8);
   bin_min = hexvstheta_CM->GetYaxis()->FindBin(80);
   bin_max = hexvstheta_CM->GetYaxis()->FindBin(90);
   hEx_slice = hexvstheta_CM->ProjectionX("hEx_80_90", bin_min, bin_max);
   hEx_slice->GetXaxis()->SetTitle("Excitation Energy (MeV)");
   hEx_slice->GetYaxis()->SetTitle("Counts");

   hEx_slice->SetTitle("80#circ < #theta_{CM} < 90#circ");
   hEx_slice->Draw("hist");

   // Pad 9: proyección en Ex para 90° < θ_cm < 100°
   hexthetaCanvas->cd(9);
   bin_min = hexvstheta_CM->GetYaxis()->FindBin(90);
   bin_max = hexvstheta_CM->GetYaxis()->FindBin(100);
   hEx_slice = hexvstheta_CM->ProjectionX("hEx_90_100", bin_min, bin_max);
   hEx_slice->GetXaxis()->SetTitle("Excitation Energy (MeV)");
   hEx_slice->GetYaxis()->SetTitle("Counts");

   hEx_slice->SetTitle("90#circ < #theta_{CM} < 100#circ");
   hEx_slice->Draw("hist");

   // Pad 10: proyección en Ex para 100° < θ_cm < 110°
   hexthetaCanvas->cd(10);
   bin_min = hexvstheta_CM->GetYaxis()->FindBin(100);
   bin_max = hexvstheta_CM->GetYaxis()->FindBin(110);
   hEx_slice = hexvstheta_CM->ProjectionX("hEx_100_110", bin_min, bin_max);
   hEx_slice->GetXaxis()->SetTitle("Excitation Energy (MeV)");
   hEx_slice->GetYaxis()->SetTitle("Counts");

   hEx_slice->SetTitle("100#circ < #theta_{CM} < 110#circ");
   hEx_slice->Draw("hist");

   // Pad 11: proyección en Ex para 110° < θ_cm < 120°
   hexthetaCanvas->cd(11);
   bin_min = hexvstheta_CM->GetYaxis()->FindBin(110);
   bin_max = hexvstheta_CM->GetYaxis()->FindBin(120);
   hEx_slice = hexvstheta_CM->ProjectionX("hEx_110_120", bin_min, bin_max);
   hEx_slice->GetXaxis()->SetTitle("Excitation Energy (MeV)");
   hEx_slice->GetYaxis()->SetTitle("Counts");

   hEx_slice->SetTitle("110#circ < #theta_{CM} < 120#circ");
   hEx_slice->Draw("hist");

   // Pad 12: proyección en Ex para 120° < θ_cm < 130°
   hexthetaCanvas->cd(12);
   bin_min = hexvstheta_CM->GetYaxis()->FindBin(120);
   bin_max = hexvstheta_CM->GetYaxis()->FindBin(130);
   hEx_slice = hexvstheta_CM->ProjectionX("hEx_120_130", bin_min, bin_max);
   hEx_slice->GetXaxis()->SetTitle("Excitation Energy (MeV)");
   hEx_slice->GetYaxis()->SetTitle("Counts");

   hEx_slice->SetTitle("120#circ < #theta_{CM} < 130#circ");
   hEx_slice->Draw("hist");

   //-------------------- PHASE SPACE ----------------------------------------------

   nbins = hexCorr->GetNbinsX();
   int binmax = hexCorr->GetMaximumBin();
   double ThetaCM_min = 0;
   double ThetaCM_max = 180;

   TFile *filePS =
      new TFile("/home/georgina/fair_install/ATTPCROOTv2/macro/a1975/PhaseSpace/PhaseSpace_16C_pd_1n.root", "READ");
   TTree *treePS = (TTree *)filePS->Get("simulated_tree");

   double Weight_sim, Ex_cal, ThetaCM_cal;

   treePS->SetBranchAddress("Weight_sim", &Weight_sim);
   treePS->SetBranchAddress("Ex_cal", &Ex_cal);
   treePS->SetBranchAddress("ThetaCM_cal", &ThetaCM_cal);

   for (int i = 0; i < treePS->GetEntries(); i++) {
      treePS->GetEntry(i);

      if (ThetaCM_cal > ThetaCM_min && ThetaCM_cal < ThetaCM_max) {
         h_PS_1n->Fill(Ex_cal, Weight_sim);
      }
   }
   h_PS_1n->Smooth();

   TH1F *hPS_prefit = (TH1F *)h_PS_1n->Clone("hPS_prefit");
   TH1F *hPS_final = (TH1F *)h_PS_1n->Clone("hPS_final");

   TGraph *graphPS = histoToTgraph(h_PS_1n);

   //---------------- Fitting the experimental data ----------------//

   ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
   TCanvas *c_prefits = new TCanvas("prefits", "", 1200, 800);
   c_prefits->cd();

   TSpectrum *sp = new TSpectrum(5);              // 5 maxima search
   int nfound = sp->Search(hexCorr, 2, "", 0.02); // 2 = sigma of smoothing, last = threshold
   Double_t *xpeaks = sp->GetPositionX();

   // copy into a vector<double>
   std::vector<double> sorted_peaks(xpeaks, xpeaks + nfound);
   // sort them
   std::sort(sorted_peaks.begin(), sorted_peaks.end());
   for (int i = 0; i < nfound; ++i) {
      // cout << "Peak " << i << " at Ex = " << sorted_peaks[i] << " MeV" << endl;
   }

   // Define a two-gaussian TF1 (ROOT built-in gaus uses amplitude = height)
   TF1 *f2g = new TF1("f2g", "gaus(0) + gaus(3)", -1.0, 1.6);

   // 1. Declaramos las variables primero
   double m1 = 0.0, m2 = 0.0, m3 = 0.0, m4 = 0.0, m5 = 0.0;

   // 2. Avisamos si faltan picos, pero QUITAMOS el 'return' para que el script siga
   if (nfound < 5) {
      std::cout << "[AVISO]: Solo encontré " << nfound << " picos de forma automática. Usando valores de respaldo.\n";
   }

   // 3. Asignamos los valores de forma segura revisando cuántos picos hay de verdad
   m1 = (nfound > 0) ? sorted_peaks[0] : 0.0;
   m2 = (nfound > 1) ? sorted_peaks[1] : 0.5;
   m3 = (nfound > 2) ? sorted_peaks[2] : 3.87;
   m4 = (nfound > 3) ? sorted_peaks[3] : 5.0; // Si nfound es 4, usa sorted_peaks[3]
   m5 = (nfound > 4) ? sorted_peaks[4] : 7.3; // Si nfound es 4, usa el respaldo 7.3 de forma segura

   cout << "Identified peaks at: " << endl;
   cout << "Peak 1: " << m1 << " MeV" << endl;
   cout << "Peak 2: " << m2 << " MeV" << endl;
   cout << "Peak 3: " << m3 << " MeV" << endl;
   cout << "Peak 4: " << m4 << " MeV" << endl;
   cout << "Peak 5: " << m5 << " MeV" << endl;

   // Peak 1
   f2g->SetParameter(1, m1);                                           // mean of gaus(0)
   f2g->SetParameter(2, 0.2);                                          // sigma guess
   f2g->SetParameter(0, hexCorr->GetBinContent(hexCorr->FindBin(m1))); // amplitude guess

   // Peak 2
   f2g->SetParameter(4, m2);                                           // mean of gaus(3)
   f2g->SetParameter(5, 0.2);                                          // sigma guess
   f2g->SetParameter(3, hexCorr->GetBinContent(hexCorr->FindBin(m2))); // amplitude guess

   hexCorr->Fit(f2g, "RQ"); // R = use the range you specified, 0 no plot, M minuit
   f2g->SetLineColor(kViolet + 2);
   f2g->SetNpx(5000); // aumenta el número de puntos para un ajuste más suave

   hexCorr->Draw();
   f2g->Draw("same hist"); // <--- REQUIRED so the fit curve is drawn

   // 2nd excited — bien definido en ~3.87 MeV
   TF1 *bwprefit1 = new TF1("bw1", "[0]*TMath::BreitWigner(x,[1],[2])", 2.5, 4.0);
   bwprefit1->SetParameters(150, m3, 0.5);
   bwprefit1->SetLineColor(kRed);
   bwprefit1->SetNpx(5000); // aumenta el número de puntos para un ajuste más suave
   hexCorr->Fit(bwprefit1, "RQ");

   // 3rd excited — muy ancho, ampliar rango
   TF1 *bwprefit2 = new TF1("bw2", "[0]*TMath::BreitWigner(x,[1],[2])", 4.7, 6);
   bwprefit2->SetParameters(50, m4, 2.0);
   bwprefit2->SetParLimits(2, 0.1, 5.0); // gamma acotado
   bwprefit2->SetLineColor(kBlue);
   bwprefit2->SetNpx(5000); // aumenta el número de puntos para un ajuste más suave
   hexCorr->Fit(bwprefit2, "RQ+");

   // 4th excited — fijar posición manualmente, TSpectrum lo ve en 7.3
   TF1 *bwprefit3 = new TF1("bw3", "[0]*TMath::BreitWigner(x,[1],[2])", 6., 8.0);
   bwprefit3->SetParameters(40, 7.3, 1.0);
   // bwprefit3->SetParLimits(2, 0.1, 3.0);
   bwprefit3->SetLineColor(kGreen);
   bwprefit3->SetNpx(5000); // aumenta el número de puntos para un ajuste más suave
   hexCorr->Fit(bwprefit3, "RQ+");

   double phaseSpace = 0.0005; // valor inicial para el ajuste, se puede ajustar según la escala de los datos

   if (hPS_prefit->GetNbinsX() > 0) {
      hPS_prefit->Scale(phaseSpace); // Normaliza a la integral deseada
   }

   hPS_prefit->SetLineColor(kGray + 2);
   hPS_prefit->SetLineWidth(2);
   hPS_prefit->SetLineStyle(2);   // ← 2 = línea discontinua
   hPS_prefit->SetMarkerSize(0);  // ← elimina los puntos
   hPS_prefit->Draw("same HIST"); // ← HIST fuerza línea, sin marcadores

   double A1 = bwprefit1->GetParameter(0);
   double mean1 = bwprefit1->GetParameter(1);
   double gamma1 = bwprefit1->GetParameter(2);

   double A2 = bwprefit2->GetParameter(0);
   double mean2 = bwprefit2->GetParameter(1);
   double gamma2 = bwprefit2->GetParameter(2);

   double A3 = bwprefit3->GetParameter(0);
   double mean3 = bwprefit3->GetParameter(1);
   double gamma3 = bwprefit3->GetParameter(2);

   //--------------------------------------------------------------------------------------------

   TCanvas *c_ExEner = new TCanvas("ExEner", "Corrected Excited Energy spectra", 1200, 800);
   c_ExEner->cd();

   SpectralModel *model = new SpectralModel(graphPS);
   TF1 *fModel = new TF1("fModel", model, -0.5, 8.0, 16, // number fitted parameters
                         "SpectralModel");

   fModel->SetNpx(5000);
   fModel->SetLineColor(kBlack);
   fModel->SetLineWidth(4);

   std::vector<double> globalParamsIni = {f2g->GetParameter(0), // Amp1
                                          f2g->GetParameter(1), // Mean1
                                          f2g->GetParameter(2), // Sigma1
                                          f2g->GetParameter(3), // Amp2
                                          f2g->GetParameter(4), // Mean2
                                          f2g->GetParameter(5), // Sigma2
                                          A1,
                                          mean1,
                                          gamma1,
                                          A2,
                                          mean2,
                                          gamma2,
                                          A3,
                                          mean3,
                                          gamma3,
                                          phaseSpace};
   // Assigns all the parameters at the same time
   for (size_t i = 0; i < globalParamsIni.size(); ++i) {
      fModel->SetParameter(i, globalParamsIni[i]);
      // cout << "i= " << i << "globalParamsIni= " << globalParamsIni[i] << endl;
   }

   fModel->SetParLimits(8, 0.5, 0.8);
   // 3rd peak (rosa)
   fModel->SetParLimits(11, 0.2, 0.8);

   // 4th peak (marrón)
   fModel->SetParLimits(14, 0.2, 0.8);

   // phase space
   fModel->SetParLimits(15, 1e-4, 1e-2); // o un límite razonable

   if (!hexCorr2) {
      std::cerr << "hexCorr2 is null\n";
      return;
   }
   if (!hPS_prefit) {
      std::cerr << "hPS_prefit is null\n";
      return;
   }
   if (!graphPS) {
      std::cerr << "graphPS is null\n";
      return;
   }

   hexCorr2->GetXaxis()->SetTitle("Excitation Energy (MeV)");
   hexCorr2->GetYaxis()->SetTitle("Counts");
   hexCorr2->Draw("E");              // Draw data points with error bars first
   hexCorr2->Fit(fModel, "RS SAME"); // Fit and superimpose

   int npar = fModel->GetNpar();
   std::vector<double> globalParamsFinals(npar);
   for (int i = 0; i < npar; ++i) {
      globalParamsFinals[i] = fModel->GetParameter(i);
   }

   for (int i = 0; i < npar; ++i) {
      double v = fModel->GetParameter(i);
      if (!std::isfinite(v)) {
         std::cerr << "Param[" << i << "] invalid: " << v << "\n";
         return;
      }
   }

   TLine *vline0 = new TLine(1.218, 0, 1.218, 330); // línea vertical en x=1.218 MeV;
   vline0->SetLineColor(kRed);                      // opcional
   vline0->SetLineStyle(2);                         // opcional: línea discontinua
   vline0->SetLineWidth(3);                         // opcional
   vline0->Draw("SAME");

   // Gaussian 1
   TF1 *gaus1 = new TF1("gaus1", "gaus(0)", Ebin_min, Ebin_max);
   gaus1->SetParameters(globalParamsFinals[0], globalParamsFinals[1], globalParamsFinals[2]);
   gaus1->SetLineColor(kOrange + 7);
   gaus1->SetNpx(5000); // aumenta el número de puntos para un ajuste más suave
   gaus1->Draw("same L");

   // Gaussian 2
   TF1 *gaus2 = new TF1("gaus2", "gaus(0)", Ebin_min, Ebin_max);
   gaus2->SetParameters(globalParamsFinals[3], globalParamsFinals[4], globalParamsFinals[5]);
   gaus2->SetLineColor(kBlue);
   gaus2->SetNpx(5000); // aumenta el número de puntos para un ajuste más suave
   gaus2->Draw("same L");

   // Breit-Wigner 1
   TF1 *bw1 = new TF1("bw1", "[0]*TMath::BreitWigner(x,[1],[2])", Ebin_min, Ebin_max);
   bw1->SetParameters(globalParamsFinals[6], globalParamsFinals[7], globalParamsFinals[8]);
   bw1->SetLineColor(kGreen + 2);
   bw1->SetNpx(5000); // aumenta el número de puntos para un ajuste más suave
   bw1->Draw("same L");

   // Breit-Wigner 2
   TF1 *bw2 = new TF1("bw2", "[0]*TMath::BreitWigner(x,[1],[2])", Ebin_min, Ebin_max);
   bw2->SetParameters(globalParamsFinals[9], globalParamsFinals[10], globalParamsFinals[11]);
   bw2->SetLineColor(kMagenta);
   bw2->SetNpx(5000); // aumenta el número de puntos para un ajuste más suave
   bw2->Draw("same L");

   // Breit-Wigner 3
   TF1 *bw3 = new TF1("bw3", "[0]*TMath::BreitWigner(x,[1],[2])", Ebin_min, Ebin_max);
   bw3->SetParameters(globalParamsFinals[12], globalParamsFinals[13], globalParamsFinals[14]);
   bw3->SetLineColor(kRed + 2);
   bw3->SetNpx(5000); // aumenta el número de puntos para un ajuste más suave
   bw3->Draw("same L");

   // Phase Space

   if (hPS_final->GetNbinsX() > 0) {
      hPS_final->Scale(globalParamsFinals[15]); // globalParamsFinals[15] es el factor de escala ajustado para el fondo
                                                // de fase espacio
   }
   hPS_final->SetLineColor(kGray + 2);
   hPS_final->SetLineWidth(4);
   // hPS_final->SetLineStyle(2);   // ← 2 = línea discontinua
   hPS_final->SetMarkerSize(0);  // ← elimina los puntos
   hPS_final->Draw("same HIST"); // ← HIST fuerza línea, sin marcadores

   // Suppose you fitted with 'fitFcn' (could be gaus1, bw1, etc.)
   double chi2 = fModel->GetChisquare();
   int ndf = fModel->GetNDF();
   double chi2Ndf = chi2 / ndf;

   TLatex latex;
   latex.SetNDC();          // normalized coordinates
   latex.SetTextSize(0.03); // smaller than legend text
   latex.DrawLatex(0.17, 0.85, Form("#chi^{2}/NDF = %.2f", chi2Ndf));

   TLegend *legend2 = new TLegend(0.6, 0.6, 0.9, 0.9); // (x1, y1, x2, y2) en coordenadas del canvas
   legend2->AddEntry(hexCorr2, "Spectrum", "l");
   legend2->AddEntry(gaus1, "Ground State", "l");
   legend2->AddEntry(gaus2, "1st Excited State", "l");
   legend2->AddEntry(bw1, "2nd Excited State", "l");
   legend2->AddEntry(bw2, "3rd Excited State", "l");
   legend2->AddEntry(bw3, "4th Excited State", "l");
   legend2->AddEntry(hPS_final, "Phase Space Bkg", "l");
   legend2->Draw("same");

   c_ExEner->Update();
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

   // --- PAD 1: Ang_Ener_Corr (cut E_ej < 20) ---
   c_combined2->cd(1);
   Ang_Ener_Corr->Draw("col");
   Ang_Ener_Corr->GetXaxis()->SetTitle("#theta_{lab} (deg)");
   Ang_Ener_Corr->GetYaxis()->SetTitle("Kinetic Energy (MeV)");

   TLine *line1 = new TLine(Ang_Ener_Corr->GetXaxis()->GetXmin(), 14, Ang_Ener_Corr->GetXaxis()->GetXmax(), 14);
   line1->SetLineColor(kRed);
   line1->SetLineStyle(2);
   line1->SetLineWidth(2);
   line1->Draw("SAME");
   for (size_t i = 0; i < graphs.size(); i++) {
      if (graphs[i] == nullptr || graphs[i]->GetN() == 0)
         continue;
      graphs[i]->Draw("L SAME");
   }

   TLine *line2 = new TLine(Ang_Ener_Corr->GetXaxis()->GetXmin(), 20, Ang_Ener_Corr->GetXaxis()->GetXmax(), 20);
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
   hexCorr1->GetXaxis()->SetTitle("Excitation Energy (MeV)");
   hexCorr1->GetYaxis()->SetTitle("Counts");
   hexCorr1->SetLineColor(kGreen + 2);
   hexCorr1->SetLineWidth(2);
   hexCorr1->Draw("HIST");
   hexCorr2->SetLineColor(kBlue);
   hexCorr2->Draw("same hist ");

   auto leg42 = new TLegend(0.45, 0.8, 0.9, 0.9);
   leg42->AddEntry(hexCorr2, "hexCorr2, ke<14 MeV", "l");
   leg42->AddEntry(hexCorr1, "hexCorr1, ke<20 MeV", "l");
   leg42->Draw("same");
   c_combined2->Update();

   TCanvas *c_ExenerCorr = new TCanvas("ExenerCorr", "Excited Energy spectra corrected", 1200, 800);
   c_ExenerCorr->cd();
   c_ExenerCorr->Divide(2, 1);
   c_ExenerCorr->cd(1);
   hexCorr->Draw("hist");
   c_ExenerCorr->cd(2);
   ExCorrvsZpos->Draw("zcol");

   TCanvas *c_AngDistr = new TCanvas("AngDistr", "Angular Distribution", 1200, 600);
   c_AngDistr->cd();
   c_AngDistr->Divide(2, 1);
   c_AngDistr->cd(1);
   AngDistr->Sumw2();
   AngDistr->Draw("hist");
   c_AngDistr->cd(2);
   AngDistrCM->Sumw2();
   AngDistrCM->Draw("hist");
   AngDistrCM->GetXaxis()->SetTitle("Angle (deg)");
   AngDistrCM->GetYaxis()->SetTitle("#frac{d#sigma}{d#Omega} (a.u.)");

   TCanvas *kin = new TCanvas("kin", "kin", 1200, 800);
   KineticEnergy->Sumw2();
   kin->cd();
   KineticEnergy->Draw("hist");
   TLine *line_ke = new TLine(20, 0, 20, 720);
   line_ke->SetLineColor(kRed);
   line_ke->SetLineStyle(2);
   line_ke->SetLineWidth(2);
   line_ke->Draw("same");

   TCanvas *hexnoEff = new TCanvas("hex_noEff", "hex_noEff", 1200, 800);
   hexnoEff->cd();
   hexCorr->SetLineColor(kRed);
   hexCorr->GetXaxis()->SetTitle("Excitation Energy (MeV)");
   hexCorr->GetYaxis()->SetTitle("Counts");
   hexCorr->Draw("hist");
   /*hex_noEff->SetLineColor(kBlue);
   hex_noEff->Draw("hist same");*/

   h_PS_1n->SetLineColor(kGreen + 2);
   h_PS_1n->Scale(0.0005); // Ajusta este factor según la escala de tus datos
   h_PS_1n->Draw("same HIST");

   TLegend *leg_noEff = new TLegend(0.6, 0.7, 0.9, 0.9);
   leg_noEff->AddEntry(hexCorr, "hexCorr (eff corrected)", "l");
   leg_noEff->AddEntry(hex_noEff, "hex_noEff (not corrected)", "l");
   leg_noEff->AddEntry(h_PS_1n, "Phase Space (scaled)", "l");
   leg_noEff->Draw("same");

   //---------------- Save plots ----------------//
   std::string nombre_pdf = "plots_C16_pd_C15.pdf";

   if (guardar_en_pdf) {
      // Usamos .Data() si nombre_pdf es TString, o .c_str() si es std::string
      // Aquí asumo que nombre_pdf es un std::string por el error que muestra ROOT

      c_ExEner->Print((nombre_pdf + "(").c_str());
      c_combined->Print(nombre_pdf.c_str());
      c_combined2->Print(nombre_pdf.c_str());
      c_ExvsZpos->Print(nombre_pdf.c_str());
      kin->Print((nombre_pdf + ")").c_str());
   }

   TFile *fOut = new TFile("analisis16Cpd_final_fits.root", "RECREATE");
   hexCorr2->Write(); // Guarda el histograma
   fModel->Write();   // Guarda el fit completo con todos sus parámetros finales
   fOut->Close();

   ofstream out("resultados16Cpd_fit_global.txt");
   // Añadimos tabulaciones para que las columnas de error queden claras
   out << "Estado\tEnergia\tErr_E\tSigma/Gamma\tErr_S/G" << endl;

   out << "GS\t" << fModel->GetParameter(1) << "\t" << fModel->GetParError(1) << "\t" << fModel->GetParameter(2) << "\t"
       << fModel->GetParError(2) << endl;
   out << "1st\t" << fModel->GetParameter(4) << "\t" << fModel->GetParError(4) << "\t" << fModel->GetParameter(5)
       << "\t" << fModel->GetParError(5) << endl;
   out << "2nd\t" << fModel->GetParameter(7) << "\t" << fModel->GetParError(7) << "\t" << fModel->GetParameter(8)
       << "\t" << fModel->GetParError(8) << endl;
   out << "3rd\t" << fModel->GetParameter(10) << "\t" << fModel->GetParError(10) << "\t" << fModel->GetParameter(11)
       << "\t" << fModel->GetParError(11) << endl;
   out << "4th\t" << fModel->GetParameter(13) << "\t" << fModel->GetParError(13) << "\t" << fModel->GetParameter(14)
       << "\t" << fModel->GetParError(14) << endl;
   out << "PS_Bkg\t" << fModel->GetParameter(15) << "\t" << fModel->GetParError(15) << endl;

   out.close();
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
