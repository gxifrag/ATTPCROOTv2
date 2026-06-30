// #include "AtELossCATIMA.h"

#include "TCanvas.h"
#include "TF1.h"
#include "TF1Convolution.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TROOT.h"

#include <fstream>
#include <iostream>

#include "/home/georgina/fair_install/ATTPCROOTv2/macro/a1975/23April_macros/penetrabilities/penetrabilities_neutron_15C_L_long0.C"
#include "/home/georgina/fair_install/ATTPCROOTv2/macro/a1975/23April_macros/penetrabilities/penetrabilities_neutron_15C_L_long1.C"
#include "/home/georgina/fair_install/ATTPCROOTv2/macro/a1975/23April_macros/penetrabilities/penetrabilities_neutron_15C_L_long2.C" // o el L que corresponda

// Crear los TGraph

TGraph *gL0 = new TGraph(num_points_L_0, energies_neutron_15C_L_0, T0_neutron_15C_values);
TGraph *gL1 = new TGraph(num_points_L_1, energies_neutron_15C_L_1, T1_neutron_15C_values);
TGraph *gL2 = new TGraph(num_points_L_2, energies_neutron_15C_L_2, T2_neutron_15C_values);

double Ebin_max = 9.;
double Ebin_min = -1.;
int NumberBins = 80; // best chi2/ndf for 80 bins, 0-8 MeV
int NumberBinsAux = 200;
const double Sn = 1.2181;

// double BreitWignerPenetrability(double E, double ER, double gamma2, TSpline3 *splinePen)
double BreitWignerPenetrability(double E, double ER, double gamma2, TGraph *gPen)
{

   double En = E - Sn;
   if (En < 0)
      return 0;

   double En_max = 7.78; // o el E_end que usaste al generar la tabla
   if (En > En_max)
      En = En_max; // clamp en vez de extrapolar sin control

   double P = gPen->Eval(En); // ← interpolación lineal
   // double P = gPen->Eval(En, nullptr, "S");

   double GammaE = 2.0 * P * gamma2;
   if (GammaE <= 0)
      return 0.0;
   double denom = (E - ER) * (E - ER) + 0.25 * GammaE * GammaE;

   if (denom < 1e-12)
      return 0.0;

   return GammaE / denom;
}

double BreitWignerPenetrabilityForTF1(double *x, double *p)
{
   double xx = x[0];
   double En = xx - Sn;
   double er = p[0];
   double gamma2 = p[1];
   if (En < 0)
      return 0;

   double En_max = 7.78; // o el E_end que usaste al generar la tabla
   if (En > En_max)
      En = En_max; // clamp en vez de extrapolar sin control

   double P = gL0->Eval(En); // ← interpolación lineal
   // double P = gPen->Eval(En, nullptr, "S");

   double GammaE = 2.0 * P * gamma2;
   if (GammaE <= 0)
      return 0.0;
   double denom = (xx - er) * (xx - er) + 0.25 * GammaE * GammaE;

   if (denom < 1e-12)
      return 0.0;

   return GammaE / denom;
}

double ConvolutedBW(double x, double ER, double sigma, double gamma2, TGraph *gPen)
{
   const int N = 1000;
   const double range = 5 * sigma;
   const double step = 2 * range / N;
   double sum = 0;

   for (int i = 0; i < N; i++) {
      double t = x - range + i * step;
      double bw = BreitWignerPenetrability(t, ER, gamma2, gPen);
      double ga = TMath::Gaus(x - t, 0, sigma, true);
      sum += bw * ga;
   }

   return sum * step;
}

class SpectralModel {
public:
   TGraph *graphPS;
   /* TSpline3 *splL0;
    TSpline3 *splL1;
    TSpline3 *splL2;*/

   TGraph *gL0;
   TGraph *gL1;
   TGraph *gL2;

   // Create pointer for convolution functions
   TF1 *fBW{};
   TF1 *fGauss{};
   TF1Convolution *fConv{};

   // SpectralModel(TGraph *ps, TSpline3 *s0, TSpline3 *s1, TSpline3 *s2) : graphPS(ps), splL0(s0), splL1(s1), splL2(s2)
   SpectralModel(TGraph *ps, TGraph *s0, TGraph *s1, TGraph *s2) : graphPS(ps), gL0(s0), gL1(s1), gL2(s2)
   {
      fBW = new TF1{"fBW", BreitWignerPenetrabilityForTF1, -10, 15, 2};
      fGauss = new TF1{"fGauss", "TMath::Gaus(x, 0, [1])", -10, 15};
      fConv = new TF1Convolution(fBW, fGauss);
   }

   double operator()(double *x, double *p)
   {
      double val = 0;
      val += p[0] * TMath::Gaus(x[0], p[1], p[2], false);
      val += p[3] * TMath::Gaus(x[0], p[4], p[5], false);

      // Set paramters for convolution of 1st unbound
      double er0 = p[7];
      double gamma0 = p[9];
      double sigma0 = p[8];
      // fConv->SetParameters(er0, gamma0, sigma0);

      // *** Componente 1: BW+penetrabilidad, convolucionada ***
      double amp1 = p[6];
      double ER1 = p[7];
      double sigma1 = p[8];
      double gamma1 = p[9];
      // val += amp1 * (*fConv)(&x[0], nullptr);
      //  val += amp1 * ConvolutedBW(x[0], ER1, sigma1, gamma1, gL1); // <- ajusta L
      val += amp1 * BreitWignerPenetrability(x[0], ER1, gamma1, gL1);

      // *** Componente 2: BW+penetrabilidad, convolucionada ***
      double amp2 = p[10];
      double ER2 = p[11];
      double sigma2 = p[12];
      double gamma2 = p[13];
      val += amp2 * ConvolutedBW(x[0], ER2, sigma2, gamma2, gL1); // <- ajusta L

      // *** Componente 3: BW+penetrabilidad, convolucionada ***
      double amp3 = p[14];
      double ER3 = p[15];
      double sigma3 = p[16];
      double gamma3 = p[17];
      val += amp3 * ConvolutedBW(x[0], ER3, sigma3, gamma3, gL2); // <- ajusta L

      val += p[18] * graphPS->Eval(x[0]);

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

//---------------------------main function---------------------------------------
void C16_pd_ana_v16_4statesFit()
{
   bool guardar_en_pdf = false; // ← cambia a false si quieres solo verlos en pantalla
   gROOT->ProcessLine(".X /home/georgina/fair_install/ATTPCROOTv2/macro/a1975/myStyle.C");

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
   auto *hexCorr1 = new TH1F("hexCorr1", "", NumberBins, Ebin_min, Ebin_max);
   auto *hexCorr2 = new TH1F("hexCorr2", "", NumberBins, Ebin_min, Ebin_max);

   auto *AngDistr = new TH1F("Ang_Distr", "Ang_Distr", 128, 0, 120);
   auto *AngDistrCM = new TH1F("Ang_Distr_CM", "Ang_Distr_CM", NumberBins, 0, 180);
   auto *ExvsZpos = new TH2F("ExvsZpos", "ExvsZpos", NumberBins, Ebin_min, Ebin_max, 100, -5, 65);
   auto *ExCorrvsZpos = new TH2F("ExCorrvsZpos", "ExCorrvsZpos", NumberBins, Ebin_min, Ebin_max, 100, -5, 65);
   auto *KineticEnergy = new TH1F("KineticEnergy", "KineticEnergy", 100, 0, 75);

   TH1F *h_PS_1n = new TH1F("h_PS_1n", "h_PS_1n", NumberBins, Ebin_min, Ebin_max);
   // auto *ExvsTrackLength = new TH2F("ExvsTrackLength", "ExvsTrackLength", 1000, -5, 15, 200, -20, 150);
   /*auto *hredchi2 = new TH1F("redchi2", "redchi2", 1000, 0, 0.0001);
   auto *hbredchi2 = new TH1F("bredchi2", "bredchi2", 1000, 0, 5);*/

   auto *hexvstheta_CM = new TH2F("hexVStheta_CM", "hexVStheta_CM", NumberBins, Ebin_min, Ebin_max, 180, 0, 180);
   auto *hexvstheta_lab = new TH2F("hexVStheta_lab", "hexVStheta_lab", NumberBins, Ebin_min, Ebin_max, 50, 0, 50);

   auto *KEvsEx = new TH2F("h_KE_Ex", "All states: KE vs E_{x};KE (MeV);E_{x} (MeV)", 200, 0, 80, 200, -3, 10);

   // Para corrección de K
   TProfile *hDeltaK_vs_Krec = new TProfile("hDeltaK_vs_Krec", "Delta K vs K_rec (GS events)", 50, 0, 15, "s");
   TH2F *h2DeltaK_vs_Krec = new TH2F("h2DeltaK_vs_Krec", "", 50, 0, 15, 100, -5, 5);

   // Para comprobar tilt
   TProfile *hEx_vs_Z = new TProfile("hEx_vs_Z", "Ex vs Z vertex", 50, -500, 500, "s");
   //--------------------------------------------------------------------------------------

   Double_t nc_tot[200];
   Double_t nc_PS_1n[200];
   Double_t x[200];
   Int_t nbins;

   // Some useful transformation constants.
   Double_t u_to_MeV = 931.49410242; // MeV/u (CODATA 2018)
   Double_t Brho_to_p = 1.602176634E-19;

   // Some masses that may be useful for the experiment.

   /*Double_t m_p = 1.007825 * u_to_MeV;
   Double_t m_d = 2.0135532 * u_to_MeV;
   Double_t m_t = 3.016049281 * u_to_MeV;
   Double_t m_He3 = 3.016029 * u_to_MeV;
   Double_t m_a = 4.00260325415 * u_to_MeV;

   Double_t m_C12 = 12.00 * u_to_MeV;
   Double_t m_C13 = 13.00335484 * u_to_MeV;
   Double_t m_C14 = 14.003242 * u_to_MeV;
   Double_t m_C15 = 15.0105999 * u_to_MeV;
   Double_t m_C16 = 16.0147 * u_to_MeV;
   Double_t m_C17 = 17.0226 * u_to_MeV;*/

   Double_t m_p = 938.272076 / 1.0;     // masa nuclear del protón
   Double_t m_d = 1875.612931 / 1.0;    // masa nuclear del deuterón
   Double_t m_C15 = 13979.218707 / 1.0; // masa nuclear del 15C
   Double_t m_C16 = 14914.533798 / 1.0; // masa nuclear del 16C

   // Correct nuclear masses (MeV/c^2)
   // = atomic mass (u) * 931.494 - Z * 0.511 (electron mass)

   // Verify Q-value
   double Q = m_C16 + m_p - m_d - m_C15;
   std::cout << "Q-value = " << Q << " MeV " << endl;

   // Beam and target parameters.
   Double_t Ebeam_buff = 11.5 * 16; // 192. yassid // MeV, energía del haz en el buffer gas 11.5MeV/u * 16 u = 184 MeV
   Double_t m_b = m_d;
   Double_t m_B = m_C15;

   // Ejectile parameters:deuterium
   int A_ej = 2;
   int Z_ej = 1;
   Double_t m_ej = m_d;

   // Set Sumw2 for histograms to properly handle errors when filling with weights
   hexCorr->Sumw2();
   hexCorr1->Sumw2();
   hexCorr2->Sumw2();

   // -----------------------------KINEMATICS FOR DIFFERENT EXCITATION ENERGIES

   std::vector<std::string> files = {"C16_pd_C15_gs_Ebeam11_5.txt", "C16_pd_C15_740keV_Ebeam11_5.txt",
                                     "C16_pd_C15_3103keV_Ebeam11_5.txt", "C16_pd_C15_4780keV_Ebeam11_5.txt"};
   std::vector<std::string> labels = {"Ground State", "1st Excited State (740keV)", "2nd Excited State (3103keV)",
                                      "3rd Excited State (4780keV)"};

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

   cout << "Ebeam at the end of the TPC = " << elossH2.GetEnergy(Ebeam_buff, 100.0 * 10.0) << " MeV" << endl;

   // --- Sanity check + fix for penetrability tables ---
   auto checkAndSort = [](TGraph *g, const char *name) {
      double prevX = -1e9;
      bool unsorted = false;
      for (int i = 0; i < g->GetN(); i++) {
         double xi, yi;
         g->GetPoint(i, xi, yi);
         if (xi <= prevX) {
            cout << "[WARNING] " << name << " no ordenado en el punto " << i << ": x=" << xi << endl;
            unsorted = true;
         }
         prevX = xi;
      }
      if (unsorted) {
         g->Sort();
         cout << "[INFO] " << name << " ordenado automáticamente." << endl;
      }
   };

   checkAndSort(gL0, "gL0");
   checkAndSort(gL1, "gL1");
   checkAndSort(gL2, "gL2");

   // --- Debug: comportamiento cerca del threshold ---
   double x0, y0;
   gL1->GetPoint(0, x0, y0);
   cout << "Primer punto tabulado en gL1: En=" << (x0 - Sn) << " MeV, P=" << y0 << endl;
   cout << "gL1->Eval(En=0.00001) = " << gL1->Eval(Sn + 0.00001) << endl;
   cout << "gL1->Eval(En=0.001)   = " << gL1->Eval(Sn + 0.001) << endl;

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

         // cout << "Ebeam at position = " << dist3D << " cm: " << Ebeam_at_z << " MeV" << endl;

         // Corrección cinemática
         double kethe = 13.;
         double theta_lab_corr_tiltCorr = theta;
         /*(theta -
          (2.0 * TMath::Pi() / 4000) * (E_ej - kethe)); // theta: rad; theta_lab_corr: rad; E_ej-kethe: MeV 29.5
          */

         double theta_lab_corr = theta; // Por ahora sin corrección de ángulo, solo para probar la implementación de la
                                        // corrección de energía en el cálculo de Ex y theta_cm
         auto [ex_energy_corr, theta_cm_corr] = kine_2b(m_C16, m_p, m_b, m_B, Ebeam_at_z, theta_lab_corr, E_ej);

         auto [ex_energy_corr_tiltCorr, theta_cm_corr_tiltCorr] =
            kine_2b(m_C16, m_p, m_b, m_B, Ebeam_at_z, theta_lab_corr_tiltCorr,
                    E_ej); // energies: MeV, angles: radians

         // Convertir theta_cm a grados para evaluar la eficiencia
         double theta_cm_deg =
            theta_cm_corr_tiltCorr; // theta_cm_corr_tiltCorr ya está en grados según la función kine_2b

         KineticEnergy->Fill(E_ej); // Usar energía calibrada para el histograma de energía cinética

         double theta_deg = theta * TMath::RadToDeg();

         Ang_Ener_tiltCorr->Fill(theta_lab_corr_tiltCorr * TMath::RadToDeg(), E_ej); // calibrado

         Ang_Ener_Corr->Fill(theta_lab_corr * TMath::RadToDeg(),
                             E_ej); // theta lab!! -> I still have to implement the correction of catima?

         // Fill corrected histogram

         if (zPos * 100 > 2.0 && zPos * 100 < 60.0 && E_ej > 5.0 && E_ej < 15.0) {
            double theta_cm_deg =
               theta_cm_corr_tiltCorr; // theta_cm_corr_tiltCorr ya está en grados según la función kine_2b

            KineticEnergy->Fill(E_ej); // Usar energía calibrada para el histograma de energía cinética

            // }
            double theta_deg = theta * TMath::RadToDeg();
            KEvsEx->Fill(E_ej, ex_energy_corr); // MeV, MeV sin corrección tilt

            // cm y MeV (zPos is in meters, E_ej is in MeV)
            ExCorrvsZpos->Fill(ex_energy_corr_tiltCorr, zPos * 100.0); // MeV, cm
            ExvsZpos->Fill(ex_energy_corr, zPos * 100.0);              // MeV, cm
            hex->Fill(ex_energy_corr);               // Llenar el histograma con corrección de eficiencia
            hexCorr->Fill(ex_energy_corr_tiltCorr);  // Llenar el histograma con corrección de eficiencia
            hexCorr2->Fill(ex_energy_corr_tiltCorr); // Llenar el histograma con corrección de eficiencia
            hexCorr1->Fill(ex_energy_corr);

            // Histograms

            Double_t vx = TMath::Sin(theta) * TMath::Sqrt(ke);
            Double_t vy = TMath::Cos(theta) * TMath::Sqrt(ke);

            hVxVy->Fill(vx, vy);

            AngDistr->Fill(theta * TMath::RadToDeg());
            AngDistrCM->Fill(theta_cm_corr_tiltCorr);
            hexvstheta_CM->Fill(ex_energy_corr_tiltCorr,
                                theta_cm_corr_tiltCorr); // theta_cm_corr_tiltCorr is already in degrees,
                                                         // ex_energy_corr_tiltCorr is in MeV
            hexvstheta_lab->Fill(
               ex_energy_corr_tiltCorr,
               theta_lab_corr_tiltCorr *
                  TMath::RadToDeg()); // theta_lab_corr_tiltCorr is in radians, convert to degrees for
                                      // the histogram ExvsTrackLength->Fill(ex_energy_corr, arclength);
         }

         // tEvents->Fill();
      } // events
      runFile->Close();
      delete runFile;
   } // Files

   AngDistrCM->Divide(new TF1("sin", "sin(x * TMath::DegToRad())", 0, 180));

   /*TCanvas *c_KEvsEx = new TCanvas("c_KEvsEx", "Excitation Energy vs Kinetic Energy", 1200, 800);
   KEvsEx->GetYaxis()->SetTitle("Excitation Energy (MeV)");
   KEvsEx->GetXaxis()->SetTitle("Kinetic Energy of Ejectile (MeV)");
   KEvsEx->Draw("colz");*/
   //-------------------- PHASE SPACE ----------------------------------------------

   nbins = hexCorr->GetNbinsX();
   int binmax = hexCorr->GetMaximumBin();
   double ThetaCM_min = 0;
   double ThetaCM_max = 180;

   TFile *filePS =
      new TFile("/home/georgina/fair_install/ATTPCROOTv2/macro/a1975/PhaseSpace/PhaseSpace_16C_pd_1n.root", "READ");
   TTree *treePS = (TTree *)filePS->Get("simulated_tree");

   if (!treePS) {
      cout << "Error!";
      return;
   }

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

   // filePS->Close();
   //---------------- Fitting the experimental data ----------------//
   ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
   TCanvas *c_prefits = new TCanvas("prefits", "", 1200, 800);
   c_prefits->cd();

   TSpectrum *sp = new TSpectrum(4);              // 5 maxima search
   int nfound = sp->Search(hexCorr, 2, "", 0.02); // 2 = sigma of smoothing, last = threshold
   Double_t *xpeaks = sp->GetPositionX();

   // copy into a vector<double>
   std::vector<double> sorted_peaks(xpeaks, xpeaks + nfound);
   // sort them
   std::sort(sorted_peaks.begin(), sorted_peaks.end());

   double m1 = sorted_peaks[0];
   double m2 = sorted_peaks[1];
   double m3 = sorted_peaks[2];
   double m4 = sorted_peaks[3];

   cout << "Picos encontrados por TSpectrum: " << nfound << endl;
   cout << "Picos ordenados:" << endl;
   for (int i = 0; i < nfound; i++) {
      cout << "Pico " << i + 1 << ": " << sorted_peaks[i] << " MeV" << endl;
   }

   if (nfound < 4) {
      cout << "[WARNING] TSpectrum no encontró los 4 picos. Usando valores por defecto para los faltantes." << endl;
   }

   // Define a two-gaussian TF1 (ROOT built-in gaus uses amplitude = height)
   TF1 *f2g = new TF1("f2g", "gaus(0) + gaus(3)", -1.0, 1.6);

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

   hexCorr->Draw("hist");
   f2g->Draw("same"); // <--- REQUIRED so the fit curve is drawn

   // Extraer valores
   double Ex1 = f2g->GetParameter(1);
   double sigma1 = f2g->GetParameter(2);

   double Ex2 = f2g->GetParameter(4);
   double sigma2 = f2g->GetParameter(5);

   // Crear arrays AHORA (cuando ya existen los valores)
   double Ex_values[2] = {Ex1, Ex2};
   double sigma_values[2] = {sigma1, sigma2};

   TF1 *bwprefit1 = new TF1("bw1", "[0]*TMath::BreitWigner(x,[1],[2])", 2.5, 3.8);
   bwprefit1->SetParameters(150, m3, 0.01);
   bwprefit1->SetLineColor(kRed);
   bwprefit1->SetNpx(5000); // aumenta el número de puntos para un ajuste más suave
   hexCorr->Fit(bwprefit1, "RQ");
   bwprefit1->Draw("same");

   TF1 *bwprefit2 = new TF1("bw2", "[0]*TMath::BreitWigner(x,[1],[2])", 3.80, 4.2);
   bwprefit2->SetParameters(50, 4.2, 1.0);
   // bwprefit2->SetParLimits(2, 0.1, 1.0); // gamma acotado
   bwprefit2->SetLineColor(kBlue);
   bwprefit2->SetNpx(5000); // aumenta el número de puntos para un ajuste más suave
   hexCorr->Fit(bwprefit2, "RQ+");
   bwprefit2->Draw("same");

   TF1 *bwprefit3 = new TF1("bw3", "[0]*TMath::BreitWigner(x,[1],[2])", 4.5, 6.0);
   bwprefit3->SetParameters(50, m4, 1.5);
   // bwprefit3->SetParLimits(2, 0., 1.8); // gamma acotado
   bwprefit3->SetLineColor(kOrange);
   bwprefit3->SetNpx(5000); // aumenta el número de puntos para un ajuste más suave
   hexCorr->Fit(bwprefit3, "RQ+");
   bwprefit3->Draw("same");

   double phaseSpace = 0.0005; // valor inicial para el ajuste, se puede ajustar según la escala de los datos

   if (hPS_prefit->GetNbinsX() > 0) {
      hPS_prefit->Scale(phaseSpace); // Normaliza a la integral deseada
   }

   hPS_prefit->SetLineColor(kGray + 2);
   hPS_prefit->SetLineWidth(2);
   hPS_prefit->SetLineStyle(2);   // ← 2 = línea discontinua
   hPS_prefit->SetMarkerSize(0);  // ← elimina los puntos
   hPS_prefit->Draw("same HIST"); // ← HIST fuerza línea, sin marcadores

   hPS_prefit->SetLineColor(kGray + 2);
   hPS_prefit->SetLineWidth(2);
   hPS_prefit->SetLineStyle(2);   // ← 2 = línea discontinua
   hPS_prefit->SetMarkerSize(0);  // ← elimina los puntos
   hPS_prefit->Draw("same HIST"); // ← HIST fuerza línea, sin marcadores

   double A1 = bwprefit1->GetParameter(0);
   double mean1 = bwprefit1->GetParameter(1);
   double gamma1 = bwprefit1->GetParameter(2);

   cout << "gamma1 =  " << gamma1 << endl;

   double A2 = bwprefit2->GetParameter(0);
   double mean2 = bwprefit2->GetParameter(1);
   double gamma2 = bwprefit2->GetParameter(2);

   double A3 = bwprefit3->GetParameter(0);
   double mean3 = bwprefit3->GetParameter(1);
   double gamma3 = bwprefit3->GetParameter(2);

   // Crear gráfico
   TGraph *gSigma = new TGraph(2, Ex_values, sigma_values);
   gSigma->SetTitle("Sigma vs Ex;E_{x} (MeV);#sigma (MeV)");
   gSigma->GetXaxis()->SetLimits(-1, 8);     // fija el rango en X
   gSigma->GetYaxis()->SetRangeUser(0, 0.5); // si quieres también rango en Y (opcional)
   gSigma->SetMarkerStyle(20);
   gSigma->SetMarkerColor(kRed + 1);
   gSigma->SetLineColor(kRed + 1);

   double Ex_mean = 0.5 * (Ex1 + Ex2);
   double sigma_mean = 0.5 * (sigma1 + sigma2);

   cout << "Ex_mean = " << Ex_mean << " MeV" << endl;
   cout << "sigma_mean = " << sigma_mean << " MeV" << endl;

   TGraph *gMean = new TGraph(1);
   gMean->SetPoint(0, Ex_mean, sigma_mean);
   gMean->SetMarkerStyle(29);
   gMean->SetMarkerSize(2.0);
   gMean->SetMarkerColor(kBlue + 2);

   // Ajuste lineal
   /*TF1 *fLin = new TF1("fLin", "pol1", -1, 8);
   gSigma->Fit(fLin, "R");

   double a = fLin->GetParameter(0);
   double b = fLin->GetParameter(1);

   TLatex *eq = new TLatex();
   eq->SetNDC(); // coordenadas normalizadas
   eq->SetTextSize(0.04);
   eq->DrawLatex(0.15, 0.85, Form("#sigma = %.3f + %.3f #times E_{x}", a, b));*/

   // Dibujar
   /*TCanvas *c_sigma = new TCanvas("c_sigma", "", 900, 600);
   gSigma->Draw("AP");
   fLin->SetLineColor(kBlue);
   fLin->Draw("same");
   gMean->Draw("P SAME");
   TLine *lineMean = new TLine(-1, sigma_mean, 8, sigma_mean);
   lineMean->SetLineColor(kGreen + 2);
   lineMean->SetLineStyle(2); // línea discontinua
   lineMean->SetLineWidth(2);
   lineMean->Draw("same");*/

   //----------------------------------------------------------------------------------------------
   TCanvas *c_ExEner = new TCanvas("ExEner", "Corrected Excited Energy spectra", 1200, 800);
   c_ExEner->cd();
   SpectralModel *model = new SpectralModel(graphPS, gL0, gL1, gL2);
   TF1 *fModel = new TF1("fModel_conv", model, -1.0, 8.5, 19, "SpectralModel");

   fModel->SetNpx(5000);
   fModel->SetLineColor(kBlack);
   fModel->SetLineWidth(4);

   std::vector<double> globalParamsIni = {
      f2g->GetParameter(0),
      f2g->GetParameter(1),
      f2g->GetParameter(2), // p0-2  Gaus1
      f2g->GetParameter(3),
      f2g->GetParameter(4),
      f2g->GetParameter(5), // p3-5  Gaus2

      A1,
      mean1,
      sigma_mean,
      gamma1, // p6-9   Comp1 (conv)

      A2,
      mean2,
      sigma_mean,
      gamma2, // p10-13 Comp2 (conv)

      A3,
      mean3,
      sigma_mean,
      gamma3, // p14-17 Comp3 (conv)

      phaseSpace // p18
   };
   // Assign parameters
   for (size_t i = 0; i < globalParamsIni.size(); ++i) {
      fModel->SetParameter(i, globalParamsIni[i]);
   }

   // Sigmas instrumentales fijas (las 3 componentes comparten la misma resolución)
   fModel->FixParameter(8, sigma_mean);  // sigma Comp1
   fModel->FixParameter(12, sigma_mean); // sigma Comp2
   fModel->FixParameter(16, sigma_mean); // sigma Comp3

   // Límites de energía/anchura
   fModel->SetParLimits(7, 2.8, 3.8);   // ER Comp1
   fModel->SetParLimits(9, 0.001, 1.0); // antes 0.001 — ahora evita el colapso a anchura ~0
   fModel->SetParLimits(10, 0.001, 1e6);

   fModel->SetParLimits(11, 3.9, 4.5);   // ER Comp2
   fModel->SetParLimits(13, 0.001, 3.0); // gamma Comp2 (opcional, antes estaba comentado)

   fModel->SetParLimits(15, 4.5, 6.0);   // ER Comp3
   fModel->SetParLimits(17, 0.001, 2.0); // gamma Comp3

   //..........................................................

   hexCorr2->Fit(fModel, "R");

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

   hexCorr2->GetXaxis()->SetTitle("Excitation Energy (MeV)");
   hexCorr2->GetYaxis()->SetTitle("Counts");
   hexCorr2->Draw("E1"); // E1
   fModel->Draw("same L");

   // Gaussian 1
   TF1 *gaus1 = new TF1("gaus1", "gaus(0)", Ebin_min, Ebin_max);
   gaus1->SetParameters(globalParamsFinals[0], globalParamsFinals[1], globalParamsFinals[2]);
   gaus1->SetLineColor(kOrange + 7);
   gaus1->SetNpx(5000);    // aumenta el número de puntos para un ajuste más suave
   gaus1->SetLineStyle(2); // línea discontinua
   gaus1->Draw("same L");

   // Gaussian 2
   TF1 *gaus2 = new TF1("gaus2", "gaus(0)", Ebin_min, Ebin_max);
   gaus2->SetParameters(globalParamsFinals[3], globalParamsFinals[4], globalParamsFinals[5]);
   gaus2->SetLineColor(kBlue);
   gaus2->SetNpx(5000);    // aumenta el número de puntos para un ajuste más suave
   gaus2->SetLineStyle(2); // línea discontinua
   gaus2->Draw("same L");

   // Componente 1
   TF1 *comp1 = new TF1(
      "comp1",
      [=](double *x, double *p) {
         return p[0] * ConvolutedBW(x[0], p[1], p[2], p[3], model->gL1);
      }, // mismo L que arriba
      Ebin_min, Ebin_max, 4);
   comp1->SetParameters(globalParamsFinals[6], globalParamsFinals[7], globalParamsFinals[8], globalParamsFinals[9]);
   comp1->SetLineColor(kGreen + 2);
   comp1->SetLineStyle(2);
   comp1->SetNpx(5000);
   comp1->Draw("same L");

   // Componente 2
   TF1 *comp2 = new TF1(
      "comp2",
      [=](double *x, double *p) {
         return p[0] * ConvolutedBW(x[0], p[1], p[2], p[3], model->gL1);
      }, // mismo L que arriba
      Ebin_min, Ebin_max, 4);
   comp2->SetParameters(globalParamsFinals[10], globalParamsFinals[11], globalParamsFinals[12], globalParamsFinals[13]);
   comp2->SetLineColor(kMagenta);
   comp2->SetLineStyle(2);
   comp2->SetNpx(5000);
   comp2->Draw("same L");

   // Componente 3
   TF1 *comp3 = new TF1(
      "comp3",
      [=](double *x, double *p) {
         return p[0] * ConvolutedBW(x[0], p[1], p[2], p[3], model->gL2);
      }, // mismo L que arriba
      Ebin_min, Ebin_max, 4);
   comp3->SetParameters(globalParamsFinals[14], globalParamsFinals[15], globalParamsFinals[16], globalParamsFinals[17]);
   comp3->SetLineColor(kCyan + 2);
   comp3->SetLineStyle(2);
   comp3->SetNpx(5000);
   comp3->Draw("same L");
   // Phase Space
   if (hPS_final->GetNbinsX() > 0) {
      hPS_final->Scale(globalParamsFinals[18]); // globalParamsFinals[14] es el factor de escala ajustado para el fondo
                                                // de fase espacio
   }
   hPS_final->SetLineColor(kGray + 2);
   hPS_final->SetLineWidth(2);
   hPS_final->SetMarkerStyle(24); // open circle
   hPS_final->SetMarkerSize(1);
   hPS_final->SetMarkerColor(kGray + 2);
   hPS_final->Draw("same P"); // P = solo marcadores, sin línea

   // Suppose you fitted with 'fitFcn' (could be gaus1, bw1, etc.)
   double chi2 = fModel->GetChisquare();
   int ndf = fModel->GetNDF();
   double chi2Ndf = chi2 / ndf;

   TPaveText *pt = new TPaveText(0.6, 0.63, 0.88, 0.85, "NDC");
   pt->SetFillColor(0);  // fondo blanco
   pt->SetFillStyle(0);  // sin relleno (transparente)
   pt->SetBorderSize(0); // sin frame
   pt->SetTextAlign(11); // left-center
   pt->SetTextSize(0.032);
   pt->SetTextFont(42); // misma fuente que ROOT por defecto
   pt->AddText("^{16}C(p,d)^{15}C  E_{beam}/A = 11.5 MeV");
   pt->AddText(Form("#sigma_{det} = %.0f keV (from g.s. + 1st)", sigma_mean * 1000.0));
   pt->AddText(Form("#chi^{2}/NDF = %.2f", chi2Ndf));
   pt->Draw("same");

   gPad->Update();                 // Asegura que el pad esté actualizado para obtener los límites correctos
   double ymax = gPad->GetUymax(); // máximo del eje Y en coordenadas del pa

   TLine *vline0 = new TLine(1.218, 0, 1.218, ymax); // línea vertical en x=1.218 MeV;
   vline0->SetLineColor(kRed);                       // opcional
   vline0->SetLineWidth(3);                          // opcional
   vline0->Draw("SAME");

   std::vector<double> Ex_theory = {0.0, 0.740, 3.103, 4.202, 4.780};
   std::vector<int> line_colors = {kOrange + 7, kBlue, kGreen + 2, kMagenta, kCyan + 1};

   /* for (int i = 0; i < Ex_theory.size(); ++i) {
       TLine *vline = new TLine(Ex_theory[i], 0, Ex_theory[i], ymax);
       vline->SetLineColor(line_colors[i]);
       // vline->SetLineStyle(2); // discontinua
       vline->SetLineWidth(2);
       vline->Draw("SAME");
    }*/

   c_ExEner->Update();

   //---------------------------------------------------
   // --- Guardar parámetros del fit global en ROOT ---
   TFile *fFitParams = new TFile("fit_params_global.root", "RECREATE");

   // Guardamos la TF1 completa (contiene parámetros + errores + chi2)
   fModel->Write("fModel_global");

   // Guardamos también graphPS si existe
   if (graphPS) {
      graphPS->Write("graphPS");
   } else {
      cout << "AVISO: graphPS es NULL, no se guarda." << endl;
   }

   // TVectorD para acceso rápido desde otros macros
   int npar_save = fModel->GetNpar(); // 13
   TVectorD params(npar_save), errors(npar_save);
   for (int i = 0; i < npar_save; ++i) {
      params[i] = fModel->GetParameter(i);
      errors[i] = fModel->GetParError(i);
   }
   params.Write("fit_parameters");
   errors.Write("fit_errors");

   // Guardamos también chi2 y NDF como TNamed para referencia
   TNamed chi2_str("chi2_ndf", Form("%.4f / %d = %.4f", fModel->GetChisquare(), fModel->GetNDF(),
                                    fModel->GetChisquare() / fModel->GetNDF()));
   chi2_str.Write();

   fFitParams->Close();
   cout << "Fit params guardados en fit_params_global.root" << endl;

   double Sn_val = 1.2181;
   for (auto [amp_idx, ER_idx, sigma_idx] :
        std::vector<std::tuple<int, int, int>>{{6, 7, 8}, {10, 11, 12}, {14, 15, 16}}) {
      double ER = fModel->GetParameter(ER_idx);
      double sigma = fModel->GetParameter(sigma_idx);
      double lower = ER - 5 * sigma;
      cout << "Comp ER=" << ER << "  rango conv. inferior = " << lower << "  Sn=" << Sn_val
           << (lower < Sn_val ? "  --> CRUZA EL THRESHOLD" : "") << endl;
   }

   //---------------- Combined Canvas ----------------//
   /* TCanvas *c_combined = new TCanvas("c_combined", "Combined Ang Energy", 1200, 800);
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
    hex->GetXaxis()->SetTitle("Excitation Energy (MeV)");
    hex->GetYaxis()->SetTitle("Counts");
    hex->SetLineColor(kGreen + 2);
    hex->SetLineWidth(2);
    hex->Draw("HIST");
    hexCorr2->SetLineColor(kBlue);
    hexCorr2->Draw("same HIST");

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
   ExCorrvsZpos->Draw("zcol");*/

   /* TCanvas *c_AngDistr = new TCanvas("AngDistr", "Angular Distribution", 1200, 600);
    c_AngDistr->cd();
    c_AngDistr->Divide(2, 1);
    c_AngDistr->cd(1);
    AngDistr->Sumw2();
    AngDistr->Draw("E1");
    c_AngDistr->cd(2);
    AngDistrCM->Sumw2();
    AngDistrCM->Draw("E1");
    AngDistrCM->GetXaxis()->SetTitle("Angle (deg)");
    AngDistrCM->GetYaxis()->SetTitle("#frac{d#sigma}{d#Omega} (a.u.)");*/

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

   /*for (int i = 0; i < gL1->GetN(); i++) {
      double x, y;
      gL1->GetPoint(i, x, y);
      cout << i << "  " << x << "  " << y << endl;
   }*/
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
