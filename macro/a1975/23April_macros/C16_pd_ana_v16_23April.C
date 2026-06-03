#include <fstream>
#include <iostream>

double Ebin_max = 9.;
double Ebin_min = -1.;
int NumberBins = 100; // 140
int NumberBinsAux = 200;

class SpectralModelSeg {
public:
   TGraph *graphPS; // Phase-space TGraph

   SpectralModelSeg(TGraph *g) : graphPS(g) {}

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

/*double TotalVerticalError(TH1F* hExp, TGraphErrors* gTheory, double scale) {
   double errorSum = 0.0;

   for (int i = 1; i <= hExp->GetNbinsX(); ++i) {
      double x = hExp->GetBinCenter(i);
      double y_exp = hExp->GetBinContent(i);
      double y_theory = gTheory->Eval(x);

      double delta = std::abs(scale * y_exp - y_theory);
      errorSum += delta;
   }

   return errorSum;
}*/

/*Bool_t compareEventName(std::string &getname, std::string &fribname)
{
   // Parsing FRIB event number
   std::regex fribregex("evt(\\d+)_\\d+");
   std::string result = std::regex_replace(fribname, fribregex, "$1\n");

   int fribnumber;
   std::istringstream iss(result);
   while (iss >> fribnumber) {
      // std::cout << fribnumber << std::endl;
   }

   // Parsing GET event name
   std::regex getregex("evt(\\d+)_data");
   result = std::regex_replace(getname, getregex, "$1\n");

   int getnumber;
   std::istringstream isss(result);
   while (isss >> getnumber) {
      // std::cout << getnumber << std::endl;
   }

   return (fribnumber == getnumber) ? 1 : 0;

   return 0;
}*/

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
void C16_pd_ana_v16_23April()
{
   bool guardar_en_pdf = false; // ← cambia a false si quieres solo verlos en pantalla
   gROOT->ProcessLine(".X /home/georgina/fair_install/ATTPCROOTv2/macro/a1975/myStyle.C");

   // Activar modo batch si estás guardando en PDF
   if (guardar_en_pdf) {
      gROOT->SetBatch(kTRUE); // ← esto evita que se abran ventanas
   } else {
      gROOT->SetBatch(kFALSE); // ← esto permite ver los canvas en pantalla
   }

   TH2F *Ang_Ener_Corr = new TH2F("Ang_Ener_Corr", "Ang_Ener_Corr", 720, 10, 40, 1000, 0, 60.0);
   TH2F *Ang_Ener_Cal = new TH2F("Ang_Ener_Cal", "Ang_Ener_Cal_JR", 720, 10, 40, 1000, 0, 60.0);

   TH2F *Ebeam_test = new TH2F("Ebeam_test", "Ebeam_test", 1000, -2, 10, 60, 0, 300);

   TH2F *ELossvsBrho = new TH2F("ELossvsBrho", "ELossvsBrho", 4000, 0, 25000, 1000, 0, 4);
   TH2F *dedxvsBrho = new TH2F("dedxvsBrho", "dedxvsBrho", 4000, 0, 10000, 1000, 0, 4);
   TH2F *hVxVy = new TH2F("hVxVy", "hVxVy", 1000, 0, 4, 1000, 0, 4);

   TH1F *hGS_AngularDistr = new TH1F("hGS_AngularDistr", "Ground State Angular Distribution", 30, 0, 180);
   TH1F *h740_AngularDistr = new TH1F("h740_AngularDistr", "740keV Angular Distribution", 30, 0, 180);
   TH1F *h3103_AngularDistr = new TH1F("h3103_AngularDistr", "3103keV Angular Distribution", 30, 0, 180);
   TH1F *h4780_AngularDistr = new TH1F("h4780_AngularDistr", "4780keV Angular Distribution", 30, 0, 180);
   TH1F *h6841_AngularDistr = new TH1F("h6841_AngularDistr", "6841keV Angular Distribution", 30, 0, 180);

   auto *hex = new TH1F("hex", "hex", NumberBins, Ebin_min, Ebin_max);
   auto *QvsEb = new TH2F("QvsEb", "QvsEb", 1000, -5, 15, 100, 0, 300);
   auto *QvsZpos = new TH2F("QvsZpos", "QvsZpos", 1000, -10, 50, 200, -100, 100);
   auto *hexCorr = new TH1F("hexCorr", "", NumberBins, Ebin_min, Ebin_max);
   auto *hexCorr1 = new TH1F("hexCorr1", "C16(p,d)", NumberBins, Ebin_min, Ebin_max);
   auto *hexCorr2 = new TH1F("hexCorr2", "C16(p,d)", NumberBins, Ebin_min, Ebin_max);
   // auto *hexCorr_seg = new TH1F("hexCorr_seg", "hexCorr_seg", NumberBins, Ebin_min, Ebin_max);

   TH1F *hGS_ExEnergy =
      new TH1F("hGS_ExEnergy", "Excitation Energy (Ground State Cut)", NumberBins, Ebin_min, Ebin_max);
   TH1F *h740_ExEnergy = new TH1F("h740_ExEnergy", "Excitation Energy (740keV Cut)", NumberBins, Ebin_min, Ebin_max);
   TH1F *h3103_ExEnergy = new TH1F("h3013_ExEnergy", "Excitation Energy (3103keV Cut)", NumberBins, Ebin_min, Ebin_max);
   TH1F *h4780_ExEnergy = new TH1F("h4780_ExEnergy", "Excitation Energy (4780KeV Cut)", NumberBins, Ebin_min, Ebin_max);
   TH1F *h6841_ExEnergy = new TH1F("h6841_ExEnergy", "Excitation Energy (6841KeV Cut)", NumberBins, Ebin_min, Ebin_max);

   auto *AngDistr = new TH1F("Ang_Distr", "Ang_Distr", 128, 0, 120);
   auto *AngDistrCM = new TH1F("Ang_Distr_CM", "Ang_Distr_CM", NumberBins, 0, 180);
   auto *ExvsZpos = new TH2F("ExvsZpos", "ExvsZpos", 1000, -5, 10, 200, -5, 80);
   // auto *ExvsTrackLength = new TH2F("ExvsTrackLength", "ExvsTrackLength", 1000, -5, 15, 200, -20, 150);
   auto *ExCorrvsZpos = new TH2F("ExCorrvsZpos", "ExCorrvsZpos", 1000, -5, 10, 200, -5, 80);
   auto *KineticEnergy = new TH1F("KineticEnergy", "KineticEnergy", 100, 0, 100);

   // After your existing graphs loading loop, add a new 2D histogram
   TH2F *h_residual = new TH2F("h_residual", "#DeltaE vs E_{ej} for GS events", 100, 0, 60, // E_ej axis
                               100, -20, 20);                                               // E_meas - E_theory

   TH2F *h_dtheta = new TH2F("h_dtheta", "#Delta#theta vs E_{ej} for GS events", 100, 0, 60, // E_ej
                             100, -0.1, 0.1);                                                // theta residual in rad

   TH2F *h_thetacorr_vs_Eej = new TH2F("h_thetacorr_vs_Eej", "#theta_{lab,corr} vs E_{ej}", 100, 0, 14, // E_ej (MeV)
                                       100, 0, 60); // theta_lab_corr (deg)

   TH2F *h_theta_vs_Eej = new TH2F("h_theta_vs_Eej", "#theta_{lab} vs E_{ej}", 100, 0, 14, // E_ej (MeV)
                                   100, 0, 60);                                            // theta_lab (deg)

   TH1F *h_PS_1n = new TH1F("h_PS_1n", "h_PS_1n", NumberBins, Ebin_min, Ebin_max);

   TH2F *h_GS_zEx = new TH2F("h_GS_zEx", "GS: Ex vs z", 50, 0, 1.0, // z en metros
                             50, -0.5, 0.5);                        // Ex en MeV

   // Perfil de Ex vs z para el 2nd excited (Ex entre 2.5 y 4.5 MeV)
   TH2F *h_2nd_zEx = new TH2F("h_2nd_zEx", "1st Excited: Ex vs z", 50, 0, 1.0, 50, 0.5, 1.5);

   /*TH2F *h_Ecal = new TH2F("h_Ecal", "E_medido vs E_teorico (GS)",
                          200, 0, 60,   // E_teorico (MeV)
                          200, 0, 60);  // E_medido (MeV)

   // Histogram to project the angle at the GS kinematic elbow
   TH1F *h_theta_elbow = new TH1F("h_theta_elbow", "GS Elbow Angle Projection",
                                40, 34, 40);*/

   /*auto *hredchi2 = new TH1F("redchi2", "redchi2", 1000, 0, 0.0001);
   auto *hbredchi2 = new TH1F("bredchi2", "bredchi2", 1000, 0, 5);*/

   auto *hexvstheta = new TH2F("hexVStheta", "hexVStheta", 100, -2, 10, 100, 0, 50);

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

   std::vector<TString> filenames;
   // Declare hHex at the beginning of your macro or function
   std::vector<TH1F *> hHex(5); // 9
   // Assuming hHex is already defined and filled
   std::vector<TF1 *> fExSpectra_vec(hHex.size(), nullptr);

   // -----------------------------KINEMATICS FOR DIFFERENT EXCITATION ENERGIES

   std::vector<std::string> files = {"C16_pd_C15_gs_Ebeam11_5.txt", "C16_pd_C15_740keV_Ebeam11_5.txt",
                                     "C16_pd_C15_3103keV_Ebeam11_5.txt", "C16_pd_C15_4780keV_Ebeam11_5.txt",
                                     "C16_pd_C15_6841keV_Ebeam11_5.txt"};
   std::vector<std::string> labels = {"Ground State", "1st Excited State (740keV)", "2nd Excited State (3103keV)",
                                      "3rd Excited State (4780keV)", "4th Excited State(6841keV)"};
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

   std::vector<std::pair<int, int>> angularBins = {{20, 30}, {30, 40}, {40, 50}, {50, 70}, {70, 90}};

   for (size_t i = 0; i < angularBins.size(); ++i) {
      int thetaMin = angularBins[i].first;
      int thetaMax = angularBins[i].second;
      TString histTitle = Form("Excitation Energy (%.1d deg - %.1d deg)", thetaMin, thetaMax);
      TString histName = Form("hex_%d_%d", thetaMin, thetaMax);
      hHex[i] = new TH1F(histName, histTitle, NumberBins, Ebin_min, Ebin_max); // 90, -5, 14
   }

   // ELoss tables.
   // AtTools::AtELossTable *elossTableH2 = new AtTools::AtELossTable();
   // elossTableH2->LoadSrimTable("StoppingPower_SRIM_C16_H2.txt"); //SRIM no me va.
   // elossTableH2->LoadLiseTable("StoppingPower_C16_H2.txt", 2.0158,3.3084e-5);

   double densityH2 = 3.3084e-5; // g/cm³
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

         double dist3D = TMath::Sqrt(vx_pos * vx_pos + vy_pos * vy_pos + zPos * zPos) * 100.0; // mm
         // Double_t Ebeam_at_z = elossH2.GetEnergy(Ebeam_buff, zPos * 100); //
         Double_t Ebeam_at_z =
            elossH2.GetEnergy(Ebeam_buff, dist3D); // Usar la distancia 3D para la corrección de energía, MeV/mm

         // Corrección cinemática
         double kethe = 13.;
         double theta_lab_corr_JR =
            (theta -
             (2.0 * TMath::Pi() / 4000) * (E_ej - kethe)); // theta: rad; theta_lab_corr: rad; E_ej-kethe: MeV 29.5

         // double theta_lab_corr_JR = theta - 0.00155 * (E_ej - 15.6); // Corrección de ángulo en radianes
         h_thetacorr_vs_Eej->Fill(E_ej, theta_lab_corr_JR * TMath::RadToDeg());
         h_theta_vs_Eej->Fill(E_ej, theta * TMath::RadToDeg());

         double theta_lab_corr = theta; // Por ahora sin corrección de ángulo, solo para probar la implementación de la
                                        // corrección de energía en el cálculo de Ex y theta_cm
         auto [ex_energy_corr, theta_cm_corr] = kine_2b(m_C16, m_p, m_b, m_B, Ebeam_at_z, theta_lab_corr, E_ej);

         auto [ex_energy_corr_JR, theta_cm_corr_JR] = kine_2b(m_C16, m_p, m_b, m_B, Ebeam_at_z, theta_lab_corr_JR,
                                                              E_ej); // energies: MeV, angles: radians

         // Fill uncorrected histogram
         hex->Fill(ex_energy);
         ExvsZpos->Fill(ex_energy_corr, zPos * 100.0); // MeV, cm
         KineticEnergy->Fill(E_ej); // Usar energía calibrada para el histograma de energía cinética

         Ang_Ener_Corr->Fill(theta_lab_corr_JR * TMath::RadToDeg(),
                             E_ej); // theta lab!! -> I still have to implement the correction of catima?
         Ang_Ener_Cal->Fill(theta_lab_corr_JR * TMath::RadToDeg(), E_ej); // calibrado

         double theta_deg = theta * TMath::RadToDeg();

         // Fill corrected histogram
         if (zPos * 100 > 2.0 && zPos * 100 < 60.0 && E_ej < 14.0) { // cm y MeV (zPos is in meters, E_ej is in MeV)
            ExCorrvsZpos->Fill(ex_energy_corr_JR, zPos * 100.0);
            hexCorr->Fill(ex_energy_corr_JR);
            hexCorr2->Fill(ex_energy_corr_JR);
         }

         if (zPos * 100 > 2.0 && zPos * 100 < 60.0 && E_ej < 20.0) { // cm y MeV
            hexCorr1->Fill(ex_energy_corr_JR);
         }

         /*if(ex_energy_corr > -0.3 && ex_energy_corr < 0.4)
            h_GS_zEx->Fill(zPos*100, ex_energy_corr);
         if(ex_energy_corr > 0.7 && ex_energy_corr < 1.1)
            h_2nd_zEx->Fill(zPos*100.0, ex_energy_corr);*/

         //-----------------------------------------------------------------------------------
         // Cross section selection with energies

         /*Double_t theta_deg_lower = 25.0;
         Double_t theta_deg_upper = 100.0;

         Double_t energyCutGS_lower = -0.6;
         Double_t energyCutGS_upper = 0.4;

         if (ex_energy_corr > energyCutGS_lower && ex_energy_corr < energyCutGS_upper &&
             theta_cm_corr > theta_deg_lower && theta_cm_corr < theta_deg_upper) {
            Double_t theta_deg = theta_cm_corr;
            hGS_AngularDistr->Fill(theta_deg);
            hGS_ExEnergy->Fill(ex_energy_corr);

         }

         if(ex_energy_corr > 0.4 && ex_energy_corr< 1.2 &&
             theta_cm_corr > theta_deg_lower && theta_cm_corr < theta_deg_upper) {
            Double_t theta_deg = theta_cm_corr;
            h740_AngularDistr->Fill(theta_deg);
            h740_ExEnergy->Fill(ex_energy_corr);
         }

         if(ex_energy_corr > 3.0 && ex_energy_corr< 3.8 &&
             theta_cm_corr > theta_deg_lower && theta_cm_corr < theta_deg_upper) {
            Double_t theta_deg = theta_cm_corr;
            h3103_AngularDistr->Fill(theta_deg);
            h3103_ExEnergy->Fill(ex_energy_corr);
         }


         if(ex_energy_corr > 4.0 && ex_energy_corr< 5.5 &&
             theta_cm_corr > theta_deg_lower && theta_cm_corr < theta_deg_upper) {
            Double_t theta_deg = theta_cm_corr;
            h4780_AngularDistr->Fill(theta_deg);
            h4780_ExEnergy->Fill(ex_energy_corr);
         }

         if(ex_energy_corr > 6.0 && ex_energy_corr< 8.0 &&
             theta_cm_corr > theta_deg_lower && theta_cm_corr < theta_deg_upper) {
            Double_t theta_deg = theta_cm_corr;
            h6841_AngularDistr->Fill(theta_deg);
            h6841_ExEnergy->Fill(ex_energy_corr);
         }*/

         //-----------------------------------------------------------------------------------

         // Histograms

         Double_t vx = TMath::Sin(theta) * TMath::Sqrt(ke);
         Double_t vy = TMath::Cos(theta) * TMath::Sqrt(ke);

         hVxVy->Fill(vx, vy);

         AngDistr->Fill(theta * TMath::RadToDeg());
         AngDistrCM->Fill(theta_cm);
         hexvstheta->Fill(ex_energy_corr, theta * TMath::RadToDeg());
         // ExvsTrackLength->Fill(ex_energy_corr, arclength);

         for (size_t i = 0; i < angularBins.size(); ++i) {
            int thetaMin = angularBins[i].first;
            int thetaMax = angularBins[i].second;

            if (theta_cm > thetaMin && theta_cm <= thetaMax && ke < 25.) {
               hHex[i]->Fill(ex_energy_corr);

               break; // Only fill one bin per event
            }
         }

         // tEvents->Fill();
      } // events
   }    // Files

   AngDistrCM->Divide(new TF1("sin", "sin(x * TMath::DegToRad())", 0, 180));

   /*TCanvas *c_prof_GS = new TCanvas("c_prof_GS", "GS: Ex vs z", 800, 600);
   c_prof_GS->cd();
   TProfile *prof_GS = h_GS_zEx->ProfileX("prof_GS");
   prof_GS->SetTitle("Ground State: Ex vs z");
   prof_GS->GetXaxis()->SetTitle("z (m)");
   prof_GS->GetYaxis()->SetTitle("Ex (MeV)");
   prof_GS->SetMarkerStyle(20);

   // 1. Fit GS with a linear function ("pol1") and extract the TF1
   prof_GS->Fit("pol1", "Q"); // "Q" makes the terminal output quiet
   TF1 *fit_GS = prof_GS->GetFunction("pol1");
   fit_GS->SetName("fit_GS"); // Rename to avoid memory conflicts
   fit_GS->SetLineColor(kRed);
   fit_GS->SetLineWidth(2);
   prof_GS->Draw("EP"); // Draw profile with its fit

   TCanvas *c_prof_2nd = new TCanvas("c_prof_2nd", "2nd Excited: Ex vs z", 800, 600);
   c_prof_2nd->cd();
   TProfile *prof_2nd = h_2nd_zEx->ProfileX("prof_2nd");
   prof_2nd->SetTitle("2nd Excited: Ex vs z");
   prof_2nd->GetXaxis()->SetTitle("z (m)");
   prof_2nd->GetYaxis()->SetTitle("Ex (MeV)");
   prof_2nd->SetMarkerStyle(20);

   prof_2nd->Fit("pol1", "Q");
   TF1 *fit_2nd = prof_2nd->GetFunction("pol1");
   fit_2nd->SetName("fit_2nd");
   fit_2nd->SetLineColor(kBlue);
   fit_2nd->SetLineWidth(2);
   prof_2nd->Draw("EP");


   TCanvas *c_ExvsZpos = new TCanvas( "ExCorrvsZpos", "Excitation Energy vs z position and track length", 1200, 800);
   c_ExvsZpos->cd();
   ExCorrvsZpos->GetXaxis()->SetTitle("Excitation Energy (MeV)");
   ExCorrvsZpos->GetYaxis()->SetTitle("z (cm)");
   gPad->SetRightMargin(0.20);
   ExCorrvsZpos->Draw("zcol");

   // --- Ground State Line ---
   // Choose the z-range in meters (based on your profile plot)
   double z_start_m = 0.0;
   double z_end_m   = 0.8; // Looking at your data, it stops around 0.8m

   // Evaluate Ex at these z positions
   double Ex_GS_start = fit_GS->Eval(z_start_m);
   double Ex_GS_end   = fit_GS->Eval(z_end_m);

   // Create a line: TLine(x1, y1, x2, y2)
   // Remember to swap axes: x = Ex, y = z (in cm!)
   TLine *line_GS = new TLine(Ex_GS_start, z_start_m * 100.0,
                              Ex_GS_end,   z_end_m * 100.0);
   line_GS->SetLineColor(kRed);
   line_GS->SetLineWidth(3);
   line_GS->Draw("same");

   // --- 2nd Excited State Line ---
   double Ex_2nd_start = fit_2nd->Eval(z_start_m);
   double Ex_2nd_end   = fit_2nd->Eval(z_end_m);

   cout << "fit params 2nd: " << fit_2nd->GetParameter(0) << " " << fit_2nd->GetParameter(1) << endl;
   cout << "fit params GS: " << fit_GS->GetParameter(0) << " " << fit_GS->GetParameter(1) << endl;

   TLine *line_2nd = new TLine(Ex_2nd_start, z_start_m * 100.0,
                               Ex_2nd_end,   z_end_m * 100.0);
   line_2nd->SetLineColor(kBlue);
   line_2nd->SetLineWidth(3);
   line_2nd->Draw("same");


// Change the existing h_residual to show the 2D distribution clearly
TCanvas *c_residual = new TCanvas("c_residual", "Kinematic Residuals GS", 800, 600);
c_residual->cd();
h_residual->Draw("colz"); // colz shows density, easier to see the trend

*/

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
   int nfound = sp->Search(hexCorr, 2, "", 0.05); // 2 = sigma of smoothing, last = threshold
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
   double m1 = sorted_peaks[0];
   double m2 = sorted_peaks[1];
   double m3 = sorted_peaks[2];
   double m4 = sorted_peaks[3];
   double m5 = sorted_peaks[4];

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
   f2g->Draw("same"); // <--- REQUIRED so the fit curve is drawn

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

   SpectralModelSeg *model = new SpectralModelSeg(graphPS);
   TF1 *fModel = new TF1("fModel", model, -0.5, 9.0, 16, // number fitted parameters
                         "SpectralModelSeg");

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
   // fModel->SetParLimits(15, 5e-4, 1e-2);  // o un límite razonable

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

   hexCorr2->Sumw2(); // activa almacenamiento de errores
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
   // hexCorr2->Draw("E1"); //E1

   TLine *vline0 = new TLine(1.218, 0, 1.218, 390); // línea vertical en x=1.218 MeV;
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
   TCanvas *c_combined = new TCanvas("c_combined", "Combined Analysis", 1200, 800);
   c_combined->Divide(2, 1); // 2 columns, 2 rows

   // --- PAD 1: Ang_Ener_Corr (uncorrected kinematics) ---
   c_combined->cd(1);
   Ang_Ener_Corr->Draw("col");
   Ang_Ener_Corr->GetXaxis()->SetTitle("#theta_{lab} (deg)");
   Ang_Ener_Corr->GetYaxis()->SetTitle("Kinetic Energy (MeV)");

   /*TLine *line1 = new TLine(Ang_Ener_Corr->GetXaxis()->GetXmin(), 20,
                             Ang_Ener_Corr->GetXaxis()->GetXmax(), 20);
   line1->SetLineColor(kRed); line1->SetLineStyle(2); line1->SetLineWidth(2);
   line1->Draw("SAME");*/
   for (size_t i = 0; i < graphs.size(); i++) {
      if (graphs[i] == nullptr || graphs[i]->GetN() == 0)
         continue;
      graphs[i]->Draw("L SAME");
   }
   auto leg0 = new TLegend(0.35, 0.7, 0.9, 0.9);
   for (size_t i = 0; i < graphs.size(); i++) {
      if (graphs[i] == nullptr)
         continue;
      leg0->AddEntry(graphs[i], labels[i].c_str(), "l");
   }
   leg0->Draw();

   // --- PAD 2: Ang_Ener_Cal (corrected angle kinematics) ---
   c_combined->cd(2);
   Ang_Ener_Cal->Draw("col");
   Ang_Ener_Cal->GetXaxis()->SetTitle("#theta_{lab} (deg)");
   Ang_Ener_Cal->GetYaxis()->SetTitle("Kinetic Energy (MeV)");

   /*TLine *line2 = new TLine(Ang_Ener_Cal->GetXaxis()->GetXmin(), 20,
                             Ang_Ener_Cal->GetXaxis()->GetXmax(), 20);
   line2->SetLineColor(kRed); line2->SetLineStyle(2); line2->SetLineWidth(2);
   line2->Draw("SAME");*/

   for (size_t i = 0; i < graphs.size(); i++) {
      if (graphs[i] == nullptr || graphs[i]->GetN() == 0)
         continue;
      graphs[i]->Draw("L SAME");
   }
   auto leg2 = new TLegend(0.55, 0.7, 0.9, 0.9);
   for (size_t i = 0; i < graphs.size(); i++) {
      if (graphs[i] == nullptr)
         continue;
      leg2->AddEntry(graphs[i], labels[i].c_str(), "l");
   }
   leg2->Draw();

   c_combined->Update();

   //---------------------combined 2
   TCanvas *c_combined2 = new TCanvas("c_combined2", "Combined Analysis", 1200, 800);
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

   auto leg12 = new TLegend(0.28, 0.75, 0.9, 0.9);
   for (size_t i = 0; i < graphs.size(); i++) {
      if (graphs[i] == nullptr)
         continue;
      leg12->AddEntry(graphs[i], labels[i].c_str(), "l");
   }
   leg12->Draw();

   // --- PAD 3: hexCorr (uncorrected excitation energy) ---
   c_combined2->cd(2);
   hexCorr1->GetXaxis()->SetTitle("Excitation Energy (MeV)");
   hexCorr1->GetYaxis()->SetTitle("Counts");
   hexCorr1->SetLineColor(kGreen + 2);
   hexCorr1->SetLineWidth(2);
   hexCorr1->Draw("HIST");
   hexCorr2->SetLineColor(kBlue);
   hexCorr2->Draw("same HIST");

   auto leg42 = new TLegend(0.45, 0.8, 0.9, 0.9);
   leg42->AddEntry(hexCorr2, "hexCorr2, ke < 14 MeV", "l");
   leg42->AddEntry(hexCorr1, "hexCorr1, ke < 20 MeV", "l");
   leg42->Draw("same");

   c_combined2->Update();

   TCanvas *c_thetacorr = new TCanvas("c_thetacorr", "#theta_{corr} vs E_{ej}", 800, 600);
   c_thetacorr->Divide(2, 1);
   c_thetacorr->cd(1);
   h_thetacorr_vs_Eej->Draw("colz");
   h_thetacorr_vs_Eej->GetXaxis()->SetTitle("E_{ej} (MeV)");
   h_thetacorr_vs_Eej->GetYaxis()->SetTitle("#theta_{lab,corr} (deg)");

   c_thetacorr->cd(2);
   h_theta_vs_Eej->Draw("colz");
   h_theta_vs_Eej->GetXaxis()->SetTitle("E_{ej} (MeV)");
   h_theta_vs_Eej->GetYaxis()->SetTitle("#theta_{lab} (deg)");

   // Excitation energy spectrum with fits --------------------------------------

   /* TCanvas* hGS_energy = new TCanvas("hGS_ExEnergy", "Ground State Excitation Energy", 1000, 600);
    hGS_energy->cd();
    hGS_ExEnergy->GetXaxis()->SetTitle("Excitation Energy (MeV)");
    hGS_ExEnergy->GetYaxis()->SetTitle("Counts");
    hGS_ExEnergy->SetLineColor(kBlue+2);
    hGS_ExEnergy->SetLineWidth(2);
    hGS_ExEnergy->Draw();

    TCanvas* cAngGS = new TCanvas("cAngGS", "Ground State Angular Distribution", 1000, 600);
    cAngGS->cd();
    //gPad->SetLogy();
    // Histograma experimental primero
    hGS_AngularDistr->Sumw2();
    hGS_AngularDistr->Scale(0.5);  // ← se aplica correctamente
    hGS_AngularDistr->SetLineColor(kBlue+2);
    hGS_AngularDistr->SetLineWidth(2);
    hGS_AngularDistr->GetXaxis()->SetTitle("#theta_cm (deg)");
    hGS_AngularDistr->GetYaxis()->SetTitle("Counts");
    hGS_AngularDistr->Draw("E1");  // ← define el marco


    TCanvas* h740_energy = new TCanvas("h740_ExEnergy", "Ground State Excitation Energy", 1000, 600);
    h740_energy->cd();
    h740_ExEnergy->GetXaxis()->SetTitle("Excitation Energy (MeV)");
    h740_ExEnergy->GetYaxis()->SetTitle("Counts");
    h740_ExEnergy->SetLineColor(kBlue+2);
    h740_ExEnergy->SetLineWidth(2);
    h740_ExEnergy->Draw();

    TCanvas* cAng740 = new TCanvas("cAng740", "740keV Angular Distribution", 1000, 600);
    cAng740->cd();
    //gPad->SetLogy();
    h740_AngularDistr->Sumw2();
    h740_AngularDistr->Scale(0.5);  // ← se aplica correctamente
    h740_AngularDistr->SetLineColor(kBlue+2);
    h740_AngularDistr->SetLineWidth(2);
    h740_AngularDistr->GetXaxis()->SetTitle("#theta_cm (deg)");
    h740_AngularDistr->GetYaxis()->SetTitle("Counts");
    h740_AngularDistr->Draw("E1");  // ← define el marco

    TCanvas* h3103_energy = new TCanvas("h3103_ExEnergy",  "Excitation Energy", 1000, 600);
    h3103_energy->cd();
    h3103_ExEnergy->GetXaxis()->SetTitle("Excitation Energy (MeV)");
    h3103_ExEnergy->GetYaxis()->SetTitle("Counts");
    h3103_ExEnergy->SetLineColor(kBlue+2);
    h3103_ExEnergy->SetLineWidth(2);
    h3103_ExEnergy->Draw();

    TCanvas* cAng3103 = new TCanvas("cAng3103", "3103keV Angular Distribution", 1000, 600);
    cAng3103->cd();
    //gPad->SetLogy();
    h3103_AngularDistr->Sumw2();
    h3103_AngularDistr->Scale(0.5);  // ← se aplica correctamente
    h3103_AngularDistr->SetLineColor(kBlue+2);
    h3103_AngularDistr->SetLineWidth(2);
    h3103_AngularDistr->GetXaxis()->SetTitle("#theta_cm (deg)");
    h3103_AngularDistr->GetYaxis()->SetTitle("Counts");
    h3103_AngularDistr->Draw("E1");  // ← define el marco

    TCanvas* h4780_energy = new TCanvas("h4780_ExEnergy",  "Excitation Energy", 1000, 600);
    h4780_energy->cd();
    h4780_ExEnergy->GetXaxis()->SetTitle("Excitation Energy (MeV)");
    h4780_ExEnergy->GetYaxis()->SetTitle("Counts");
    h4780_ExEnergy->SetLineColor(kBlue+2);
    h4780_ExEnergy->SetLineWidth(2);
    h4780_ExEnergy->Draw();

    TCanvas* cAng4780 = new TCanvas("cAng4780", "4780keV Angular Distribution", 1000, 600);
    cAng4780->cd();
    //gPad->SetLogy();
    h4780_AngularDistr->Sumw2();
    h4780_AngularDistr->Scale(0.5);  // ← se aplica correctamente
    h4780_AngularDistr->SetLineColor(kBlue+2);
    h4780_AngularDistr->SetLineWidth(2);
    h4780_AngularDistr->GetXaxis()->SetTitle("#theta_cm (deg)");
    h4780_AngularDistr->GetYaxis()->SetTitle("Counts");
    h4780_AngularDistr->Draw("E1");  // ← define el marco

 */
   //-----------------------------------------------------------------------------------------
   /*auto* g1 = new TGraphErrors("/home/georgina/twofnr/21.groundState", "%lg %lg");
   double bestScale = 1.0;
   double minError = 1e9;

   for (double testScale = 0.01; testScale <= 5.0; testScale += 0.01) {
      double error = TotalVerticalError(hGS_AngularDistr, g1, testScale);
      if (error < minError) {
         minError = error;
         bestScale = testScale;
      }
   }
   std::cout << "Mejor escala por distancia vertical: " << bestScale << std::endl;

   TCanvas* cAngGS_overlay = new TCanvas("cAngGS_overlay", "Ground State Angular Distribution Overlay", 1000, 600);
   cAngGS_overlay->cd();

   cAngGS_overlay->SetLogy();
    // Histograma experimental primero
   hGS_AngularDistr->Sumw2();
   hGS_AngularDistr->SetLineColor(kRed);
   hGS_AngularDistr->SetLineWidth(2);
   hGS_AngularDistr->SetMarkerStyle(20);       // marcador redondo sólido
   hGS_AngularDistr->SetMarkerColor(kRed);     // color del marcador
   hGS_AngularDistr->SetMarkerSize(1.2);       // tamaño del marcador

   hGS_AngularDistr->GetXaxis()->SetTitle("#theta_CM (deg)");
   hGS_AngularDistr->GetYaxis()->SetTitle("d#sigma/d#Omega (mb/sr)");
   hGS_AngularDistr->Scale(bestScale);  // ← se aplica el ajuste automático
   hGS_AngularDistr->Draw("E1");

   g1->SetStats(0);
   g1->SetLineWidth(2);
   g1->SetMarkerStyle(20);
   g1->GetYaxis()->SetRangeUser(0.1, 1000);
   g1->SetMarkerColor(kBlack);
   g1->SetLineColor(kBlack);
   g1->Draw("same P L");  // "P" for markers, "L" for line

   auto *legend8 = new TLegend(0.75, 0.75, 0.88, 0.88);
   legend8->AddEntry(hGS_AngularDistr, "GS: s_{1/2}^{+}", "p");  // use h740_AngularDistr here
   legend8->SetTextSize(0.04);
   legend8->AddEntry(g1, "twofnr", "lp");
   legend8->Draw();*/

   //----------------------------------------------------------------------------

   /*auto* g1= new TGraphErrors("/home/georgina/twofnr/21.1stExcited", "%lg %lg");
   double bestScale = 1.0;
   double minError = 1e9;

   for (double testScale = 0.01; testScale <= 5.0; testScale += 0.01) {
      double error = TotalVerticalError(hGS_AngularDistr, g1, testScale);
      if (error < minError) {
         minError = error;
         bestScale = testScale;
      }
   }
   std::cout << "Mejor escala por distancia vertical: " << bestScale << std::endl;*/

   /*TCanvas* cAng740_overlay = new TCanvas("cAng740_overlay", "c", 1000, 600);
   cAng740_overlay->cd();

   cAng740_overlay->SetLogy();
    // Histograma experimental primero
   h740_AngularDistr->Sumw2();
   h740_AngularDistr->SetLineColor(kRed);
   h740_AngularDistr->SetLineWidth(2);
   h740_AngularDistr->SetMarkerStyle(20);       // marcador redondo sólido
   h740_AngularDistr->SetMarkerColor(kRed);     // color del marcador
   h740_AngularDistr->SetMarkerSize(1.2);       // tamaño del marcador

   h740_AngularDistr->GetXaxis()->SetTitle("#theta_CM (deg)");
   h740_AngularDistr->GetYaxis()->SetTitle("d#sigma/d#Omega (mb/sr)");
   h740_AngularDistr->Scale(bestScale);  // ← se aplica el ajuste automático
   h740_AngularDistr->Draw("E1");

   g1->SetStats(0);
   g1->SetLineWidth(2);
   g1->SetMarkerStyle(20);
   g1->GetYaxis()->SetRangeUser(0.1, 1000);
   g1->SetMarkerColor(kBlack);
   g1->SetLineColor(kBlack);
   g1->Draw("same P L");  // "P" for markers, "L" for line

   auto *legend9 = new TLegend(0.75, 0.75, 0.88, 0.88);
   legend9->AddEntry(h740_AngularDistr, "1st: d_{5/2}^{+}", "p");  // use h740_AngularDistr here
   legend9->SetTextSize(0.04);
   legend9->AddEntry(g1, "twofnr", "lp");
   legend9->Draw();

   cout << "best scale 1st excited: " << bestScale << endl;*/
   //----------------------------------------------------------------------------
   /*auto* g2= new TGraphErrors("/home/georgina/twofnr/21.2ndExcitedP", "%lg %lg");
   double bestScale = 1.0;
   double minError = 1e9;

   for (double testScale = 0.01; testScale <= 5.0; testScale += 0.01) {
      double error = TotalVerticalError(hGS_AngularDistr, g2, testScale);
      if (error < minError) {
         minError = error;
         bestScale = testScale;
      }
   }*/

   /*TCanvas* cAng3103_overlay = new TCanvas("cAng3103_overlay", "c", 1000, 600);
   cAng3103_overlay->cd();

   cAng3103_overlay->SetLogy();
    // Histograma experimental primero
   h3103_AngularDistr->Sumw2();
   h3103_AngularDistr->SetLineColor(kRed);
   h3103_AngularDistr->SetLineWidth(2);
   h3103_AngularDistr->SetMarkerStyle(20);       // marcador redondo sólido
   h3103_AngularDistr->SetMarkerColor(kRed);     // color del marcador
   h3103_AngularDistr->SetMarkerSize(1.2);       // tamaño del marcador

   h3103_AngularDistr->GetXaxis()->SetTitle("#theta_CM (deg)");
   h3103_AngularDistr->GetYaxis()->SetTitle("d#sigma/d#Omega (mb/sr)");
   h3103_AngularDistr->Scale(0.05);  // ← se aplica el ajuste automático
   h3103_AngularDistr->Draw("E1");

   g2->SetLineColor(kRed);
   g2->SetStats(0);
   g2->SetLineWidth(2);
   g2->SetMarkerStyle(20);
   g2->GetYaxis()->SetRangeUser(0.1, 1000);
   g2->SetMarkerColor(kBlack);
   g2->Draw("same P");

   auto *legend9 = new TLegend(0.75, 0.75, 0.88, 0.88);
   legend9->AddEntry(h740_AngularDistr, "2nd: p_{1/2}^{-}", "p");  // use h740_AngularDistr here
   legend9->SetTextSize(0.04);
   legend9->AddEntry(g2, "twofnr", "lp");
   legend9->Draw();*/

   //--------------------------------------------------------------------------------------------------

   /*TCanvas *c_ExenerCorr = new TCanvas("ExenerCorr", "Excited Energy spectra corrected", 1200, 800);
   c_ExenerCorr->cd();
   c_ExenerCorr->Divide(2, 1);
   c_ExenerCorr->cd(1);
   hexCorr->Draw("E1");
   c_ExenerCorr->cd(2);
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
   AngDistrCM->GetYaxis()->SetTitle("#frac{d#sigma}{d#Omega} (a.u.)");*/

   /*TCanvas *c_ExvsZpos = new TCanvas( "ExCorrvsZpos", "Excitation Energy vs z position and track length", 800, 1200);
   c_ExvsZpos->cd();
   c_ExvsZpos->Divide(2, 1);
   c_ExvsZpos->cd(1);
   ExvsZpos->Draw("zcol");
   ExvsZpos->GetXaxis()->SetTitle("Excitation Energy (MeV)");
   ExvsZpos->GetYaxis()->SetTitle("z (cm)");
   c_ExvsZpos->cd(2);
   ExCorrvsZpos->Draw("zcol");*/
   /*ExvsTrackLength->Draw("zcol");
   ExvsTrackLength->GetXaxis()->SetTitle("Excitation Energy (MeV)");
   ExvsTrackLength->GetYaxis()->SetTitle("Track Length (cm)");*/

   /* TCanvas *kin = new TCanvas( "kin", "kin", 1200, 800);
    KineticEnergy->Sumw2();
    kin->cd();
    KineticEnergy->Draw("E1");*/

   // Para dibujar líneas verticales en ese plot 2D puedes usar TLine. Añade esto después de dibujar el histograma:
   // Líneas verticales para cada estado
   /* TLine *l_GS   = new TLine(0.0,  0, 0.0,  80);   // Ground State
    TLine *l_1st  = new TLine(0.74, 0, 0.74, 80);   // 1st excited: 740 keV
    TLine *l_2nd  = new TLine(3.10, 0, 3.10, 80);   // 2nd excited: 3103 keV
    TLine *l_3rd  = new TLine(4.78, 0, 4.78, 80);   // 3rd excited: 4780 keV
    TLine *l_4th  = new TLine(6.84, 0, 6.84, 80);   // 4th excited: 6841 keV

    // Estilo
    for(auto l : {l_GS, l_1st, l_2nd, l_3rd, l_4th}) {
       l->SetLineColor(kRed);
       l->SetLineWidth(5);
       l->SetLineStyle(2);  // discontinua
       l->Draw("same");
    }
 TLine *l_GS_fit  = new TLine(fModel->GetParameter(1),  0, fModel->GetParameter(1),  80);
 TLine *l_1st_fit = new TLine(fModel->GetParameter(4),  0, fModel->GetParameter(4),  80);
 TLine *l_2nd_fit = new TLine(fModel->GetParameter(7),  0, fModel->GetParameter(7),  80);
 TLine *l_3rd_fit = new TLine(fModel->GetParameter(10), 0, fModel->GetParameter(10), 80);
 TLine *l_4th_fit = new TLine(fModel->GetParameter(13), 0, fModel->GetParameter(13), 80);

 for(auto l1 : {l_GS_fit, l_1st_fit, l_2nd_fit, l_3rd_fit, l_4th_fit}) {
     l1->SetLineColor(kYellow+2);
     l1->SetLineWidth(5);
     l1->Draw("same");
 }*/

   cout << "\n";
   cout << "Starting fits for " << hHex.size() << " histograms.\n";

   // Número de parámetros que quieres guardar por fit
   const int nParams = 15;

   // Contenedor: un vector de hHex.size() vectores, cada uno con nParams
   /*std::vector<std::vector<double>> all_fit_params(hHex.size(), std::vector<double>(nParams, 0.0));

   for (int i = 0; i < hHex.size(); i++) {

      TSpectrum *sp = new TSpectrum(5); // 5 maxima search
      int nfound = sp->Search(hexCorr, 2, "", 0.05); // 2 = sigma of smoothing, last = threshold
      Double_t *xpeaks = sp->GetPositionX();

      // copy into a vector<double>
      std::vector<double> sorted_peaks(xpeaks, xpeaks + nfound);
      // sort them
      std::sort(sorted_peaks.begin(), sorted_peaks.end());
      for(int i=0;i<nfound; ++i){
         cout << sorted_peaks[i] << endl;
      }
       fExSpectra_vec[i] = new TF1(Form("fExSpectra%d", i+1),
           "gaus(0) + gaus(3) + [6]*TMath::BreitWigner(x,[7],[8]) + [9]*TMath::BreitWigner(x,[10],[11]) +
   [12]*TMath::BreitWigner(x,[13],[14])", -1, 9.);

           for (int p = 0; p < nParams; ++p) {
               fExSpectra_vec[i]->SetParameter(p, globalParamsFinals[p]);
            }
            fExSpectra_vec[i]->SetParLimits(11, 0.7, 1.9);
           // Set width constraints BEFORE fitting

      if (hHex[i]->GetEntries() < 100) {
           fExSpectra_vec[i] = nullptr; // Optional: mark as skipped
           continue; // Skip fitting this histogram
       }
       auto *c_temp = new TCanvas(Form("c_temp_%d", i), Form("Fit for hHex%d", i+1), 800, 600);

      hHex[i]->Fit(fExSpectra_vec[i]);

      int status = hHex[i]->Fit(fExSpectra_vec[i], "R");

      // Guardar parámetros del fit en el vector correspondiente
       for (int p = 0; p < nParams; ++p) {
           all_fit_params[i][p] = fExSpectra_vec[i]->GetParameter(p);
       }
   }

   std::cout << hHex.size() << " histograms processed.\n";

   // Print fit parameters for each angular bin
   auto *c_hex_segmented1 = new TCanvas("c_hex_segmented1", "Hex Spectra", 1200, 800);
   c_hex_segmented1->Divide(3, 2); //(3,3)

   for (int i = 0; i < 5; ++i) { //9
       if (i >= hHex.size()) continue; // Safety check

       double thetaMin = angularBins[i].first;
       double thetaMax = angularBins[i].second;

       c_hex_segmented1->cd(i + 1); // Switch to pad (pads are 1-indexed)
       hHex[i]->Sumw2();
       hHex[i]->Draw(); //E1

      // Gaussian 1
      TF1 *gaus1 = new TF1("gaus1", "gaus(0)", -5, 14);
      gaus1->SetParameters(globalParamsFinals[0], globalParamsFinals[1], globalParamsFinals[2]);
      gaus1->SetLineColor(kViolet);
      gaus1->Draw("same");

      // Gaussian 2
      TF1 *gaus2 = new TF1("gaus2", "gaus(0)", -5, 14);
      gaus2->SetParameters(globalParamsFinals[3], globalParamsFinals[4], globalParamsFinals[5]);
      gaus2->SetLineColor(kBlue);
      gaus2->Draw("same");

      // Breit-Wigner 1
      TF1 *bw1 = new TF1("bw1", "[0]*TMath::BreitWigner(x,[1],[2])", -5, 14);
      bw1->SetParameters(globalParamsFinals[6], globalParamsFinals[7], globalParamsFinals[8]);
      bw1->SetLineColor(kGreen+2);
      bw1->Draw("same");

      // Breit-Wigner 2
      TF1 *bw2 = new TF1("bw2", "[0]*TMath::BreitWigner(x,[1],[2])", -5, 14);
      bw2->SetParameters(globalParamsFinals[9], globalParamsFinals[10], globalParamsFinals[11]);
      bw2->SetLineColor(kMagenta);
      bw2->Draw("same");

      // Breit-Wigner 3
      TF1 *bw3 = new TF1("bw3", "[0]*TMath::BreitWigner(x,[1],[2])", -5, 14);
      bw3->SetParameters(globalParamsFinals[12], globalParamsFinals[13], globalParamsFinals[14]);
      bw3->SetLineColor(kOrange+7);
      bw3->Draw("same l");

      // Phase Space
      h_PS_1n->Scale(globalParamsFinals[15]); // Escala el histograma con el parámetro del fit
      h_PS_1n->SetLineColor(kBlack);
      h_PS_1n->SetLineWidth(2);
      h_PS_1n->Draw("same");

      TLegend* legend2 = new TLegend(0.7, 0.15, 0.9, 0.3); // (x1, y1, x2, y2) en coordenadas del canvas
      //legend2->SetBorderSize(0); // sin borde
      legend2->SetFillStyle(0);  // fondo transparente
      legend2->AddEntry(hexCorr, "hexCorr", "l");
      legend2->AddEntry(gaus1, "gaus(0)", "l");
      legend2->AddEntry(gaus2, "gaus(1)", "l");
      legend2->AddEntry(bw1, "Breit-Wigner 1", "l");
      legend2->AddEntry(bw2, "Breit-Wigner 2", "l");
      legend2->AddEntry(bw3, "Breit-Wigner 3", "l");
      legend2->AddEntry(h_PS_1n, "Phase Space Background", "l");
      legend2->Draw("same");


       //if (fExSpectra_vec[i]) {
         //  fExSpectra_vec[i]->SetLineColor(kRed); // Optional: distinguish fit
         //  fExSpectra_vec[i]->Draw("same");
       //}

       //gPad->Update(); // Force update of current pad
   }

   c_hex_segmented1->Update(); // ← actualiza todo el canvas
   */

   //---------------- Save plots ----------------//
   std::string nombre_pdf = "plots_C16_pd_C15.pdf";

   if (guardar_en_pdf) {
      // c_ExEner->Print((nombre_pdf + "(").c_str()); // abre el PDF multipágina
      // c_AngEner_Corr->Print(nombre_pdf.c_str());
      // c_AngEner->Print(nombre_pdf.c_str());
      // c_redchi2->Print(nombre_pdf.c_str());
      // c_ExenerCorr->Print(nombre_pdf.c_str());
      // c_AngDistr->Print(nombre_pdf.c_str());
      // c_ExvsZpos->Print(nombre_pdf.c_str());
      // c_hex_segmented1->Print(nombre_pdf.c_str());
      // kin->Print((nombre_pdf + ")").c_str()); // cierra el PDF multipágina

      // gSystem->Exec(("xdg-open " + nombre_pdf).c_str()); // abre el PDF automáticamente
   }

   /*fout->Write();
   fout->Close();*/

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
