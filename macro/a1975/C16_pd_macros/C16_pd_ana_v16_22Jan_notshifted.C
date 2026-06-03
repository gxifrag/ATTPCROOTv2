#include <fstream>
#include <iostream>

double Ebin_max = 10.0;
double Ebin_min = -2.0;
int NumberBins = 140; // 140
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

double TotalVerticalError(TH1F *hExp, TGraphErrors *gTheory, double scale)
{
   double errorSum = 0.0;

   for (int i = 1; i <= hExp->GetNbinsX(); ++i) {
      double x = hExp->GetBinCenter(i);
      double y_exp = hExp->GetBinContent(i);
      double y_theory = gTheory->Eval(x);

      double delta = std::abs(scale * y_exp - y_theory);
      errorSum += delta;
   }

   return errorSum;
}

Bool_t compareEventName(std::string &getname, std::string &fribname)
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
   fBW->SetNpx(1000);
   return fBW;
}

//---------------------------main function---------------------------------------
void C16_pd_ana_v16_22Jan_notshifted()
{
   bool guardar_en_pdf = false; // ← cambia a false si quieres solo verlos en pantalla
   gROOT->ProcessLine(".X /home/georgina/fair_install/ATTPCROOTv2/macro/a1975/myStyle.C");

   // Activar modo batch si estás guardando en PDF
   if (guardar_en_pdf) {
      gROOT->SetBatch(kTRUE); // ← esto evita que se abran ventanas
   } else {
      gROOT->SetBatch(kFALSE); // ← esto permite ver los canvas en pantalla
   }

   // FairRunAna *run = new FairRunAna();

   TH2F *Ang_Ener_Corr = new TH2F("Ang_Ener_Corr", "Ang_Ener_Corr", 720, 10, 60, 1000, 0, 60.0); // 14.0
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
   auto *hexCorr2 = new TH1F("hexCorr2", "", NumberBins, Ebin_min, Ebin_max);
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

   TH1F *h_PS_1n_plot = new TH1F("h_PS_1n_plot", "h_PS_1n_plot", NumberBins, Ebin_min, Ebin_max);

   /*auto *hredchi2 = new TH1F("redchi2", "redchi2", 1000, 0, 0.0001);
   auto *hbredchi2 = new TH1F("bredchi2", "bredchi2", 1000, 0, 5);*/

   auto *hexvstheta = new TH2F("hexVStheta", "hexVStheta", 100, -2, 10, 100, 0, 50);

   Double_t nc_tot[200];
   Double_t nc_PS_1n[200];
   Double_t x[200];
   Int_t nbins;

   // Some useful transformation constants.
   Double_t u_to_MeV = 931.49401;
   Double_t Brho_to_p = 1.602176634E-19;

   // Some masses that may be useful for the experiment.
   Double_t m_p = 1.007825 * u_to_MeV;
   Double_t m_d = 2.0135532 * u_to_MeV;
   Double_t m_t = 3.016049281 * u_to_MeV;
   Double_t m_He3 = 3.016029 * u_to_MeV;
   Double_t m_a = 4.00260325415 * u_to_MeV;

   Double_t m_C12 = 12.00 * u_to_MeV;
   Double_t m_C13 = 13.00335484 * u_to_MeV;
   Double_t m_C14 = 14.003242 * u_to_MeV;
   Double_t m_C15 = 15.0105993 * u_to_MeV;
   Double_t m_C16 = 16.0147 * u_to_MeV;
   Double_t m_C17 = 17.0226 * u_to_MeV;

   // Beam and target parameters.
   Double_t Ebeam_buff = 11.5 * 16; // 11.5
   cout << " Beam energy in buffer gas : " << Ebeam_buff << "\n";
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

   /*std::vector<std::pair<int, int>> angularBins = {
    {20, 30},
    {30, 40},
    {40, 50},
    {50, 70},
    {70, 90}
   };

  for (size_t i = 0; i < angularBins.size(); ++i) {
    int thetaMin = angularBins[i].first;
    int thetaMax = angularBins[i].second;
    TString histTitle = Form("Excitation Energy (%.1d deg - %.1d deg)", thetaMin, thetaMax);
    TString histName = Form("hex_%d_%d", thetaMin, thetaMax);
    //std::cout << "Procesando ángulos entre " << thetaMin << " y " << thetaMax << " grados.\n";
    //std::cout << NumberBins << " " << Ebin_min << " " << Ebin_max << "\n";

      hHex[i] = new TH1F(histName, histTitle, NumberBins, Ebin_min, Ebin_max); //90, -5, 14
   }*/

   // ELoss tables.
   // AtTools::AtELossTable *elossTableH2 = new AtTools::AtELossTable();
   // elossTableH2->LoadSrimTable("StoppingPower_SRIM_C16_H2.txt"); //SRIM no me va.
   // elossTableH2->LoadLiseTable("StoppingPower_C16_H2.txt", 2.0158,3.3084e-5);

   double densityH2 = 3.3084e-5; // g/cm³
   AtTools::AtELossCATIMA elossH2(densityH2);
   double mass{16.0147};                        // Mass of C16 in u
   elossH2.SetMaterial(catima::Material(1, 1)); // Set material to H2
   elossH2.SetProjectile(16, 6, mass);          // Set projectile to proton

   double kethe = 13.;

   filenames.push_back("run_0104_2H.root");
   filenames.push_back("run_0105_2H.root");
   filenames.push_back("run_0106_2H.root");
   filenames.push_back("run_0107_2H.root");
   filenames.push_back("run_0108_2H.root");
   filenames.push_back("run_0109_2H.root");
   filenames.push_back("run_0110_2H.root");
   // filenames.push_back("run_0111_2H.root");
   filenames.push_back("run_0112_2H.root");
   filenames.push_back("run_0113_2H.root");
   filenames.push_back("run_0114_2H.root");
   filenames.push_back("run_0115_2H.root");
   filenames.push_back("run_0116_2H.root");
   filenames.push_back("run_0117_2H.root");
   filenames.push_back("run_0118_2H.root");
   filenames.push_back("run_0119_2H.root");
   filenames.push_back("run_0120_2H.root");
   // filenames.push_back("run_0121_2H.root");
   filenames.push_back("run_0122_2H.root");
   filenames.push_back("run_0123_2H.root");
   filenames.push_back("run_0124_2H.root");
   filenames.push_back("run_0125_2H.root");
   filenames.push_back("run_0126_2H.root");
   filenames.push_back("run_0127_2H.root");
   filenames.push_back("run_0128_2H.root");
   filenames.push_back("run_0129_2H.root");
   filenames.push_back("run_0130_2H.root");
   filenames.push_back("run_0131_2H.root");
   filenames.push_back("run_0132_2H.root");
   filenames.push_back("run_0133_2H.root");
   filenames.push_back("run_0134_2H.root");
   filenames.push_back("run_0135_2H.root");
   filenames.push_back("run_0136_2H.root");
   filenames.push_back("run_0137_2H.root");
   filenames.push_back("run_0138_2H.root");
   filenames.push_back("run_0139_2H.root");
   filenames.push_back("run_0140_2H.root");
   filenames.push_back("run_0141_2H.root");
   filenames.push_back("run_0142_2H.root");
   filenames.push_back("run_0143_2H.root");
   filenames.push_back("run_0144_2H.root");
   filenames.push_back("run_0145_2H.root");
   filenames.push_back("run_0146_2H.root");
   filenames.push_back("run_0147_2H.root");
   // filenames.push_back("run_0148_2H.root");
   // filenames.push_back("run_0149_2H.root");
   filenames.push_back("run_0150_2H.root");
   filenames.push_back("run_0151_2H.root");
   filenames.push_back("run_0152_2H.root");
   filenames.push_back("run_0153_2H.root");
   filenames.push_back("run_0154_2H.root");
   filenames.push_back("run_0155_2H.root");
   filenames.push_back("run_0156_2H.root");
   filenames.push_back("run_0157_2H.root");
   filenames.push_back("run_0158_2H.root");
   filenames.push_back("run_0159_2H.root");
   filenames.push_back("run_0160_2H.root");
   filenames.push_back("run_0161_2H.root");
   filenames.push_back("run_0162_2H.root");
   filenames.push_back("run_0163_2H.root");
   filenames.push_back("run_0164_2H.root");
   filenames.push_back("run_0165_2H.root");
   filenames.push_back("run_0166_2H.root");
   filenames.push_back("run_0167_2H.root");
   filenames.push_back("run_0168_2H.root");
   filenames.push_back("run_0169_2H.root");
   filenames.push_back("run_0170_2H.root");
   filenames.push_back("run_0171_2H.root");
   filenames.push_back("run_0172_2H.root");
   filenames.push_back("run_0173_2H.root");
   filenames.push_back("run_0174_2H.root");
   filenames.push_back("run_0175_2H.root");
   filenames.push_back("run_0176_2H.root");
   filenames.push_back("run_0177_2H.root");
   filenames.push_back("run_0178_2H.root");
   filenames.push_back("run_0179_2H.root");
   filenames.push_back("run_0180_2H.root");
   filenames.push_back("run_0181_2H.root");
   filenames.push_back("run_0182_2H.root");
   filenames.push_back("run_0183_2H.root");
   filenames.push_back("run_0184_2H.root");
   filenames.push_back("run_0185_2H.root");
   filenames.push_back("run_0186_2H.root");
   filenames.push_back("run_0187_2H.root");
   filenames.push_back("run_0188_2H.root");
   filenames.push_back("run_0189_2H.root");

   for (auto filename : filenames) {
      TFile *runFile = new TFile("/home/georgina/C16_analysis/C16_H2/C16_pd_v16_root/" + filename, "R");
      // new TFile("/home/georgina/C16_analysis/C16_H2/C16_pd_v16_tb440_MMG30/InterpSolver_root/" + filename);
      TTree *Tphysics = (TTree *)runFile->Get("parquettree");

      Double_t theta{};
      Double_t phi{};
      Double_t Brho{};
      Double_t redchi{};
      Double_t zPos{};
      Double_t ke{};
      Tphysics->SetBranchAddress("polar", &theta);
      Tphysics->SetBranchAddress("azimuthal", &phi);
      Tphysics->SetBranchAddress("brho", &Brho);
      Tphysics->SetBranchAddress("redchisq", &redchi);
      Tphysics->SetBranchAddress("vertex_z", &zPos);
      Tphysics->SetBranchAddress("ke", &ke);

      /*std::string estimation_filename = std::string(filename.Data());
      size_t pos = estimation_filename.find("_2H");
      if (pos != std::string::npos) {
         estimation_filename.erase(pos, 3); // elimina "_2H"
      }

      TFile *estimationFile = new TFile(
      //("/home/georgina/C16_analysis/C16_H2/Estimation_C16_H2_v16/root/" + estimation_filename).c_str(), "R");
      ("/home/georgina/C16_analysis/C16_H2/Estimation_C16_H2_v16_tb440_MMG30/" + estimation_filename).c_str(), "R");

      TTree *Testimation = (TTree *)estimationFile->Get("parquettree");

      Double_t arclength{};
      Testimation->SetBranchAddress("arclength", &arclength);*/

      for (int i = 0; i < Tphysics->GetEntries(); i++) {
         Tphysics->GetEntry(i);
         // Testimation->GetEntry(i);
         // eventID = i;

         Double_t p_ej = Brho * Z_ej * 2.99792458 / 10 * 1000;
         Double_t E_ej = TMath::Sqrt(p_ej * p_ej + m_ej * m_ej) - m_ej;

         auto [ex_energy, theta_cm] = kine_2b(m_C16, m_p, m_b, m_B, Ebeam_buff, theta, ke);

         Double_t Ebeam_at_z = elossH2.GetEnergy(Ebeam_buff, zPos * 100); //

         // Corrección cinemática

         // double theta_lab_corr=(theta-(2.0*TMath::Pi()/4000) * (E_ej - kethe)); //theta: rad; theta_lab_corr: rad;
         // E_ej-kethe: MeV 29.5
         double theta_lab_corr = theta; // theta: rad; theta_lab_corr: rad; E_ej-kethe: MeV 29.5

         auto [ex_energy_corr, theta_cm_corr] = kine_2b(m_C16, m_p, m_b, m_B, Ebeam_at_z, theta_lab_corr, ke);

         // Fill uncorrected histogram
         hex->Fill(ex_energy);
         // cout << "Excitation Energy uncorrected: " << ex_energy << " MeV\n";
         ExvsZpos->Fill(ex_energy, zPos * 100.0);
         KineticEnergy->Fill(ke);
         Ang_Ener_Corr->Fill(theta_lab_corr * TMath::RadToDeg(),
                             ke); // theta lab!! -> I still have to implement the correction of catima?

         double Ebeam_buff_test = 0.0;
         double ex_energy_test = 0.0;

         for (int j = 0; j < 60; j++) {
            Ebeam_buff_test = 5. * j; // MeV
            auto [ex_energy_test, theta_cm_test] = kine_2b(m_C16, m_p, m_b, m_B, Ebeam_buff_test, theta_lab_corr, ke);
            Ebeam_test->Fill(ex_energy_test, Ebeam_buff_test);
            // cout << ex_energy_test << " " << "ebeam_buffer" << Ebeam_buff_test<< "\n";
         }

         // Fill corrected histogram
         if (zPos * 100 > 2.0 && zPos * 100 < 60.0 && ke < 25.0) {
            ExCorrvsZpos->Fill(ex_energy_corr, zPos * 100.0);
            hexCorr->Fill(ex_energy_corr);  //-0.5 MeV shift to match the known C15 gs energy
            hexCorr2->Fill(ex_energy_corr); //-0.5 MeV shift to match the known C15 gs energy
         }

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

         /*for (size_t i = 0; i < angularBins.size(); ++i) {
          int thetaMin = angularBins[i].first;
          int thetaMax = angularBins[i].second;

             if (theta_cm > thetaMin && theta_cm <= thetaMax && ke <25.) {
                hHex[i]->Fill(ex_energy_corr-0.55);

                break; // Only fill one bin per event
             }
          }*/

         // tEvents->Fill();
      } // events
   } // Files

   AngDistrCM->Divide(new TF1("sin", "sin(x * TMath::DegToRad())", 0, 180));

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

   TH1F *h_PS_1n = new TH1F("h_PS_1n", "h_PS_1n", NumberBins, Ebin_min, Ebin_max);

   for (int i = 0; i < treePS->GetEntries(); i++) {
      treePS->GetEntry(i);
      if (ThetaCM_cal > ThetaCM_min && ThetaCM_cal < ThetaCM_max) {
         h_PS_1n->Fill(Ex_cal, Weight_sim);
      }
   }
   h_PS_1n->Smooth();

   TGraph *graphPS = histoToTgraph(h_PS_1n);

   // -----------------------------KINEMATICS FOR DIFFERENT EXCITATION ENERGIES

   std::vector<std::string> files = {"C16_pd_C15_gs_Ebeam11_5.txt", "C16_pd_C15_740keV_Ebeam11_5.txt",
                                     "C16_pd_C15_3103keV_Ebeam11_5.txt", "C16_pd_C15_4780keV_Ebeam11_5.txt",
                                     "C16_pd_C15_6841keV_Ebeam11_5.txt"};
   std::vector<std::string> labels = {"Ground State", "1st Excited State (740keV)", "2nd Excited State (3103keV)",
                                      "3rd Excited State (4780keV)", "4th Excited State(6841keV)"};

   // Colors for each line (ROOT color codes: 2=red,4=blue,8=green, etc.)
   std::vector<int> colors = {kOrange + 7, kBlue, kGreen + 2, kMagenta, kRed + 2};
   std::vector<TGraph *> graphs;

   for (size_t i = 0; i < files.size(); i++) {
      TString fileKine =
         Form("/home/georgina/fair_install/ATTPCROOTv2/macro/Kinematics/Decay_kinematics/%s", files[i].c_str());
      std::ifstream kineStr(fileKine.Data());

      if (kineStr.fail()) {
         std::cout << " Warning : No Kinematics file found for " << labels[i] << "!" << std::endl;
         continue;
      }

      // Temporary storage
      std::vector<Double_t> ThetaCMS, ThetaLabRec, EnerLabRec, ThetaLabSca, EnerLabSca;

      Double_t tCMS, tLabRec, eLabRec, tLabSca, eLabSca;
      while (kineStr >> tCMS >> tLabRec >> eLabRec >> tLabSca >> eLabSca) {
         ThetaCMS.push_back(tCMS);
         ThetaLabRec.push_back(tLabRec);
         EnerLabRec.push_back(eLabRec);
         ThetaLabSca.push_back(tLabSca);
         EnerLabSca.push_back(eLabSca);
      }

      // Build graph
      /*TGraph *g = new TGraph(ThetaLabRec.size(), ThetaLabRec.data(), EnerLabRec.data());
      g->SetLineColor(colors[i]);
      g->SetLineWidth(2);
      g->SetTitle(labels[i].c_str());

      graphs.push_back(g);*/
   }

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
      cout << sorted_peaks[i] << endl;
   }

   // Define a two-gaussian TF1 (ROOT built-in gaus uses amplitude = height)
   TF1 *f2g = new TF1("f2g", "gaus(0) + gaus(3)", -1., 1.4);
   double m1 = sorted_peaks[0];
   double m2 = sorted_peaks[1];

   // Peak 1
   f2g->SetParameter(1, m1);                                           // mean of gaus(0)
   f2g->SetParameter(2, 0.2);                                          // sigma guess
   f2g->SetParameter(0, hexCorr->GetBinContent(hexCorr->FindBin(m1))); // amplitude guess

   // Peak 2
   f2g->SetParameter(4, m2);                                           // mean of gaus(3)
   f2g->SetParameter(5, 0.2);                                          // sigma guess
   f2g->SetParameter(3, hexCorr->GetBinContent(hexCorr->FindBin(m2))); // amplitude guess

   hexCorr->Fit(f2g, "R"); // R = use the range you specified, 0 no plot, M minuit
   f2g->SetLineColor(kViolet + 2);

   // hexCorr->Draw();
   f2g->Draw("same"); // <--- REQUIRED so the fit curve is drawn

   // Breit-Wigner 1
   TF1 *bwprefit1 = new TF1("bw1", "[0]*TMath::BreitWigner(x,[1],[2])", 3.3 - 0.55, 4.3 - 0.55);
   TF1 *bwprefit2 = new TF1("bw2", "[0]*TMath::BreitWigner(x,[1],[2])", 4.5 - 0.55, 6.25 - 0.55);
   TF1 *bwprefit3 = new TF1("bw3", "[0]*TMath::BreitWigner(x,[1],[2])", 6., 7.5); // 6., 8.

   double m3 = sorted_peaks[3];
   double m4 = sorted_peaks[4];
   double m5 = sorted_peaks[6];

   bwprefit1->SetParameters(230, m3, 0.5);
   // bwprefit1->SetParameters(230,3.5-0.55,0.5);
   hexCorr->Fit(bwprefit1, "R0");
   bwprefit1->SetLineColor(kRed);
   bwprefit1->Draw("same");

   bwprefit2->SetParameters(150, m4, 2.);
   // bwprefit2->SetParameters(150,5.2-0.55,2.);
   hexCorr->Fit(bwprefit2, "R0+");
   bwprefit2->SetLineColor(kBlue);
   bwprefit2->Draw("same");

   bwprefit3->SetParameters(60, m5, 0.5); // 0.5
   // bwprefit3->SetParameters(60,6.8-0.55,0.5);//0.5
   hexCorr->Fit(bwprefit3, "R0+");
   bwprefit3->SetLineColor(kGreen);
   bwprefit3->Draw("same");

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
   TF1 *fModel = new TF1("fModel", model, -1., 9., 16, // number fitted parameters
                         "SpectralModelSeg");

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
                                          0.0001};
   // Assigns all the parameters at the same time
   for (size_t i = 0; i < globalParamsIni.size(); ++i) {
      fModel->SetParameter(i, globalParamsIni[i]);
      // cout << "i= " << i << "globalParamsIni= " << globalParamsIni[i] << endl;
   }

   fModel->SetParLimits(11, 0.7 - 0.55, 1.9 - 0.55);
   // fModel->SetParLimits(14,6.25,7.5);

   if (!hexCorr2) {
      std::cerr << "hexCorr2 is null\n";
      return;
   }
   if (!h_PS_1n) {
      std::cerr << "h_PS_1n is null\n";
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
   hexCorr2->Draw("E1"); // E1

   TLine *vline0 = new TLine(1.218, 0, 1.218, 380);
   vline0->SetLineColor(kRed); // opcional
   vline0->SetLineStyle(2);    // opcional: línea discontinua
   vline0->SetLineWidth(3);    // opcional
   vline0->Draw("SAME");

   // Gaussian 1
   TF1 *gaus1 = new TF1("gaus1", "gaus(0)", Ebin_min, Ebin_max);
   gaus1->SetParameters(globalParamsFinals[0], globalParamsFinals[1], globalParamsFinals[2]);
   gaus1->SetLineColor(kOrange + 7);
   gaus1->Draw("same");

   // Gaussian 2
   TF1 *gaus2 = new TF1("gaus2", "gaus(0)", Ebin_min, Ebin_max);
   gaus2->SetParameters(globalParamsFinals[3], globalParamsFinals[4], globalParamsFinals[5]);
   gaus2->SetLineColor(kBlue);
   gaus2->Draw("same");

   // Breit-Wigner 1
   TF1 *bw1 = new TF1("bw1", "[0]*TMath::BreitWigner(x,[1],[2])", Ebin_min, Ebin_max);
   bw1->SetParameters(globalParamsFinals[6], globalParamsFinals[7], globalParamsFinals[8]);
   bw1->SetLineColor(kGreen + 2);
   bw1->Draw("same");

   // Breit-Wigner 2
   TF1 *bw2 = new TF1("bw2", "[0]*TMath::BreitWigner(x,[1],[2])", Ebin_min, Ebin_max);
   bw2->SetParameters(globalParamsFinals[9], globalParamsFinals[10], globalParamsFinals[11]);
   bw2->SetLineColor(kMagenta);
   bw2->Draw("same");

   // Breit-Wigner 3
   TF1 *bw3 = new TF1("bw3", "[0]*TMath::BreitWigner(x,[1],[2])", Ebin_min, Ebin_max);
   bw3->SetParameters(globalParamsFinals[12], globalParamsFinals[13], globalParamsFinals[14]);
   bw3->SetLineColor(kRed + 2);
   bw3->Draw("same");

   // Phase Space
   if (h_PS_1n->GetNbinsX() > 0) {
      h_PS_1n->Scale(globalParamsFinals[15]);
   }
   h_PS_1n->SetLineColor(kGray + 2);
   h_PS_1n->Draw("same");

   // Suppose you fitted with 'fitFcn' (could be gaus1, bw1, etc.)
   double chi2 = fModel->GetChisquare();
   int ndf = fModel->GetNDF();
   double chi2Ndf = chi2 / ndf;

   TLatex latex;
   latex.SetNDC();          // normalized coordinates
   latex.SetTextSize(0.03); // smaller than legend text
   latex.DrawLatex(0.2, 0.85, Form("#chi^{2}/NDF = %.2f", chi2Ndf));

   TLegend *legend2 = new TLegend(0.6, 0.6, 0.9, 0.9); // (x1, y1, x2, y2) en coordenadas del canvas
   legend2->AddEntry(hexCorr2, "Spectrum", "l");
   legend2->AddEntry(gaus1, "Ground State", "l");
   legend2->AddEntry(gaus2, "1st Excited State", "l");
   legend2->AddEntry(bw1, "2nd Excited State", "l");
   legend2->AddEntry(bw2, "3rd Excited State", "l");
   legend2->AddEntry(bw3, "4th Excited State", "l");
   legend2->AddEntry(h_PS_1n, "Phase Space Bkg", "l");
   legend2->Draw("same");

   c_ExEner->Update();

   //---------------- Plots ----------------//

   TCanvas *c_AngEner = new TCanvas("AngEner", "Energy as a function of #theta", 800, 1200);
   c_AngEner->cd();
   Ang_Ener_Corr->Draw("col");
   Ang_Ener_Corr->GetXaxis()->SetTitle("#theta_lab (deg)");
   Ang_Ener_Corr->GetYaxis()->SetTitle("Kinetic Energy (MeV)");

   double xMin = Ang_Ener_Corr->GetXaxis()->GetXmin();
   double xMax = Ang_Ener_Corr->GetXaxis()->GetXmax();

   // Create and draw a horizontal line at y=25
   TLine *line = new TLine(xMin, 25, xMax, 25);
   line->SetLineColor(kRed); // optional: set line color
   line->SetLineStyle(2);    // optional: dashed line
   line->SetLineWidth(2);    // optional: thicker line
   line->Draw("SAME");

   // Draw all kinematics graphs from the loop
   for (size_t i = 0; i < graphs.size(); i++) {
      graphs[i]->Draw("L SAME"); // "L SAME" draws as a line on the same canvas
   }

   // Optional: add a legend
   auto legend = new TLegend(0.35, 0.7, 0.9, 0.9);
   for (size_t i = 0; i < graphs.size(); i++) {
      legend->AddEntry(graphs[i], labels[i].c_str(), "l");
   }
   legend->Draw();

   TCanvas *c_check2 = new TCanvas("check2", "check2", 1200, 800);
   c_check2->cd();
   hexvstheta->GetXaxis()->SetTitle("Excitation Energy (MeV)");
   hexvstheta->GetYaxis()->SetTitle("#theta_{lab} (deg)");
   hexvstheta->Draw("colz");

   TCanvas *test = new TCanvas("test", "test", 1200, 800);
   test->cd();
   Ebeam_test->GetXaxis()->SetTitle("Excitation Energy (MeV)");
   Ebeam_test->GetYaxis()->SetTitle("E beam (MeV)");
   Ebeam_test->Draw("colz");

   // Excitation energy spectrum with fits --------------------------------------

   /*TCanvas* hGS_energy = new TCanvas("hGS_ExEnergy", "Ground State Excitation Energy", 1000, 600);
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
   auto *g2 = new TGraphErrors("/home/georgina/twofnr/21.2ndExcitedP", "%lg %lg");
   double bestScale = 1.0;
   double minError = 1e9;

   for (double testScale = 0.01; testScale <= 5.0; testScale += 0.01) {
      double error = TotalVerticalError(hGS_AngularDistr, g2, testScale);
      if (error < minError) {
         minError = error;
         bestScale = testScale;
      }
   }
   std::cout << "Mejor escala por distancia vertical: " << bestScale << std::endl;

   TCanvas *cAng3103_overlay = new TCanvas("cAng3103_overlay", "c", 1000, 600);
   cAng3103_overlay->cd();

   cAng3103_overlay->SetLogy();
   // Histograma experimental primero
   h3103_AngularDistr->Sumw2();
   h3103_AngularDistr->SetLineColor(kRed);
   h3103_AngularDistr->SetLineWidth(2);
   h3103_AngularDistr->SetMarkerStyle(20);   // marcador redondo sólido
   h3103_AngularDistr->SetMarkerColor(kRed); // color del marcador
   h3103_AngularDistr->SetMarkerSize(1.2);   // tamaño del marcador

   h3103_AngularDistr->GetXaxis()->SetTitle("#theta_CM (deg)");
   h3103_AngularDistr->GetYaxis()->SetTitle("d#sigma/d#Omega (mb/sr)");
   h3103_AngularDistr->Scale(0.05); // ← se aplica el ajuste automático
   h3103_AngularDistr->Draw("E1");

   g2->SetLineColor(kRed);
   g2->SetStats(0);
   g2->SetLineWidth(2);
   g2->SetMarkerStyle(20);
   g2->GetYaxis()->SetRangeUser(0.1, 1000);
   g2->SetMarkerColor(kBlack);
   g2->Draw("same P");

   auto *legend9 = new TLegend(0.75, 0.75, 0.88, 0.88);
   legend9->AddEntry(h740_AngularDistr, "2nd: p_{1/2}^{-}", "p"); // use h740_AngularDistr here
   legend9->SetTextSize(0.04);
   legend9->AddEntry(g2, "twofnr", "lp");
   legend9->Draw();

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

   TCanvas *c_ExvsZpos = new TCanvas("ExCorrvsZpos", "Excitation Energy vs z position and track length", 1200, 800);
   c_ExvsZpos->cd();
   ExCorrvsZpos->GetXaxis()->SetTitle("Excitation Energy (MeV)");
   ExCorrvsZpos->GetYaxis()->SetTitle("z (cm)");

   gPad->SetRightMargin(0.20);
   ExCorrvsZpos->Draw("zcol");

   cout << "\n";
   cout << "Starting fits for " << hHex.size() << " histograms.\n";

   // Número de parámetros que quieres guardar por fit
   const int nParams = 15;

   // Contenedor: un vector de hHex.size() vectores, cada uno con nParams
   std::vector<std::vector<double>> all_fit_params(hHex.size(), std::vector<double>(nParams, 0.0));

   /*for (int i = 0; i < hHex.size(); i++) {

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
               fExSpectra_vec[i]->SetParameter(p, globalParams[p]);
            }
            fExSpectra_vec[i]->SetParLimits(11, 0.7-0.55, 1.9-0.55);
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
       //hHex[i]->Sumw2();
       hHex[i]->Draw(); //E1

      // Gaussian 1
      TF1 *gaus1 = new TF1("gaus1", "gaus(0)", -5, 14);
      gaus1->SetParameters(globalParams[0], globalParams[1], globalParams[2]);
      gaus1->SetLineColor(kViolet);
      gaus1->Draw("same");

      // Gaussian 2
      TF1 *gaus2 = new TF1("gaus2", "gaus(0)", -5, 14);
      gaus2->SetParameters(globalParams[3], globalParams[4], globalParams[5]);
      gaus2->SetLineColor(kBlue);
      gaus2->Draw("same");

      // Breit-Wigner 1
      TF1 *bw1 = new TF1("bw1", "[0]*TMath::BreitWigner(x,[1],[2])", -5, 14);
      bw1->SetParameters(globalParams[6], globalParams[7], globalParams[8]);
      bw1->SetLineColor(kGreen+2);
      bw1->Draw("same");

      // Breit-Wigner 2
      TF1 *bw2 = new TF1("bw2", "[0]*TMath::BreitWigner(x,[1],[2])", -5, 14);
      bw2->SetParameters(globalParams[9], globalParams[10], globalParams[11]);
      bw2->SetLineColor(kMagenta);
      bw2->Draw("same");

      // Breit-Wigner 3
      TF1 *bw3 = new TF1("bw3", "[0]*TMath::BreitWigner(x,[1],[2])", -5, 14);
      bw3->SetParameters(globalParams[12], globalParams[13], globalParams[14]);
      bw3->SetLineColor(kOrange+7);
      bw3->Draw("same l");

      // Phase Space
      h_PS_1n->Scale(globalParams[15]); // Escala el histograma con el parámetro del fit
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
   }*/

   // c_hex_segmented1->Update(); // ← actualiza todo el canvas

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
