#include <fstream>
#include <iostream>

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

   // for inverse kinematics Note: this angle corresponds to the recoil
   theta_cm = TMath::Pi() - acos((pow(s, 2) + s * (2 * t - pow(m1, 2) - pow(m2, 2) - pow(m3, 2) - pow(m4_ex, 2)) +
                                  (pow(m1, 2) - pow(m2, 2)) * (pow(m3, 2) - pow(m4_ex, 2))) /
                                 (omega(s, pow(m1, 2), pow(m2, 2)) * omega(s, pow(m3, 2), pow(m4_ex, 2))));

   theta_cm = theta_cm * TMath::RadToDeg();
   return std::make_tuple(Ex, theta_cm);
}

void GetEnergy(Double_t M, Double_t IZ, Double_t BRO, Double_t &E);

TF1 *CreateSpectralModelWithPS(const char *name, TH1F *h_PS_1n)
{
   TF1 *fModel = new TF1(
      name,
      [=](double *x, double *p) {
         double val = 0;
         val += p[0] * TMath::Gaus(x[0], p[1], p[2], true);
         val += p[3] * TMath::Gaus(x[0], p[4], p[5], true);
         val += p[6] * TMath::BreitWigner(x[0], p[7], p[8]);
         val += p[9] * TMath::BreitWigner(x[0], p[10], p[11]);
         val += p[12] * TMath::BreitWigner(x[0], p[13], p[14]);

         int bin = h_PS_1n->FindBin(x[0]);
         if (bin >= 1 && bin <= h_PS_1n->GetNbinsX())
            val += p[15] * h_PS_1n->GetBinContent(bin);

         return val;
      },
      -5, 14, 16);

   fModel->SetNpx(1000);
   return fModel;
}

void C16_pd_ana_v16_1Nov()
{
   bool guardar_en_pdf = false; // ← cambia a false si quieres solo verlos en pantalla

   // Activar modo batch si estás guardando en PDF
   if (guardar_en_pdf) {
      gROOT->SetBatch(kTRUE); // ← esto evita que se abran ventanas
   } else {
      gROOT->SetBatch(kFALSE); // ← esto permite ver los canvas en pantalla
   }

   /*TFile *fout = new TFile("event_data.root", "RECREATE");
   if (!fout || fout->IsZombie()) {
       std::cerr << "ERROR: No se pudo crear el archivo event_data.root" << std::endl;
       return;
   }
   fout->cd();

   TTree *tEvents = new TTree("eventTree", "Datos completos por evento");

   Double_t thetacm_corr, ex_corr, thetalab_corr;
   Int_t eventID;

   tEvents->Branch("thetacm_corr", &thetacm_corr, "thetacm_corr/D");
   tEvents->Branch("ex_corr", &ex_corr, "ex_corr/D");
   tEvents->Branch("thetalab_corr", &thetalab_corr, "thetalab_corr/D");
   tEvents->Branch("eventID", &eventID, "eventID/I");*/

   // FairRunAna *run = new FairRunAna();

   double Ebin_max = 10.0;
   double Ebin_min = -3.0;
   int NumberBins = 100; // 100
   int NumberBinsAux = 200;

   TH2F *Ang_Ener = new TH2F("Ang_Ener", "Ang_Ener", 720, 10, 60, 1000, 0, 60.0); // 14.0
   TH2F *Ang_Ener_Corr = new TH2F("Ang_Ener_Corr", "Ang_Ener_Corr", 720, 10, 60, 1000, 0, 60.0);

   TH2F *ELossvsBrho = new TH2F("ELossvsBrho", "ELossvsBrho", 4000, 0, 25000, 1000, 0, 4);
   TH2F *dedxvsBrho = new TH2F("dedxvsBrho", "dedxvsBrho", 4000, 0, 10000, 1000, 0, 4);
   TH2F *hVxVy = new TH2F("hVxVy", "hVxVy", 1000, 0, 4, 1000, 0, 4);
   TH1F *hGS_AngularDistr = new TH1F("hGS_AngularDistr", "Ground State Angular Distribution", 30, 0, 180);

   auto *hex = new TH1F("hex", "hex", NumberBins, Ebin_min, Ebin_max);
   auto *QvsEb = new TH2F("QvsEb", "QvsEb", 1000, -5, 15, 100, 0, 300);
   auto *QvsZpos = new TH2F("QvsZpos", "QvsZpos", 1000, -10, 50, 200, -100, 100);
   auto *hexCorr = new TH1F("hexCorr", "hexCorr", NumberBins, Ebin_min, Ebin_max);
   TH1F *hGS_ExEnergy =
      new TH1F("hGS_ExEnergy", "Excitation Energy (Ground State Cut)", NumberBins, Ebin_min, Ebin_max);

   auto *AngDistr = new TH1F("Ang_Distr", "Ang_Distr", 128, 0, 120);
   auto *AngDistrCM = new TH1F("Ang_Distr_CM", "Ang_Distr_CM", 90, 0, 180);
   auto *ExvsZpos = new TH2F("ExvsZpos", "ExvsZpos", 1000, -5, 15, 200, -20, 150);
   auto *ExvsTrackLength = new TH2F("ExvsTrackLength", "ExvsTrackLength", 1000, -5, 15, 200, -20, 150);
   auto *ExCorrvsZpos = new TH2F("ExCorrvsZpos", "ExCorrvsZpos", 1000, -10, 10, 200, -100, 100);
   auto *KineticEnergy = new TH1F("KineticEnergy", "KineticEnergy", 100, 0, 100);

   TH1F *h_PS_1n_plot = new TH1F("h_PS_1n_plot", "h_PS_1n_plot", NumberBins, Ebin_min, Ebin_max);

   /*auto *hredchi2 = new TH1F("redchi2", "redchi2", 1000, 0, 0.0001);
   auto *hbredchi2 = new TH1F("bredchi2", "bredchi2", 1000, 0, 5);*/

   auto *hexvstheta = new TH2F("hexVStheta", "hexVStheta", 1000, -5, 15, 90, 0, 90);

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

   std::vector<std::pair<int, int>> angularBins = {{20, 30}, {30, 40}, {40, 50}, {50, 70}, {70, 90}};

   for (size_t i = 0; i < angularBins.size(); ++i) {
      int thetaMin = angularBins[i].first;
      int thetaMax = angularBins[i].second;
      TString histTitle = Form("Excitation Energy (%.1d deg - %.1d deg)", thetaMin, thetaMax);
      TString histName = Form("hex_%d_%d", thetaMin, thetaMax);
      std::cout << "Procesando ángulos entre " << thetaMin << " y " << thetaMax << " grados.\n";

      hHex[i] = new TH1F(histName, histTitle, 90, -5, 14); // 90, -5, 14
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

      std::string estimation_filename = std::string(filename.Data());
      size_t pos = estimation_filename.find("_2H");
      if (pos != std::string::npos) {
         estimation_filename.erase(pos, 3); // elimina "_2H"
      }

      TFile *estimationFile = new TFile(
         ("/home/georgina/C16_analysis/C16_H2/Estimation_C16_pd_v16/root/" + estimation_filename).c_str(), "R");

      TTree *Testimation = (TTree *)estimationFile->Get("parquettree");

      Double_t arclength{};
      Testimation->SetBranchAddress("arclength", &arclength);

      for (int i = 0; i < Tphysics->GetEntries(); i++) {
         Tphysics->GetEntry(i);
         Testimation->GetEntry(i);
         // eventID = i;

         Double_t p_ej = Brho * Z_ej * 2.99792458 / 10 * 1000;
         Double_t E_ej = TMath::Sqrt(p_ej * p_ej + m_ej * m_ej) - m_ej;

         auto [ex_energy, theta_cm] = kine_2b(m_C16, m_p, m_b, m_B, Ebeam_buff, theta, ke);

         Double_t Ebeam_at_z = elossH2.GetEnergy(Ebeam_buff, zPos * 100); //

         // Correccion cinematica

         double theta_lab_corr =
            (theta - (2.0 * TMath::Pi() / 4000) * (E_ej - kethe)); // theta: rad; theta_lab_corr: rad; E_ej-kethe: MeV

         auto [ex_energy_corr, theta_cm_corr] = kine_2b(m_C16, m_p, m_b, m_B, Ebeam_at_z, theta_lab_corr, ke);

         /*thetacm_corr = theta_cm_corr;
         ex_corr = ex_energy_corr;
         thetalab_corr = theta_lab_corr;*/

         // Fill uncorrected histogram
         hex->Fill(ex_energy);
         ExvsZpos->Fill(ex_energy, zPos * 100.0);
         KineticEnergy->Fill(ke);

         // Fill corrected histogram
         if (zPos * 100 > 2.0 && zPos * 100 < 60.0) {
            // cout << " Ex corrected : " << ex_energy_corr << "  " << " zpos: " << zPos*100.0 << " Ebeam at z: " <<
            // Ebeam_at_z << "\n";
            ExCorrvsZpos->Fill(ex_energy_corr, zPos * 100.0);
            hexCorr->Fill(ex_energy_corr); //-0.5 MeV shift to match the known C15 gs energy
         }

         Double_t theta_deg_lower = 25.0;
         Double_t theta_deg_upper = 100.0;

         Double_t energyCutGS_lower = -0.3;
         Double_t energyCutGS_upper = 0.7;

         if (ex_energy_corr > energyCutGS_lower && ex_energy_corr < energyCutGS_upper &&
             theta_cm_corr > theta_deg_lower && theta_cm_corr < theta_deg_upper) {
            Double_t theta_deg = theta_cm_corr;
            hGS_AngularDistr->Fill(theta_deg);
            hGS_ExEnergy->Fill(ex_energy_corr);
            // cout << " Ex corrected (GS cut) : " << ex_energy_corr << "  " << " \n";
            // cout << " theta cm corr (GS cut) : " << theta_deg << "  " << " \n";
         }

         /*if(ex_ener_corr > 5.0 && ex_energy_corr < 7.0 &&
             theta_cm_corr > theta_deg_lower && theta_cm_corr < theta_deg_upper) {
            Double_t theta_deg = theta_cm_corr;
            //h740keV_AngularDistr->Fill(theta_deg);
         }*/

         // Histograms
         // hredchi2->Fill(redchi);

         Ang_Ener->Fill(theta * TMath::RadToDeg(),
                        ke); // theta lab!! -> I still have to implement the correction of catima?
         Ang_Ener_Corr->Fill(theta_lab_corr * TMath::RadToDeg(), ke); //

         Double_t vx = TMath::Sin(theta) * TMath::Sqrt(ke);
         Double_t vy = TMath::Cos(theta) * TMath::Sqrt(ke);

         hVxVy->Fill(vx, vy);

         AngDistr->Fill(theta * TMath::RadToDeg());
         AngDistrCM->Fill(theta_cm);
         hexvstheta->Fill(ex_energy, theta * TMath::RadToDeg());
         ExvsTrackLength->Fill(ex_energy, arclength);

         for (size_t i = 0; i < angularBins.size(); ++i) {
            int thetaMin = angularBins[i].first;
            int thetaMax = angularBins[i].second;

            if (theta_cm > thetaMin && theta_cm <= thetaMax) {
               hHex[i]->Fill(ex_energy_corr);

               break; // Only fill one bin per event
            }
         }

         // tEvents->Fill();
      } // events
   } // Files

   AngDistrCM->Divide(new TF1("sin", "sin(x * TMath::DegToRad())", 0, 180));

   //--------------------fit------------------------------------------------------

   nbins = hexCorr->GetNbinsX();
   int binmax = hexCorr->GetMaximumBin();
   double ThetaCM_min = 0;
   double ThetaCM_max = 180;

   TString PhaseSpace_FileName = "PhaseSpace_16C_pd_1n.root";

   TFile *f_PS = new TFile(PhaseSpace_FileName, "READ");
   TTree *t_PS = (TTree *)f_PS->Get("simulated_tree");

   double Weight_sim, Ex_cal, ThetaCM_cal;

   t_PS->SetBranchAddress("Weight_sim", &Weight_sim);
   t_PS->SetBranchAddress("Ex_cal", &Ex_cal);
   t_PS->SetBranchAddress("ThetaCM_cal", &ThetaCM_cal);

   std::cout << "Entries PS tree = " << t_PS->GetEntries() << std::endl;

   TH1F *h_PS_1n = new TH1F("h_PS_1n", "h_PS_1n", NumberBins, Ebin_min, Ebin_max);

   for (int i = 0; i < t_PS->GetEntries(); i++) {
      t_PS->GetEntry(i);
      if (ThetaCM_cal > ThetaCM_min && ThetaCM_cal < ThetaCM_max) {
         h_PS_1n->Fill(Ex_cal, Weight_sim);
      }
      if (i < 10) {
         std::cout << "Weight = " << Weight_sim << std::endl;
      }
   }
   h_PS_1n->Smooth();

   std::cout << "---- PhaseSpace file ----" << std::endl;

   if (!f_PS || f_PS->IsZombie()) {
      std::cerr << "ERROR: no se abre el fichero PS\n";
      return;
   }

   f_PS->ls(); // 🔥 CLAVE: ver qué hay dentro

   if (!t_PS) {
      std::cerr << "ERROR: tree no encontrado\n";
      return;
   }

   std::cout << "Entries en PS tree = " << t_PS->GetEntries() << std::endl;

   // -----------------------------KINEMATICS FOR DIFFERENT EXCITATION ENERGIES

   std::vector<std::string> files = {"C16_pd_C15_gs_Ebeam11_5.txt", "C16_pd_C15_740keV_Ebeam11_5.txt",
                                     "C16_pd_C15_3103keV_Ebeam11_5.txt", "C16_pd_C15_4780keV_Ebeam11_5.txt",
                                     "C16_pd_C15_6841keV_Ebeam11_5.txt"};
   std::vector<std::string> labels = {"gs", "740keV", "3103keV", "4780keV", "6841keV"};

   // Colors for each line (ROOT color codes: 2=red,4=blue,8=green, etc.)
   std::vector<int> colors = {kBlack, kRed, kBlue, kGreen + 2, kMagenta};
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
      TGraph *g = new TGraph(ThetaLabRec.size(), ThetaLabRec.data(), EnerLabRec.data());
      g->SetLineColor(colors[i]);
      g->SetLineWidth(2);
      g->SetTitle(labels[i].c_str());

      graphs.push_back(g);
   }

   //---------------- Fitting the experimental data ----------------//
   ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
   TFile *filePS = new TFile("PhaseSpace_16C_pd_1n.root", "READ");
   TTree *treePS = (TTree *)filePS->Get("simulated_tree");

   // Búsqueda de picos con TSpectrum
   auto *cps = new TCanvas("cps", "cps", 800, 600);
   cps->cd();
   int nPeaksToSearch = 5;
   TSpectrum *spec = new TSpectrum(nPeaksToSearch);
   int nFound = spec->Search(hexCorr, 2, "", 0.1);
   delete cps;

   // Vector para almacenar los parámetros de los pre-fits
   std::vector<std::vector<double>> prefit_params;

   // Para cada pico encontrado, hacer un fit gaussiano local
   for (int i = 0; i < nFound; i++) {
      double peak_pos = spec->GetPositionX()[i];
      double peak_height = hexCorr->GetBinContent(hexCorr->FindBin(peak_pos));

      // Crear una gaussiana para el pre-fit
      TF1 *gaus_prefit = new TF1(Form("gaus_prefit_%d", i), "gaus", peak_pos - 0.5, peak_pos + 0.5);

      // Parámetros iniciales para la gaussiana
      gaus_prefit->SetParameters(peak_height, peak_pos, 0.2);

      // Hacer el fit en un rango limitado alrededor del pico
      auto *ctest0 = new TCanvas("ctest0", "ctest0", 800, 600);
      hexCorr->Fit(gaus_prefit, "RQN+", "", peak_pos - 0.5, peak_pos + 0.5);
      delete ctest0;

      // Almacenar los parámetros del fit
      std::vector<double> params = {
         gaus_prefit->GetParameter(0), // amplitud
         gaus_prefit->GetParameter(1), // media
         gaus_prefit->GetParameter(2)  // sigma
      };
      prefit_params.push_back(params);
      delete gaus_prefit; // Limpieza de memoria
   }

   // First sort prefit_params by mean value (index 1)
   std::sort(prefit_params.begin(), prefit_params.end(), [](const auto &a, const auto &b) { return a[1] < b[1]; });

   double initParams[16]; // Declarar el array fuera del condicional

   // Ordenar los picos por posición
   std::sort(prefit_params.begin(), prefit_params.end(), [](const auto &a, const auto &b) { return a[1] < b[1]; });

   // Verificar el número de picos encontrados e inicializar initParams
   // Usar los parámetros del pre-fit
   initParams[0] = prefit_params[0][0]; // Primera gaussiana
   initParams[1] = prefit_params[0][1];
   initParams[2] = prefit_params[0][2];
   initParams[3] = prefit_params[1][0]; // Segunda gaussiana
   initParams[4] = prefit_params[1][1];
   initParams[5] = prefit_params[1][2];
   initParams[6] = prefit_params[2][0]; // Primer BW
   initParams[7] = prefit_params[2][1];
   initParams[8] = prefit_params[2][2]; // 0.2
   initParams[9] = prefit_params[3][0]; // Segundo BW
   initParams[10] = prefit_params[3][1];
   initParams[11] = prefit_params[3][2];
   initParams[12] = prefit_params[4][0]; // Tercer BW
   initParams[13] = prefit_params[4][1];
   initParams[14] = prefit_params[4][2];
   initParams[15] = 0.0005; // Phase space 0.0005

   /*double initParams[16] = {400, 0.57,0.218,
       500, 1.31, 0.218,
       220, 3.9, 0.2,
       50, 5.5, 0.2,
       25, 6.9, 0.2, 0.0005}; //BW: amplitude mean width 0.0001
   */
   auto *ctest = new TCanvas("ctest", "ctest", 800, 600);
   TF1 *model = CreateSpectralModelWithPS("fExSpectra", h_PS_1n);

   /*model->SetParLimits(2, -0.3, 0.3);    // Límites para anchura de Gaussian 1
   model->SetParLimits(10, 5.0, 6.);   // Límites para mean de BW 2
   model->SetParLimits(11, 1.5, 2.3);   // Límites para width de BW 2
   model->SetParLimits(13, 6.5, 7.2);*/   // Límites para mean de BW 3
   /// model->SetParLimits(14, 0.1, 0.5);   // Límites para width de BW 3
   model->SetParLimits(2, -0.1, 0.25);

   // Parámetros iniciales (ajústalos según tu caso)

   model->SetParameters(initParams);             // fondo PS
   hexCorr->Fit("fExSpectra", "", "", 0.4, 8.5); // "R" para rango, "Q" para modo silencioso
   delete ctest;

   // Extraer parámetros ajustados
   std::vector<double> finalParams;
   for (int i = 0; i < model->GetNpar(); ++i) {
      finalParams.push_back(model->GetParameter(i));
      // cout << "Fitted p[" << i << "] = " << model->GetParameter(i) << std::endl;
   }

   std::ofstream out("fit_params_16Cpd15C.txt");
   for (size_t i = 0; i < finalParams.size(); ++i) {
      out << finalParams[i] << "\n";
   }
   out.close();

   //---------------- Plots ----------------//

   /*TCanvas *c_AngEner = new TCanvas("AngEner", "Energy as a function of #theta", 1200, 800);
   Ang_Ener->SetMarkerStyle(20);
   Ang_Ener->SetMarkerSize(0.5);
   Ang_Ener->Draw("col");
   Ang_Ener->GetXaxis()->SetTitle("#theta_lab (deg)");
   Ang_Ener->GetYaxis()->SetTitle("Energy (MeV)"); //Energy (MeV)
   // Draw all kinematics graphs from the loop

   for (size_t i = 0; i < graphs.size(); i++) {
       graphs[i]->Draw("L SAME");  // "L SAME" draws as a line on the same canvas
   }

   // Optional: add a legend
   auto legend = new TLegend(0.65, 0.65, 0.88, 0.88);
   for (size_t i = 0; i < graphs.size(); i++) {
       legend->AddEntry(graphs[i], labels[i].c_str(), "l");
   }
   legend->Draw();

   TCanvas *c_AngEner_Corr = new TCanvas("AngEnerCorr", "Energy as a function of #theta_Corr", 1200, 800);
   Ang_Ener_Corr->SetMarkerStyle(20);
   Ang_Ener_Corr->SetMarkerSize(0.5);
   Ang_Ener_Corr->Draw("col");
   Ang_Ener_Corr->GetXaxis()->SetTitle("#theta_lab (deg)");
   Ang_Ener_Corr->GetYaxis()->SetTitle("Energy (MeV)"); //Energy (MeV)

   // Draw all kinematics graphs from the loop
   for (size_t i = 0; i < graphs.size(); i++) {
       graphs[i]->Draw("L SAME");  // "L SAME" draws as a line on the same canvas
   }

   // Optional: add a legend
   auto legend0 = new TLegend(0.65, 0.65, 0.88, 0.88);
   for (size_t i = 0; i < graphs.size(); i++) {
       legend0->AddEntry(graphs[i], labels[i].c_str(), "l");
   }
   legend0->Draw();*/

   // Excitation energy spectrum with fits --------------------------------------

   /*TCanvas *c_hex = new TCanvas("hex", "Excited Energy spectra", 1200, 800);
      c_hex->cd();
      hex->GetXaxis()->SetTitle("Excitation Energy (MeV)");
      hex->GetYaxis()->SetTitle("Counts");
      hex->SetLineColor(kRed);
      hex->Draw("hist");
      //hexCorr->SetLineColor(kRed);
      hexCorr->Draw("hist same");*/

   TCanvas *hGS_energy = new TCanvas("hGS_ExEnergy", "Ground State Excitation Energy", 1000, 600);
   hGS_energy->cd();
   hGS_ExEnergy->GetXaxis()->SetTitle("Excitation Energy (MeV)");
   hGS_ExEnergy->GetYaxis()->SetTitle("Counts");
   hGS_ExEnergy->SetLineColor(kBlue + 2);
   hGS_ExEnergy->SetLineWidth(2);
   hGS_ExEnergy->Draw();

   /*TCanvas* cAngGS = new TCanvas("cAngGS", "Ground State Angular Distribution", 1000, 600);
   cAngGS->cd();
   //gPad->SetLogy();
   // Histograma experimental primero
   hGS_AngularDistr->Sumw2();
   hGS_AngularDistr->Scale(0.5);  // ← se aplica correctamente
   hGS_AngularDistr->SetLineColor(kBlue+2);
   hGS_AngularDistr->SetLineWidth(2);
   hGS_AngularDistr->GetXaxis()->SetTitle("Theta_cm (deg)");
   hGS_AngularDistr->GetYaxis()->SetTitle("Counts");
   hGS_AngularDistr->Draw("E1");  // ← define el marco
   */

   //-----------------------------------------------------------------------------------------
   auto *g = new TGraphErrors("/home/georgina/twofnr/21.groundState", "%lg %lg");
   double bestScale = 1.0;
   double minError = 1e9;

   for (double testScale = 0.01; testScale <= 5.0; testScale += 0.01) {
      double error = TotalVerticalError(hGS_AngularDistr, g, testScale);
      if (error < minError) {
         minError = error;
         bestScale = testScale;
      }
   }
   std::cout << "Mejor escala por distancia vertical: " << bestScale << std::endl;

   TCanvas *cAngGS_overlay = new TCanvas("cAngGS_overlay", "Ground State Angular Distribution Overlay", 1000, 600);
   cAngGS_overlay->cd();

   cAngGS_overlay->SetLogy();
   // Histograma experimental primero
   hGS_AngularDistr->Sumw2();
   hGS_AngularDistr->SetLineColor(kRed);
   hGS_AngularDistr->SetLineWidth(2);
   hGS_AngularDistr->SetMarkerStyle(20);   // marcador redondo sólido
   hGS_AngularDistr->SetMarkerColor(kRed); // color del marcador
   hGS_AngularDistr->SetMarkerSize(1.2);   // tamaño del marcador

   hGS_AngularDistr->GetXaxis()->SetTitle("Theta_cm (deg)");
   hGS_AngularDistr->GetYaxis()->SetTitle("Counts");
   hGS_AngularDistr->Scale(bestScale); // ← se aplica el ajuste automático
   hGS_AngularDistr->Draw("E1");

   g->SetLineColor(kRed);
   g->SetStats(0);
   g->SetLineWidth(2);
   g->SetMarkerStyle(20);
   g->GetYaxis()->SetRangeUser(0.1, 1000);
   g->SetMarkerColor(kBlack);
   g->Draw("same P");

   auto *legend = new TLegend(0.6, 0.7, 0.88, 0.88);
   legend->AddEntry(hGS_AngularDistr, "Angular distribution ground state", "l");
   legend->AddEntry(g, "twofnr, OMP: Chapel and Daehnick", "lp");
   legend->Draw();

   //-----------------------------------------------------------------------------------------

   TCanvas *c_ExEner = new TCanvas("ExEner", "Corrected Excited Energy spectra", 1200, 800);
   c_ExEner->cd();
   hexCorr->SetStats(0);
   // hexCorr->Sumw2(); // activa almacenamiento de errores
   // hexCorr->Fit(model, "RQ");
   hexCorr->GetXaxis()->SetTitle("Excitation Energy (MeV)");
   hexCorr->GetYaxis()->SetTitle("Counts");
   hexCorr->Draw(); // E1

   // Gaussian 1
   TF1 *gaus1 = new TF1("gaus1", "gaus(0)", -5, 14);
   gaus1->SetParameters(finalParams[0], finalParams[1], finalParams[2]);
   gaus1->SetLineColor(kViolet);
   gaus1->SetNpx(1000);
   gaus1->Draw("same");

   // Gaussian 2
   TF1 *gaus2 = new TF1("gaus2", "gaus(0)", -5, 14);
   gaus2->SetParameters(finalParams[3], finalParams[4], finalParams[5]);
   gaus2->SetLineColor(kBlue);
   gaus2->SetNpx(1000);
   gaus2->Draw("same");

   //---------------------------------------------------------------------

   cout << " Integral total fit: " << model->Integral(0, 2) / hexCorr->GetBinWidth(0) << "\n";
   cout << " Integral gaus 1: " << gaus1->Integral(-5, 14) / hexCorr->GetBinWidth(0) << "\n";
   cout << " Integral gaus 2: " << gaus2->Integral(-5, 14) / hexCorr->GetBinWidth(0) << "\n";
   cout << "hexCorr Integral: " << hexCorr->Integral(48, 66) << "\n";
   //---------------------------------------------------------------------

   // Breit-Wigner 1
   TF1 *bw1 = new TF1("bw1", "[0]*TMath::BreitWigner(x,[1],[2])", -5, 14);
   bw1->SetParameters(finalParams[6], finalParams[7], finalParams[8]);
   bw1->SetLineColor(kGreen + 2);
   bw1->SetNpx(1000);
   bw1->Draw("same");

   // Breit-Wigner 2
   TF1 *bw2 = new TF1("bw2", "[0]*TMath::BreitWigner(x,[1],[2])", -5, 14);
   bw2->SetParameters(finalParams[9], finalParams[10], finalParams[11]);
   bw2->SetLineColor(kMagenta);
   bw2->SetNpx(1000);
   bw2->Draw("same");

   // Breit-Wigner 3
   TF1 *bw3 = new TF1("bw3", "[0]*TMath::BreitWigner(x,[1],[2])", -5, 14);
   bw3->SetParameters(finalParams[12], finalParams[13], finalParams[14]);
   bw3->SetLineColor(kOrange + 7);
   bw3->SetNpx(1000);
   bw3->Draw("same l");

   // Phase Space
   h_PS_1n->Scale(finalParams[15]); // Escala el histograma con el parámetro del fit
   h_PS_1n->SetLineColor(kBlack);
   h_PS_1n->SetLineWidth(2);
   h_PS_1n->Draw("same");

   TLegend *legend2 = new TLegend(0.7, 0.7, 0.9, 0.9); // (x1, y1, x2, y2) en coordenadas del canvas
   // legend2->SetBorderSize(0); // sin borde
   legend2->SetFillStyle(0); // fondo transparente
   legend2->AddEntry(hexCorr, "hexCorr", "l");
   legend2->AddEntry(gaus1, "gaus(0)", "l");
   legend2->AddEntry(gaus2, "gaus(1)", "l");
   legend2->AddEntry(bw1, "Breit-Wigner 1", "l");
   legend2->AddEntry(bw2, "Breit-Wigner 2", "l");
   legend2->AddEntry(bw3, "Breit-Wigner 3", "l");
   legend2->AddEntry(h_PS_1n, "Phase Space Background", "l");
   legend2->Draw();

   c_ExEner->Update();

   //--------------------------------------------------------------------------------------------------

   TCanvas *c_ExenerCorr = new TCanvas("ExenerCorr", "Excited Energy spectra corrected", 1200, 800);
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
   AngDistrCM->GetYaxis()->SetTitle("#frac{d#sigma}{d#Omega} (a.u.)");

   TCanvas *c_ExvsZpos = new TCanvas("ExvsZpos", "Excitation Energy vs z position and track length", 1200, 800);
   c_ExvsZpos->cd();
   c_ExvsZpos->Divide(2, 1);
   c_ExvsZpos->cd(1);
   ExvsZpos->Draw("zcol");
   ExvsZpos->GetXaxis()->SetTitle("Excitation Energy (MeV)");
   ExvsZpos->GetYaxis()->SetTitle("z (cm)");
   c_ExvsZpos->cd(2);
   ExvsTrackLength->Draw("zcol");
   ExvsTrackLength->GetXaxis()->SetTitle("Excitation Energy (MeV)");
   ExvsTrackLength->GetYaxis()->SetTitle("Track Length (cm)");

   TCanvas *kin = new TCanvas("kin", "kin", 1200, 800);
   KineticEnergy->Sumw2();
   kin->cd();
   KineticEnergy->Draw("E1");

   /*std::vector<std::vector<double>> all_fit_params(hHex.size(), std::vector<double>(15, 0.0));

   cout << "\n";
   cout << "Starting fits for " << hHex.size() << " histograms.\n";

   for (int i = 0; i < hHex.size(); i++) {
       fExSpectra_vec[i] = new TF1(Form("fExSpectra%d", i+1),
           "gaus(0) + gaus(3) + [6]*TMath::BreitWigner(x,[7],[8]) + [9]*TMath::BreitWigner(x,[10],[11]) +
   [12]*TMath::BreitWigner(x,[13],[14])", -5, 14); fExSpectra_vec[i]->SetParameters(initParams);
           // Set width constraints BEFORE fitting

      if (hHex[i]->GetEntries() < 100) {
           fExSpectra_vec[i] = nullptr; // Optional: mark as skipped
           continue; // Skip fitting this histogram
       }
       auto *c_temp = new TCanvas(Form("c_temp_%d", i), Form("Fit for hHex%d", i+1), 800, 600);

         hHex[i]->Fit(fExSpectra_vec[i]);

         int status = hHex[i]->Fit(fExSpectra_vec[i], "RQ");

         delete c_temp;

         for (int p = 0; p < 15; ++p) {
            all_fit_params[i][p] = fExSpectra_vec[i]->GetParameter(p);

         }
   }

   std::cout << hHex.size() << " histograms processed.\n";*/

   // Print fit parameters for each angular bin
   /*auto *c_hex_segmented1 = new TCanvas("c_hex_segmented1", "Hex Spectra", 1200, 800);
   c_hex_segmented1->Divide(3, 2); //(3,3)

   for (int i = 0; i < 5; ++i) { //9
       if (i >= hHex.size()) continue; // Safety check

       double thetaMin = angularBins[i].first;
       double thetaMax = angularBins[i].second;

       c_hex_segmented1->cd(i + 1); // Switch to pad (pads are 1-indexed)
       hHex[i]->Sumw2();
       hHex[i]->Draw("E1");

      // Gaussian 1
      TF1 *gaus1 = new TF1("gaus1", "gaus(0)", -5, 14);
      gaus1->SetParameters(finalParams[0], finalParams[1], finalParams[2]);
      gaus1->SetLineColor(kViolet);
      gaus1->Draw("same");

      // Gaussian 2
      TF1 *gaus2 = new TF1("gaus2", "gaus(0)", -5, 14);
      gaus2->SetParameters(finalParams[3], finalParams[4], finalParams[5]);
      gaus2->SetLineColor(kBlue);
      gaus2->Draw("same");

      // Breit-Wigner 1
      TF1 *bw1 = new TF1("bw1", "[0]*TMath::BreitWigner(x,[1],[2])", -5, 14);
      bw1->SetParameters(finalParams[6], finalParams[7], finalParams[8]);
      bw1->SetLineColor(kGreen+2);
      bw1->Draw("same");

      // Breit-Wigner 2
      TF1 *bw2 = new TF1("bw2", "[0]*TMath::BreitWigner(x,[1],[2])", -5, 14);
      bw2->SetParameters(finalParams[9], finalParams[10], finalParams[11]);
      bw2->SetLineColor(kMagenta);
      bw2->Draw("same");

      // Breit-Wigner 3
      TF1 *bw3 = new TF1("bw3", "[0]*TMath::BreitWigner(x,[1],[2])", -5, 14);
      bw3->SetParameters(finalParams[12], finalParams[13], finalParams[14]);
      bw3->SetLineColor(kOrange+7);
      bw3->Draw("same l");

      // Phase Space
      h_PS_1n->Scale(finalParams[15]); // Escala el histograma con el parámetro del fit
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