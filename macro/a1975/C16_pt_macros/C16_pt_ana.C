#include <fstream>
#include <iostream>

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

void C16_pt_ana()
{
   // FairRunAna *run = new FairRunAna();

   TH2F *Ang_Ener = new TH2F("Ang_Ener", "Ang_Ener", 720, 0, 179, 1000, 0, 20.0);
   TH2F *Ang_Ener_PRAC = new TH2F("Ang_Ener_PRAC", "Ang_Ener_PRAC", 1000, 0, 100, 1000, 0, 200.0);
   TH2F *ELossvsBrho = new TH2F("ELossvsBrho", "ELossvsBrho", 4000, 0, 25000, 1000, 0, 4);
   TH2F *dedxvsBrho = new TH2F("dedxvsBrho", "dedxvsBrho", 4000, 0, 10000, 1000, 0, 4);
   TH2F *hVxVy = new TH2F("hVxVy", "hVxVy", 1000, 0, 4, 1000, 0, 4);
   TH1F *henergyIC = new TH1F("henergyIC", "henergyIC", 2048, 0, 2047);

   auto *hex = new TH1F("hex", "hex", 80, -2, 15);
   auto *QvsEb = new TH2F("QvsEb", "QvsEb", 1000, -5, 15, 100, 0, 300);
   auto *QvsZpos = new TH2F("QvsZpos", "QvsZpos", 1000, -10, 50, 200, -100, 100);
   auto *hexCorr = new TH1F("hexCorr", "hexCorr", 80, -2, 15);

   auto *AngDistr = new TH1F("Ang_Distr", "Ang_Distr", 2048, 0, 120);
   auto *AngDistrCM = new TH1F("Ang_Distr_CM", "Ang_Distr_CM", 360, 0, 180);
   auto *ExvsZpos = new TH2F("ExvsZpos", "ExvsZpos", 1000, -5, 15, 200, -20, 150);
   auto *ExvsTrackLength = new TH2F("ExvsTrackLength", "ExvsTrackLength", 1000, -5, 15, 200, -20, 150);
   auto *ExCorrvsZpos = new TH2F("ExCorrvsZpos", "ExCorrvsZpos", 1000, -10, 10, 200, -100, 100);

   auto *hredchi2 = new TH1F("redchi2", "redchi2", 1000, 0, 0.0001);
   auto *hbredchi2 = new TH1F("bredchi2", "bredchi2", 1000, 0, 5);
   auto *hex11 = new TH1F("hex11", "hex11", 100, -5, 20);
   auto *hex12 = new TH1F("hex12", "hex12", 100, -5, 20);
   auto *hex13 = new TH1F("hex13", "hex13", 100, -5, 20);
   auto *hex21 = new TH1F("hex21", "hex21", 100, -5, 20);
   auto *hex22 = new TH1F("hex22", "hex22", 100, -5, 20);
   auto *hex23 = new TH1F("hex23", "hex23", 100, -5, 20);
   auto *hex31 = new TH1F("hex31", "hex31", 100, -5, 20);
   auto *hex32 = new TH1F("hex32", "hex32", 100, -5, 20);
   auto *hex33 = new TH1F("hex33", "hex33", 100, -5, 20);

   auto *hex41 = new TH1F("hex41", "hex41", 100, -5, 20);
   auto *hex42 = new TH1F("hex42", "hex42", 100, -5, 20);
   auto *hex43 = new TH1F("hex43", "hex43", 100, -5, 20);
   auto *hex51 = new TH1F("hex51", "hex54", 100, -5, 20);
   auto *hex52 = new TH1F("hex52", "hex52", 100, -5, 20);
   auto *hex53 = new TH1F("hex53", "hex53", 100, -5, 20);
   auto *hex61 = new TH1F("hex61", "hex61", 100, -5, 20);
   auto *hex62 = new TH1F("hex62", "hex62", 100, -5, 20);
   auto *hex63 = new TH1F("hex63", "hex63", 100, -5, 20);

   auto *hexvstheta = new TH2F("hexVStheta", "hexVStheta", 1000, -5, 15, 90, 0, 90);

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
   // For the 16C(p,t)15C reaction: A(a,b)B
   // B refers to the residual nucleus (14C), and b refers to the outgoing ejectile (triton, t).
   Double_t Ebeam_buff = 10.5 * 16; // 11.5 MeV per nucleon for C16 beam
   Double_t m_b = m_t;              // mass of triton (ejectile)
   Double_t m_B = m_C14;            // mass of 14C (residual nucleus)

   // Ejectile parameters.
   int A_ej = 3;
   int Z_ej = 1;
   Double_t m_ej = m_t;

   std::vector<TString> filenames;

   filenames.push_back("run_0110_3H.root");
   // filenames.push_back("run_0111_3H.root");
   filenames.push_back("run_0112_3H.root");
   filenames.push_back("run_0113_3H.root");
   filenames.push_back("run_0114_3H.root");
   filenames.push_back("run_0115_3H.root");
   filenames.push_back("run_0116_3H.root");
   filenames.push_back("run_0117_3H.root");
   filenames.push_back("run_0118_3H.root");
   filenames.push_back("run_0119_3H.root");
   filenames.push_back("run_0120_3H.root");
   filenames.push_back("run_0122_3H.root");
   filenames.push_back("run_0123_3H.root");
   filenames.push_back("run_0124_3H.root");
   filenames.push_back("run_0125_3H.root");
   filenames.push_back("run_0126_3H.root");
   filenames.push_back("run_0127_3H.root");
   filenames.push_back("run_0128_3H.root");
   filenames.push_back("run_0129_3H.root");
   filenames.push_back("run_0130_3H.root");
   filenames.push_back("run_0131_3H.root");
   filenames.push_back("run_0132_3H.root");
   filenames.push_back("run_0133_3H.root");
   filenames.push_back("run_0134_3H.root");
   filenames.push_back("run_0135_3H.root");
   filenames.push_back("run_0136_3H.root");
   filenames.push_back("run_0137_3H.root");
   filenames.push_back("run_0138_3H.root");
   filenames.push_back("run_0139_3H.root");
   filenames.push_back("run_0140_3H.root");
   filenames.push_back("run_0141_3H.root");
   filenames.push_back("run_0142_3H.root");
   filenames.push_back("run_0143_3H.root");
   // filenames.push_back("run_0144_3H.root");
   filenames.push_back("run_0145_3H.root");
   filenames.push_back("run_0146_3H.root");
   filenames.push_back("run_0147_3H.root");
   // filenames.push_back("run_0148_3H.root");
   filenames.push_back("run_0149_3H.root");
   filenames.push_back("run_0150_3H.root");
   filenames.push_back("run_0151_3H.root");
   filenames.push_back("run_0152_3H.root");
   filenames.push_back("run_0153_3H.root");
   filenames.push_back("run_0154_3H.root");
   filenames.push_back("run_0155_3H.root");
   filenames.push_back("run_0156_3H.root");
   filenames.push_back("run_0157_3H.root");
   filenames.push_back("run_0158_3H.root");
   filenames.push_back("run_0159_3H.root");
   filenames.push_back("run_0160_3H.root");
   filenames.push_back("run_0161_3H.root");
   filenames.push_back("run_0162_3H.root");
   filenames.push_back("run_0163_3H.root");
   filenames.push_back("run_0164_3H.root");
   filenames.push_back("run_0165_3H.root");
   filenames.push_back("run_0166_3H.root");
   filenames.push_back("run_0167_3H.root");
   filenames.push_back("run_0168_3H.root");
   filenames.push_back("run_0169_3H.root");
   filenames.push_back("run_0170_3H.root");
   filenames.push_back("run_0171_3H.root");
   filenames.push_back("run_0172_3H.root");
   filenames.push_back("run_0173_3H.root");
   filenames.push_back("run_0174_3H.root");
   filenames.push_back("run_0175_3H.root");
   filenames.push_back("run_0176_3H.root");
   filenames.push_back("run_0177_3H.root");
   // filenames.push_back("run_0178_3H.root");
   filenames.push_back("run_0179_3H.root");
   filenames.push_back("run_0180_3H.root");
   filenames.push_back("run_0181_3H.root");
   filenames.push_back("run_0182_3H.root");
   filenames.push_back("run_0183_3H.root");
   filenames.push_back("run_0184_3H.root");
   filenames.push_back("run_0185_3H.root");
   filenames.push_back("run_0186_3H.root");
   filenames.push_back("run_0187_3H.root");
   filenames.push_back("run_0188_3H.root");
   filenames.push_back("run_0189_3H.root");

   for (auto filename : filenames) {
      TFile *runFile =
         new TFile("/home/georgina/C16_analysis/C16_H2/C16_pt/InterpolationSolver_pt_root/" + filename, "R");
      TTree *Tphysics = (TTree *)runFile->Get("parquettree");

      Double_t theta{};
      Double_t phi{};
      Double_t Brho{};
      Double_t redchi{};
      Double_t zPos{};
      Tphysics->SetBranchAddress("polar", &theta);
      Tphysics->SetBranchAddress("azimuthal", &phi);
      Tphysics->SetBranchAddress("brho", &Brho);
      Tphysics->SetBranchAddress("redchisq", &redchi);
      Tphysics->SetBranchAddress("vertex_z", &zPos);

      for (int i = 0; i < Tphysics->GetEntries(); i++) {
         Tphysics->GetEntry(i);

         Double_t p_ej = Brho * Z_ej * 2.99792458 / 10 * 1000;

         Double_t E_ej = TMath::Sqrt(p_ej * p_ej + m_ej * m_ej) - m_ej;

         auto [ex_energy, theta_cm] = kine_2b(m_C16, m_p, m_b, m_B, Ebeam_buff, theta, E_ej);

         /*if(E_ej < 0.0 || E_ej > 15.0)
             continue;*/

         // if(theta > 90 || theta < 10.0)
         // continue;

         // if(redchi > 0.000004)
         //  continue;

         /*if(theta > 90 || theta < 10.0)
             continue;

         if(energy*Am < 0.0 || energy*Am > 20.0)
             continue;*/

         if (zPos * 100.0 < 1.0)
            continue;

         // if(ex_energy < -1.5 || ex_energy > 1.5)
         // continue;

         Double_t p0 = 0.618451;
         Double_t p1 = -0.00409149;
         Double_t mFactor = 1.00;
         Double_t offSet = 0.0;
         Double_t QcorrZ = 0.0;
         QcorrZ = ex_energy - mFactor * p1 * (zPos * 100.0) - p0;
         hexCorr->Fill(QcorrZ);

         // Histograms
         // hredchi2->Fill(redchi);
         Ang_Ener->Fill(theta * TMath::RadToDeg(), E_ej);
         hex->Fill(ex_energy);

         Double_t vx = TMath::Sin(theta) * TMath::Sqrt(E_ej);
         Double_t vy = TMath::Cos(theta) * TMath::Sqrt(E_ej);

         hVxVy->Fill(vx, vy);

         AngDistr->Fill(theta * TMath::RadToDeg());
         AngDistrCM->Fill(theta_cm);

         ExvsZpos->Fill(ex_energy, zPos * 100.0);
         ExCorrvsZpos->Fill(QcorrZ, zPos * 100.0);

         if (theta_cm > 10 && theta_cm <= 12.5)
            hex11->Fill(ex_energy);
         if (theta_cm > 12.5 && theta_cm <= 15)
            hex12->Fill(ex_energy);
         if (theta_cm > 15 && theta_cm <= 17.5)
            hex13->Fill(ex_energy);
         if (theta_cm > 17.5 && theta_cm <= 20)
            hex21->Fill(ex_energy);
         if (theta_cm > 20 && theta_cm <= 22.5)
            hex22->Fill(ex_energy);
         if (theta_cm > 22.5 && theta_cm <= 25)
            hex23->Fill(ex_energy);
         if (theta_cm > 25 && theta_cm <= 27.5)
            hex31->Fill(ex_energy);
         if (theta_cm > 27.5 && theta_cm <= 30)
            hex32->Fill(ex_energy);
         if (theta_cm > 30 && theta_cm <= 32.5)
            hex33->Fill(ex_energy);
         if (theta_cm > 32.5 && theta_cm <= 35)
            hex41->Fill(ex_energy);
         if (theta_cm > 35 && theta_cm <= 37.5)
            hex42->Fill(ex_energy);
         if (theta_cm > 37.5 && theta_cm <= 40)
            hex43->Fill(ex_energy);
         if (theta_cm > 40 && theta_cm <= 42.5)
            hex51->Fill(ex_energy);
         if (theta_cm > 42.5 && theta_cm <= 45)
            hex52->Fill(ex_energy);
         if (theta_cm > 45 && theta_cm <= 47.5)
            hex53->Fill(ex_energy);
         if (theta_cm > 47.5 && theta_cm <= 50)
            hex61->Fill(ex_energy);
         if (theta_cm > 50 && theta_cm <= 52.5)
            hex62->Fill(ex_energy);
         if (theta_cm > 52.5 && theta_cm <= 55)
            hex63->Fill(ex_energy);

         hexvstheta->Fill(ex_energy, theta * TMath::RadToDeg());

      } // events
   }    // Files

   AngDistrCM->Divide(new TF1("sin", "sin(x * TMath::DegToRad())", 0, 180));

   // Kinematics
   Double_t *ThetaCMS = new Double_t[20000];
   Double_t *ThetaLabRec = new Double_t[20000];
   Double_t *EnerLabRec = new Double_t[20000];
   Double_t *ThetaLabSca = new Double_t[20000];
   Double_t *EnerLabSca = new Double_t[20000];
   Double_t *MomLabRec = new Double_t[20000];

   TString fileKine =
      "/home/georgina/fair_install/ATTPCROOTv2/macro/Kinematics/Decay_kinematics/C16_pt_14C_gs_Ebeam11_5.txt";
   std::ifstream *kineStr = new std::ifstream(fileKine.Data());
   Int_t numKin = 0;

   if (!kineStr->fail()) {
      while (!kineStr->eof()) {
         *kineStr >> ThetaCMS[numKin] >> ThetaLabRec[numKin] >> EnerLabRec[numKin] >> ThetaLabSca[numKin] >>
            EnerLabSca[numKin];
         // numKin++;

         // MomLabRec[numKin] =( pow(EnerLabRec[numKin] + M_Ener,2) - TMath::Power(M_Ener, 2))/1000.0;
         // std::cout<<" Momentum : " <<MomLabRec[numKin]<<"\n";
         // Double_t E = TMath::Sqrt(TMath::Power(p, 2) + TMath::Power(M_Ener, 2)) - M_Ener;
         numKin++;
      }
   } else if (kineStr->fail())
      std::cout << " Warning : No Kinematics file found for this reaction GS!" << std::endl;

   TGraph *kine_gs = new TGraph(numKin, ThetaLabRec, EnerLabRec);

   // Kinematics
   Double_t *ThetaCMS1 = new Double_t[20000];
   Double_t *ThetaLabRec1 = new Double_t[20000];
   Double_t *EnerLabRec1 = new Double_t[20000];
   Double_t *ThetaLabSca1 = new Double_t[20000];
   Double_t *EnerLabSca1 = new Double_t[20000];
   Double_t *MomLabRec1 = new Double_t[20000];

   TString fileKine1 =
      "/home/georgina/fair_install/ATTPCROOTv2/macro/Kinematics/Decay_kinematics/C16_pt_14C_1st_Ebeam11_5.txt";
   std::ifstream *kineStr1 = new std::ifstream(fileKine1.Data());
   Int_t numKin1 = 0;

   if (!kineStr1->fail()) {
      while (!kineStr1->eof()) {
         *kineStr1 >> ThetaCMS1[numKin1] >> ThetaLabRec1[numKin1] >> EnerLabRec1[numKin1] >> ThetaLabSca1[numKin1] >>
            EnerLabSca1[numKin1];
         // numKin++;

         // MomLabRec[numKin] =( pow(EnerLabRec[numKin] + M_Ener,2) - TMath::Power(M_Ener, 2))/1000.0;
         // std::cout<<" Momentum : " <<MomLabRec[numKin]<<"\n";
         // Double_t E = TMath::Sqrt(TMath::Power(p, 2) + TMath::Power(M_Ener, 2)) - M_Ener;
         numKin1++;
      }
   } else if (kineStr1->fail())
      std::cout << " Warning : No Kinematics file found for this reaction 1!" << std::endl;

   TGraph *kine_1st = new TGraph(numKin1, ThetaLabRec1, EnerLabRec1);
   kine_1st->SetLineColor(kRed);
   kine_1st->SetLineWidth(2);
   kine_1st->SetLineStyle(2); // dashed line

   Double_t *ThetaCMS2 = new Double_t[20000];
   Double_t *ThetaLabRec2 = new Double_t[20000];
   Double_t *EnerLabRec2 = new Double_t[20000];
   Double_t *ThetaLabSca2 = new Double_t[20000];
   Double_t *EnerLabSca2 = new Double_t[20000];
   Double_t *MomLabRec2 = new Double_t[20000];

   TString fileKine2 =
      "/home/georgina/fair_install/ATTPCROOTv2/macro/Kinematics/Decay_kinematics/C16_pt_14C_2nd_Ebeam11_5.txt";

   std::ifstream *kineStr2 = new std::ifstream(fileKine2.Data());
   Int_t numKin2 = 0;

   if (!kineStr2->fail()) {
      while (!kineStr2->eof()) {
         *kineStr2 >> ThetaCMS2[numKin2] >> ThetaLabRec2[numKin2] >> EnerLabRec2[numKin2] >> ThetaLabSca2[numKin2] >>
            EnerLabSca2[numKin2];
         // numKin++;

         // MomLabRec[numKin] =( pow(EnerLabRec[numKin] + M_Ener,2) - TMath::Power(M_Ener, 2))/1000.0;
         // std::cout<<" Momentum : " <<MomLabRec[numKin]<<"\n";
         // Double_t E = TMath::Sqrt(TMath::Power(p, 2) + TMath::Power(M_Ener, 2)) - M_Ener;
         numKin2++;
      }
   } else if (kineStr2->fail())
      std::cout << " Warning : No Kinematics file found for this reaction 2!" << std::endl;

   TGraph *kine_2nd = new TGraph(numKin2, ThetaLabRec2, EnerLabRec2);
   kine_2nd->SetLineColor(kBlue);
   kine_2nd->SetLineWidth(2);
   kine_2nd->SetLineStyle(2); // dashed line

   /*Double_t *ThetaCMS3 = new Double_t[20000];
   Double_t *ThetaLabRec3 = new Double_t[20000];
   Double_t *EnerLabRec3 = new Double_t[20000];
   Double_t *ThetaLabSca3 = new Double_t[20000];
   Double_t *EnerLabSca3 = new Double_t[20000];
   Double_t *MomLabRec3 = new Double_t[20000];

   TString fileKine3 =
   "/home/georgina/fair_install/ATTPCROOTv2/macro/Kinematics/Decay_kinematics/C15_dt_3rd_Ebeam11_5.txt";

   std::ifstream *kineStr3 = new std::ifstream(fileKine3.Data());
   Int_t numKin3 = 0;

   if (!kineStr3->fail()) {
      while (!kineStr3->eof()) {
         *kineStr3 >> ThetaCMS3[numKin3] >> ThetaLabRec3[numKin3] >> EnerLabRec3[numKin3] >> ThetaLabSca3[numKin3] >>
            EnerLabSca3[numKin3];
         // numKin++;

         // MomLabRec[numKin] =( pow(EnerLabRec[numKin] + M_Ener,2) - TMath::Power(M_Ener, 2))/1000.0;
         // std::cout<<" Momentum : " <<MomLabRec[numKin]<<"\n";
         // Double_t E = TMath::Sqrt(TMath::Power(p, 2) + TMath::Power(M_Ener, 2)) - M_Ener;
         numKin3++;
      }
   } else if (kineStr3->fail())
      std::cout << " Warning : No Kinematics file found for this reaction!" << std::endl;

   TGraph *kine_2nd = new TGraph(numKin3, ThetaLabRec3, EnerLabRec3);
   kine_3rd->SetLineColor(kOrange);
   kine_3rd->SetLineWidth(2);
   kine_3rd->SetLineStyle(2); // dashed line
   */

   Double_t *ThetaCMS4 = new Double_t[20000];
   Double_t *ThetaLabRec4 = new Double_t[20000];
   Double_t *EnerLabRec4 = new Double_t[20000];
   Double_t *ThetaLabSca4 = new Double_t[20000];
   Double_t *EnerLabSca4 = new Double_t[20000];
   Double_t *MomLabRec4 = new Double_t[20000];

   TString fileKine4 =
      "/home/georgina/fair_install/ATTPCROOTv2/macro/Kinematics/Decay_kinematics/C16_pt_14C_3rd_Ebeam11_5.txt";

   std::ifstream *kineStr4 = new std::ifstream(fileKine4.Data());
   Int_t numKin4 = 0;

   if (!kineStr4->fail()) {
      while (!kineStr4->eof()) {
         *kineStr4 >> ThetaCMS4[numKin4] >> ThetaLabRec4[numKin4] >> EnerLabRec4[numKin4] >> ThetaLabSca4[numKin4] >>
            EnerLabSca4[numKin4];
         // numKin++;

         // MomLabRec[numKin] =( pow(EnerLabRec[numKin] + M_Ener,2) - TMath::Power(M_Ener, 2))/1000.0;
         // std::cout<<" Momentum : " <<MomLabRec[numKin]<<"\n";
         // Double_t E = TMath::Sqrt(TMath::Power(p, 2) + TMath::Power(M_Ener, 2)) - M_Ener;
         numKin4++;
      }
   } else if (kineStr4->fail())
      std::cout << " Warning : No Kinematics file found for this reaction 3!" << std::endl;

   TGraph *kine_3rd = new TGraph(numKin4, ThetaLabRec4, EnerLabRec4);
   kine_3rd->SetLineColor(kGreen);
   kine_3rd->SetLineWidth(2);
   kine_3rd->SetLineStyle(2); // dashed line

   // End Kinematics

   TCanvas *c1 = new TCanvas();
   c1->Divide(2, 2);
   c1->Draw();
   c1->cd(1);

   if (Ang_Ener->GetEntries() > 0) {
      Ang_Ener->SetMarkerStyle(20);
      Ang_Ener->SetMarkerSize(0.5);
      Ang_Ener->Draw("col");
      Ang_Ener->GetXaxis()->SetTitle("Angle (deg)");
      Ang_Ener->GetYaxis()->SetTitle("Energy (MeV)");
      kine_gs->Draw("SAME");
      kine_1st->Draw("SAME");
      kine_2nd->Draw("SAME");
      kine_3rd->Draw("SAME");
   }
   c1->cd(2);
   if (Ang_Ener_PRAC->GetEntries() > 0) {
      Ang_Ener_PRAC->Draw("col");
   }
   c1->cd(3);
   if (hVxVy->GetEntries() > 0) {
      hVxVy->Draw("zcol");
   }

   TCanvas *c_ExEner = new TCanvas();
   if (hex->GetEntries() > 0) {
      hex->Draw();
      hex->GetXaxis()->SetTitle("Excitation Energy (MeV)");
      hex->GetYaxis()->SetTitle("Counts");
   }

   TCanvas *c_ExenerCorr = new TCanvas();
   c_ExenerCorr->Divide(2, 1);
   c_ExenerCorr->Draw();
   c_ExenerCorr->cd(1);
   hexCorr->Draw();
   c_ExenerCorr->cd(2);
   ExCorrvsZpos->Draw("zcol");

   TCanvas *c_AngDistr = new TCanvas();
   c_AngDistr->Divide(2, 1);
   c_AngDistr->cd(1);
   AngDistr->Draw();
   c_AngDistr->cd(2);
   AngDistrCM->Draw();
   AngDistrCM->GetXaxis()->SetTitle("Angle (deg)");
   AngDistrCM->GetYaxis()->SetTitle("#frac{d#sigma}{d#Omega} (a.u.)");

   TCanvas *c_ExvsZpos = new TCanvas();
   c_ExvsZpos->Divide(2, 1);
   c_ExvsZpos->cd(1);
   ExvsZpos->Draw("zcol");
   ExvsZpos->GetXaxis()->SetTitle("Excitation Energy (MeV)");
   ExvsZpos->GetYaxis()->SetTitle("z (cm)");
   c_ExvsZpos->cd(2);
   ExvsTrackLength->Draw("zcol");
   ExvsTrackLength->GetXaxis()->SetTitle("Excitation Energy (MeV)");
   ExvsTrackLength->GetYaxis()->SetTitle("Track Length (cm)");

   /*TCanvas *c_redchi2 = new TCanvas();
   c_redchi2->Divide(2, 1);
   c_redchi2->cd(1);
   hredchi2->Draw("zcol");
   c_redchi2->cd(2);
   hbredchi2->Draw("zcol");*/

   TCanvas *c_hex_segmented = new TCanvas();
   c_hex_segmented->Divide(3, 3);
   c_hex_segmented->cd(1);
   hex11->Draw();
   c_hex_segmented->cd(2);
   hex12->Draw();
   c_hex_segmented->cd(3);
   hex13->Draw();
   c_hex_segmented->cd(4);
   hex21->Draw();
   c_hex_segmented->cd(5);
   hex22->Draw();
   c_hex_segmented->cd(6);
   hex23->Draw();
   c_hex_segmented->cd(7);
   hex31->Draw();
   c_hex_segmented->cd(8);
   hex32->Draw();
   c_hex_segmented->cd(9);
   hex33->Draw();

   TCanvas *c_hex_segmented2 = new TCanvas();
   c_hex_segmented2->Divide(3, 3);
   c_hex_segmented2->cd(1);
   hex41->Draw();
   c_hex_segmented2->cd(2);
   hex42->Draw();
   c_hex_segmented2->cd(3);
   hex43->Draw();
   c_hex_segmented2->cd(4);
   hex51->Draw();
   c_hex_segmented2->cd(5);
   hex52->Draw();
   c_hex_segmented2->cd(6);
   hex53->Draw();
   c_hex_segmented2->cd(7);
   hex61->Draw();
   c_hex_segmented2->cd(8);
   hex62->Draw();
   c_hex_segmented2->cd(9);
   hex63->Draw();

   TCanvas *c_hex_vs_theta = new TCanvas();
   hexvstheta->Draw("zcol");
   hexvstheta->GetXaxis()->SetTitle("Energy (MeV)");
   hexvstheta->GetYaxis()->SetTitle("Angle (deg)");

   // last_20cm.close();
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
