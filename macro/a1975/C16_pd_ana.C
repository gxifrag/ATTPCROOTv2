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

void C16_pd_ana()
{
   // FairRunAna *run = new FairRunAna();

   TH2F *Ang_Ener = new TH2F("Ang_Ener", "Ang_Ener", 600, 0, 60, 1000, 0, 20.0);
   TH2F *Ang_Ener_PRAC = new TH2F("Ang_Ener_PRAC", "Ang_Ener_PRAC", 1000, 0, 100, 1000, 0, 200.0);
   TH2F *ELossvsBrho = new TH2F("ELossvsBrho", "ELossvsBrho", 4000, 0, 25000, 1000, 0, 4);
   TH2F *dedxvsBrho = new TH2F("dedxvsBrho", "dedxvsBrho", 4000, 0, 10000, 1000, 0, 4);
   TH2F *hVxVy = new TH2F("hVxVy", "hVxVy", 1000, 0, 4, 1000, 0, 4);
   TH1F *henergyIC = new TH1F("henergyIC", "henergyIC", 2048, 0, 2047);

   auto *hex = new TH1F("hex", "hex", 80, -2, 8);
   auto *QvsEb = new TH2F("QvsEb", "QvsEb", 1000, -5, 15, 100, 0, 300);
   auto *QvsZpos = new TH2F("QvsZpos", "QvsZpos", 1000, -10, 50, 200, -100, 100);
   auto *hexCorr = new TH1F("hexCorr", "hexCorr", 80, -2, 8);

   auto *AngDistr = new TH1F("Ang_Distr", "Ang_Distr", 128, 0, 120);
   auto *AngDistrCM = new TH1F("Ang_Distr_CM", "Ang_Distr_CM", 90, 0, 180);
   auto *ExvsZpos = new TH2F("ExvsZpos", "ExvsZpos", 1000, -5, 15, 200, -20, 150);
   auto *ExvsTrackLength = new TH2F("ExvsTrackLength", "ExvsTrackLength", 1000, -5, 15, 200, -20, 150);
   auto *ExCorrvsZpos = new TH2F("ExCorrvsZpos", "ExCorrvsZpos", 1000, -10, 10, 200, -100, 100);

   /*auto *hredchi2 = new TH1F("redchi2", "redchi2", 1000, 0, 0.0001);
   auto *hbredchi2 = new TH1F("bredchi2", "bredchi2", 1000, 0, 5);*/
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

   /*auto *hex20 = new TH1F("hex0-19", "hex0-19", 100, -5, 20);
double densityH2 = 3.3084e-5; // g/cm³
   AtTools::AtELossCATIMA elossH2(densityH2);
   double mass{16.0147}; // Mass of C16 in u
   elossH2.SetMaterial(catima::Material(1, 1)); // Set material to H2
   elossH2.SetProjectile(16 6, mass);           // Set projectile to proton
   auto *hex40 = new TH1F("hex20-39", "hex20-39", 100, -5, 20);
   auto *hex60 = new TH1F("hex40-59", "hex40-59", 100, -5, 20);
   auto *hex80 = new TH1F("hex60-79", "hex60-79", 100, -5, 20);
   auto *hex100 = new TH1F("hex80-99", "hex80-99", 100, -5, 20);*/

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
   Double_t Ebeam_buff = 11.5 * 16;
   cout<<" Beam energy in buffer gas : "<<Ebeam_buff<<"\n";
   //Double_t Ebeam_buff = 10.6* 16; // MeV
   Double_t m_b = m_d;
   Double_t m_B = m_C15;

   // Ejectile parameters:deuterium
   int A_ej = 2;
   int Z_ej = 1;
   Double_t m_ej = m_d;

   std::vector<TString> filenames;

   // ELoss tables.
   //AtTools::AtELossTable *elossTableH2 = new AtTools::AtELossTable();
   //elossTableH2->LoadSrimTable("StoppingPower_SRIM_C16_H2.txt"); //SRIM no me va.
   //elossTableH2->LoadLiseTable("StoppingPower_C16_H2.txt", 2.0158,3.3084e-5);

   double densityH2 = 3.3084e-5; // g/cm³
   AtTools::AtELossCATIMA elossH2(densityH2);
   double mass{16.0147}; // Mass of C16 in u
   elossH2.SetMaterial(catima::Material(1, 1)); // Set material to H2
   elossH2.SetProjectile(16, 6, mass);           // Set projectile to proton

   double kethe= 13.;     
   
   filenames.push_back("run_0110_2H.root");
   //filenames.push_back("run_0111_1H.root");
   filenames.push_back("run_0112_2H.root");
   filenames.push_back("run_0113_2H.root");
   filenames.push_back("run_0114_2H.root");
   filenames.push_back("run_0115_2H.root");
   filenames.push_back("run_0116_2H.root");
   filenames.push_back("run_0117_2H.root");
   filenames.push_back("run_0118_2H.root");
   filenames.push_back("run_0119_2H.root");
   filenames.push_back("run_0120_2H.root");
   // filenames.push_back("run_0121_1H.root");
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
   // filenames.push_back("run_0144_1H.root");
   filenames.push_back("run_0145_2H.root");
   filenames.push_back("run_0146_2H.root");
   filenames.push_back("run_0147_2H.root");
   // filenames.push_back("run_0148_1H.root");
   filenames.push_back("run_0149_2H.root");
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
   // filenames.push_back("run_0165_1H.root");
   filenames.push_back("run_0166_2H.root");
   filenames.push_back("run_0167_2H.root");
   filenames.push_back("run_0168_2H.root");
   filenames.push_back("run_0169_2H.root");
   filenames.push_back("run_0170_2H.root");
   // filenames.push_back("run_0171_1H.root");
   filenames.push_back("run_0172_2H.root");
   filenames.push_back("run_0173_2H.root");
   filenames.push_back("run_0174_2H.root");
   filenames.push_back("run_0175_2H.root");
   filenames.push_back("run_0176_2H.root");
   filenames.push_back("run_0177_2H.root");
   // filenames.push_back("run_0178_1H.root");
   filenames.push_back("run_0179_2H.root");
   // filenames.push_back("run_0180_1H.root");
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
      TFile *runFile =
         new TFile("/home/georgina/C16_analysis/C16_H2/C16_pd/InterpolationSolver_pd_root/" + filename, "R");
      TTree *Tphysics = (TTree *)runFile->Get("parquettree");

      Double_t theta{};
      Double_t phi{};
      Double_t Brho{};
      //Double_t redchi{};
      Double_t zPos{};
      Tphysics->SetBranchAddress("polar", &theta);
      Tphysics->SetBranchAddress("azimuthal", &phi);
      Tphysics->SetBranchAddress("brho", &Brho);
      //Tphysics->SetBranchAddress("redchisq", &redchi);
      Tphysics->SetBranchAddress("vertex_z", &zPos);

      for (int i = 0; i < Tphysics->GetEntries(); i++) {
         Tphysics->GetEntry(i);

         Double_t p_ej = Brho * Z_ej * 2.99792458 / 10 * 1000;
         Double_t E_ej = TMath::Sqrt(p_ej * p_ej + m_ej * m_ej) - m_ej;

         auto [ex_energy, theta_cm] = kine_2b(m_C16, m_p, m_b, m_B, Ebeam_buff, theta, E_ej);

         // Beam energy corrected for energy loss along the target at vertex position
         Double_t Ebeam_at_z = elossH2.GetEnergy(Ebeam_buff, zPos * 100); // 
         
         std::cout<<" Beam energy at z : "<<Ebeam_at_z<< "  " <<  "z:" << zPos*100 << "\n";

         auto [ex_energy_corr, theta_cm_corr]= kine_2b(m_C16, m_p, m_b, m_B,Ebeam_at_z, theta, E_ej);

//         double theta_corr_JR=(theta-(2.0*3.14159/4000) * (E_ej - kethe)) * TMath::RadToDeg(); //theta: rad, E_ej-kethe: MeV
         //cout << " theta_corr_JR: " << theta_corr_JR << "\n";

         /*if (theta * TMath::RadToDeg() > 90 || theta * TMath::RadToDeg() < 12.0)
            continue;*/

         // Fill uncorrected histogram
         hex->Fill(ex_energy);      
         ExvsZpos->Fill(ex_energy, zPos*100.0);

          // Fill corrected histogram
         /*if (zPos*100> 2.0 && zPos*100 < 60.0)  // only consider reactions occuring within the target region
         {
            cout << " Ex corrected : " << ex_energy_corr << "  " << " zpos: " << zPos*100.0 << " Ebeam at z: " << Ebeam_at_z << "\n";
            
            ExCorrvsZpos->Fill(ex_energy_corr, zPos*100.0);
            hexCorr->Fill(ex_energy_corr-0.8);
         }*/

         // Histograms
         //hredchi2->Fill(redchi);

         Ang_Ener->Fill(theta* TMath::RadToDeg(), E_ej); //theta lab!! -> I still have to implement the correction of catima?

         Double_t vx = TMath::Sin(theta) * TMath::Sqrt(E_ej);
         Double_t vy = TMath::Cos(theta) * TMath::Sqrt(E_ej);

         hVxVy->Fill(vx, vy);

         AngDistr->Fill(theta * TMath::RadToDeg());
         AngDistrCM->Fill(theta_cm);

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

         /*if(theta_cm > 0 && theta_cm <= 20) hex20->Fill(ex_energy);
         if(theta_cm > 20 && theta_cm <= 40) hex40->Fill(ex_energy);
         if(theta_cm > 40 && theta_cm <= 60) hex60->Fill(ex_energy);
         if(theta_cm > 60 && theta_cm <= 80) hex80->Fill(ex_energy);
         if(theta_cm > 80 && theta_cm <= 100) hex100->Fill(ex_energy);*/

         hexvstheta->Fill(ex_energy, theta * TMath::RadToDeg());
      } // events
   } // Files

   AngDistrCM->Divide(new TF1("sin", "sin(x * TMath::DegToRad())", 0, 180));

   // Kinematics
   /*Double_t *ThetaCMS = new Double_t[20000];
   Double_t *ThetaLabRec = new Double_t[20000];
   Double_t *EnerLabRec = new Double_t[20000];
   Double_t *ThetaLabSca = new Double_t[20000];
   Double_t *EnerLabSca = new Double_t[20000];
   Double_t *MomLabRec = new Double_t[20000];

C16_pd_C15_3103keV_Ebeam11_5.txt
C16_pd_C15_4780keV_Ebeam11_5.txt
C16_pd_C15_6841keV_Ebeam11_5.txt



   TString fileKine =
      "/home/georgina/fair_install/ATTPCROOTv2/macro/Kinematics/Decay_kinematics/C16_pd_C15_gs_Ebeam11_5.txt";
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

   Double_t *ThetaCMS2 = new Double_t[20000];
   Double_t *ThetaLabRec2 = new Double_t[20000];
   Double_t *EnerLabRec2 = new Double_t[20000];
   Double_t *ThetaLabSca2 = new Double_t[20000];
   Double_t *EnerLabSca2 = new Double_t[20000];
   Double_t *MomLabRec2 = new Double_t[20000];

   TString fileKine2 =
      "/home/georgina/fair_install/ATTPCROOTv2/macro/Kinematics/Decay_kinematics/C16_pd_C15_740keV_Ebeam11_5.txt";
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
      std::cout << " Warning : No Kinematics file found for this reaction!" << std::endl;

   TGraph *kine_1st = new TGraph(numKin2, ThetaLabRec2, EnerLabRec2);*/

   // End Kinematics


   // KINEMATICS FOR DIFFERENT EXCITATION ENERGIES
std::vector<std::string> files = {
   "C16_pd_C15_gs_Ebeam11_5.txt",
   "C16_pd_C15_740keV_Ebeam11_5.txt",
   "C16_pd_C15_3103keV_Ebeam11_5.txt",
   "C16_pd_C15_4780keV_Ebeam11_5.txt",
   "C16_pd_C15_6841keV_Ebeam11_5.txt"
};

std::vector<std::string> labels = {
   "gs", "740keV", "3103keV", "4780keV", "6841keV"
};

// Colors for each line (ROOT color codes: 2=red,4=blue,8=green, etc.)
std::vector<int> colors = {kBlack, kRed, kBlue, kGreen+2, kMagenta};

std::vector<TGraph*> graphs;

for (size_t i = 0; i < files.size(); i++) {
   TString fileKine = Form("/home/georgina/fair_install/ATTPCROOTv2_3/macro/Kinematics/Decay_kinematics/%s", files[i].c_str());
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


   /*TCanvas *c1 = new TCanvas();
   c1->Divide(2, 2);
   c1->Draw();
   c1->cd(1);
   Ang_Ener->SetMarkerStyle(20);
   Ang_Ener->SetMarkerSize(0.5);
   Ang_Ener->Draw("col");
   Ang_Ener->GetXaxis()->SetTitle("Angle (deg)");
   Ang_Ener->GetYaxis()->SetTitle("Energy (MeV)");
   kine_gs->Draw("SAME");
   kine_1st->Draw("SAME");*/


   /*c1->cd(2);
   Ang_Ener_PRAC->Draw("col");
   c1->cd(3);
   hVxVy->Draw("zcol");*/

TCanvas *c1 = new TCanvas();
Ang_Ener->SetMarkerStyle(20);
Ang_Ener->SetMarkerSize(0.5);
Ang_Ener->Draw("col");
Ang_Ener->GetXaxis()->SetTitle("Angle (deg)");
Ang_Ener->GetYaxis()->SetTitle("Energy (MeV)");

// Draw all kinematics graphs from the loop
/*for (size_t i = 0; i < graphs.size(); i++) {
    graphs[i]->Draw("L SAME");  // "L SAME" draws as a line on the same canvas
}

// Optional: add a legend
auto legend = new TLegend(0.65, 0.65, 0.88, 0.88);
for (size_t i = 0; i < graphs.size(); i++) {
    legend->AddEntry(graphs[i], labels[i].c_str(), "l");
}
legend->Draw();*/


   TCanvas *c_ExEner = new TCanvas();
   hex->Draw();
   hex->GetXaxis()->SetTitle("Excitation Energy (MeV)");
   hex->GetYaxis()->SetTitle("Counts");

   /*TCanvas *c_ExenerCorr = new TCanvas();
   c_ExenerCorr->Divide(2, 1);
   c_ExenerCorr->Draw();
   c_ExenerCorr->cd(1);
   hexCorr->Draw();
   c_ExenerCorr->cd(2);
   ExCorrvsZpos->Draw("zcol");

   TCanvas *c_ExenerCorr1 = new TCanvas();
   c_ExenerCorr->Draw();
   hexCorr->Draw();*/


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

   /*TCanvas *c_hex_segmented = new TCanvas();
   c_hex_segmented->Divide(3, 2);
   c_hex_segmented->cd(1);
   hex20->Draw();
   c_hex_segmented->cd(2);
   hex40->Draw();
   c_hex_segmented->cd(3);
   hex60->Draw();
   c_hex_segmented->cd(4);
   hex80->Draw();
   c_hex_segmented->cd(5);
   hex100->Draw();*/

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
