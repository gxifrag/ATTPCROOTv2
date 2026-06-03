#include <TRandom.h>

TH1F *hpstot;
TH1F *hps1;
TH1F *hps2;
TH1F *hps3;
TH1F *hps4;
TH1F *hps5;
TH2F *hp2d;
TRandom1 rGen(0);
TRandom1 rGen1(0);

double Ex_calculation(double ELab, double ThetaLab, double mass_in_MeV[4], double T_b_A)
{ // thetalab: degrees, Tb=11,5MeV/A

   double mbeam, mtarget, mejectile, mrecoil, Tbeam;

   double MeV_to_uma = 1.0 / 931.494061;

   mbeam = mass_in_MeV[0];
   mtarget = mass_in_MeV[1];
   mejectile = mass_in_MeV[2];
   mrecoil = mass_in_MeV[3];

   Tbeam = T_b_A * mbeam * MeV_to_uma;

   double Ebeam, pbeam, pejectile, C1, Ex;

   Ebeam = mbeam + Tbeam;
   pbeam = sqrt(2 * mbeam * Tbeam + Tbeam * Tbeam);
   pejectile = sqrt(2 * mejectile * ELab + ELab * ELab);
   C1 = mbeam * mbeam + mtarget * mtarget + mejectile * mejectile + 2 * Ebeam * mtarget;

   Ex = sqrt(C1 - 2 * (Ebeam + mtarget) * (mejectile + ELab) +
             2 * pbeam * pejectile * cos(ThetaLab * 3.14159265359 / 180.0)) -
        mrecoil;

   return Ex;
}

void Ex_and_ThetaCM_calculation_v2(double ELab, double ThetaLab, double mass_in_MeV[4], double T_b_A,
                                   double Ex_ThetaCM[2])
{

   double ThetaCM;

   double Eex_recoil = Ex_calculation(ELab, ThetaLab, mass_in_MeV, T_b_A);

   double T_b, T3, E3, E4, E1, Et, Et_cm, E3_cm, p1, p3_cm, p4, p3_x, p3_y, p3, theta, phi;
   double m1, m2, m3, m4, gamma, beta, cos_theta_CM;

   double MeV_to_uma = 1.0 / 931.494061;

   m1 = mass_in_MeV[0];
   m2 = mass_in_MeV[1];
   m3 = mass_in_MeV[2];
   m4 = mass_in_MeV[3] + Eex_recoil; // we need ti use the one excited

   E3 = ELab + m3;

   T_b = T_b_A * m1 * MeV_to_uma;
   E1 = m1 + T_b;
   p1 = sqrt(pow(E1, 2) - pow(m1, 2));
   Et = E1 + m2;
   Et_cm = sqrt(2 * E1 * m2 + pow(m2, 2) + pow(m1, 2));
   E3_cm = 0.5 * (Et_cm + (pow(m3, 2) - pow(m4, 2)) / Et_cm);
   p3_cm = sqrt(pow(E3_cm, 2) - pow(m3, 2));

   gamma = Et / Et_cm;
   beta = p1 / Et;

   cos_theta_CM = ((E3 / gamma) - E3_cm) / (beta * p3_cm);

   ThetaCM = 180.0 - acos(cos_theta_CM) * 180.0 / TMath::Pi();

   Ex_ThetaCM[0] = Eex_recoil;
   Ex_ThetaCM[1] = ThetaCM;
}

void PhaseSpace_1n_16Cpt14C_binWidth(int Nbin, double MinLim, double MaxLim, double sigma_Ex, double thetaCM_min,
                                     double thetaCM_max, double theta_min, double theta_max, double T_min, double T_max,
                                     int SavingDataFile, TString OutputFileName_1n)
{

   if (!gROOT->GetClass("TGenPhaseSpace"))
      gSystem->Load("libPhysics");

   // Change masses + beam kinetic energy
   Double_t m1, m2, m3, m4, m6, m7, mn, bE, T1, p1, m5, md;

   m1 = 16.014701 * 931.494;             // Beam mass 16C
   m2 = 1.007825 * 931.494;              // proton mass: target
   m3 = 3.016049 * 931.494;              // Ejectile mass t
   m4 = 16.014701 * 931.494;             // Beam-like fragment 16C
   m5 = 15.010599256 * 931.494;          // Beam-like fragment 15C
   m6 = 14.00324198843 * 931.494;        // Beam mass 14C
   m7 = 13.0033548378 * 931.494;         // Beam mass 13C
   md = 2.014101 * 931.494;              // proton mass
   mn = 1.008665 * 931.494;              // Neutron mass
   bE = 11.5;                            // Beam energy MeV/A
   T1 = 16.014701 * bE;                  // Beam kinetic energy middle of target
   p1 = TMath::Sqrt(T1 * (T1 + 2 * m1)); // Beam momentum

   double mass_in_MeV[4];
   double Ex_ThetaCM[2];

   mass_in_MeV[0] = m1;
   mass_in_MeV[1] = m2;
   mass_in_MeV[2] = m3;
   mass_in_MeV[3] = m6;

   // px,py,pz,E in MeV
   TLorentzVector target(0.0, 0.0, 0.0, m2);
   TLorentzVector beam(0.0, 0.0, p1, T1 + m1);
   TLorentzVector W = beam + target; // LAB reference frame

   std::cout << "W mass = " << W.M() << ", sum of final masses = " << (md + m5 + mn) << std::endl;

   Double_t masses[3] = {m3, m7, mn};
   std::cout << "md = " << md << ", m5 = " << m5 << ", mn = " << mn << std::endl;

   // Set your Phase Space to 3 particles with masses given by array masses.
   TGenPhaseSpace event;
   event.SetDecay(W, 3, masses);

   //------------------------------------------------------------------------

   // Q-value de la reacción
   double Q = (m1 + m2) - (m3 + m7 + mn);
   cout << "Q-value: " << Q << " MeV" << endl;

   //------------------------------------------------------------------------------
   double etot, etot_cm, gam, beta, e3, p3, e3_cm, e3n, t3, poxy, poxyth, brho, m4ex, theta;

   etot = T1 + m1 + m2;
   etot_cm = TMath::Sqrt(m1 * m1 + m2 * m2 + 2 * m2 * (T1 + m1));
   gam = etot / etot_cm;
   beta = TMath::Sqrt(1 - 1 / gam / gam);

   double weight, theta_deg, ep;

   TH1F *h_PS_1n = new TH1F("h_PS_1n", "Phase Space: ^{16}C + p #rightarrow t + ^{13}C + 1n", Nbin, MinLim, MaxLim);
   TH2F *h_kn_1n = new TH2F("h_kn_1n", "h_kn: ^{16}C + p #rightarrow t + ^{13}C + 1n", 180, 0, 50, 600, 0, 60);

   TH1F *h_theta_debug1n = new TH1F("h_theta_debug1n", "#theta_{lab} de eventos generados", 180, 0, 180);
   TH1F *h_t3_debug1n = new TH1F("h_t3_debug1n", "kinetic energy deuteron", 600, 0, 100);

   TFile *EfficiencySimFile;
   TTree *out_tree;
   double ThetaCM_cal, Ex_cal, Weight_sim;

   if (SavingDataFile) {
      EfficiencySimFile = new TFile(OutputFileName_1n, "RECREATE");
      out_tree = new TTree("simulated_tree", " mysim tree ");
      TBranch *ThetaCM_cal_b = out_tree->Branch("ThetaCM_cal", &ThetaCM_cal, "ThetaCM_cal/D");
      TBranch *Ex_cal_b = out_tree->Branch("Ex_cal", &Ex_cal, "Ex_cal/D");
      TBranch *Weight_sim_b = out_tree->Branch("Weight_sim", &Weight_sim, "Weight_sim/D");
   }
   int I_percent = 0;
   int Nsim = 10000000;
   int factor_100 = Nsim / 100;
   for (Int_t n = 0; n < Nsim; n++) {

      if (n % factor_100 == 0) {
         cout << I_percent << " /100 progress " << endl;
         I_percent += 1;
      }
      // Generates a Random final state.
      weight = event.Generate();
      if (weight <= 0)
         continue;

      // Returns the Lorentz Vector for a final state.
      TLorentzVector *pTriton = event.GetDecay(0);
      TLorentzVector *p13C = event.GetDecay(1);
      TLorentzVector *pNeutron = event.GetDecay(2);

      e3n = pTriton->Energy(); // From Lorentz vector
      t3 = e3n - m3;           // MeV resolution calculated from the 3 alpha
      // cout <<  " e3n: " << e3n << " t3: " << t3 << endl ;

      theta_deg = pTriton->Theta() * 180. / TMath::Pi(); // lab angle in degrees
      // cout << " theta_deg: " << theta_deg << endl ;

      Ex_and_ThetaCM_calculation_v2(t3, theta_deg, mass_in_MeV, bE, Ex_ThetaCM);

      if (theta_deg >= theta_min && theta_deg < theta_max && t3 > T_min && t3 < T_max && Ex_ThetaCM[1] < thetaCM_max &&
          Ex_ThetaCM[1] > thetaCM_min) {
         h_theta_debug1n->Fill(theta_deg);
         h_t3_debug1n->Fill(t3);

         double Ex_res = rGen.Gaus(Ex_ThetaCM[0], sigma_Ex);
         h_PS_1n->Fill(Ex_res, weight);
         h_kn_1n->Fill(theta_deg, t3, weight);

         if (SavingDataFile) {
            ThetaCM_cal = Ex_ThetaCM[1];
            Ex_cal = Ex_res;
            Weight_sim = weight;
            out_tree->Fill();
         }
      }
   }

   TCanvas *C1_1n = new TCanvas();
   C1_1n->cd();
   h_kn_1n->GetXaxis()->SetTitle("#theta_{lab} [deg]");
   h_kn_1n->GetYaxis()->SetTitle("Kinetic Energy T_{d} [MeV]");
   h_kn_1n->Draw("colz");

   TCanvas *C3_1n = new TCanvas();
   C3_1n->cd();
   h_PS_1n->GetXaxis()->SetTitle("Excitation Energy E_{x} [MeV]");
   h_PS_1n->GetYaxis()->SetTitle("Weighted Counts");
   h_PS_1n->Draw("hist");

   TCanvas *C_t3_debug1n = new TCanvas("C_t3_debug1n", "kinetic energy triton", 800, 600);
   h_t3_debug1n->GetXaxis()->SetTitle("T_{3} [MeV]");
   h_t3_debug1n->Draw("hist");

   TCanvas *C_theta_debug1n = new TCanvas("C_theta_debug1n", "theta", 800, 600);
   h_theta_debug1n->GetXaxis()->SetTitle("theta");
   h_theta_debug1n->Draw("hist");

   if (SavingDataFile) {
      cout << "writting the data file" << endl;
      out_tree->Write();
      EfficiencySimFile->Close();
   }
}

void PhaseSpace_2n_16Cpt14C_binWidth(int Nbin, double MinLim, double MaxLim, double sigma_Ex, double thetaCM_min,
                                     double thetaCM_max, double theta_min, double theta_max, double T_min, double T_max,
                                     int SavingDataFile, TString OutputFileName_2n)
{

   if (!gROOT->GetClass("TGenPhaseSpace"))
      gSystem->Load("libPhysics");

   // Change masses + beam kinetic energy
   Double_t m1, m2, m3, m4, m6, m7, m8, mn, md, bE, T1, p1, m5;

   m1 = 16.014701 * 931.494;             // Beam mass 16C
   m2 = 1.007825 * 931.494;              // proton
   m3 = 3.016049 * 931.494;              // Ejectile mass t
   m4 = 16.014701 * 931.494;             // Beam-like fragment 16C
   m5 = 15.010599256 * 931.494;          // Beam-like fragment 15C
   m6 = 14.00324198843 * 931.494;        // Beam mass 14C
   m7 = 13.0033548378 * 931.494;         // Beam mass 13C
   m8 = 12.000000 * 931.494;             // Beam mass 12C
   mn = 1.008665 * 931.494;              // Neutron mass
   md = 2.014101 * 931.494;              // deuteron mass: target
   bE = 11.5;                            // Beam energy MeV/A
   T1 = 16.014701 * bE;                  // Beam kinetic energy middle of target
   p1 = TMath::Sqrt(T1 * (T1 + 2 * m1)); // Beam momentum

   double mass_in_MeV[4];
   double Ex_ThetaCM[2];

   mass_in_MeV[0] = m1;
   mass_in_MeV[1] = m2;
   mass_in_MeV[2] = m3;
   mass_in_MeV[3] = m6;

   // px,py,pz,E in MeV
   TLorentzVector target(0.0, 0.0, 0.0, m2);
   TLorentzVector beam(0.0, 0.0, p1, T1 + m1);
   TLorentzVector W = beam + target; // LAB reference frame

   Double_t masses[4] = {m3, m8, mn, mn};
   //------------------------------------------------------------------------

   // Q-value de la reacción
   double Q = (m1 + m2) - (m8 + m3 + 2 * mn);
   cout << "Q-value: " << Q << " MeV" << endl;

   //------------------------------------------------------------------------------

   // Set your Phase Space to 3 particles with masses given by array masses.
   TGenPhaseSpace event;
   event.SetDecay(W, 4, masses);

   double etot, etot_cm, gam, beta, e3, p3, e3_cm, e3n, t3, poxy, poxyth, brho, m4ex, theta;

   etot = T1 + m1 + m2;
   etot_cm = TMath::Sqrt(m1 * m1 + m2 * m2 + 2 * m2 * (T1 + m1));
   gam = etot / etot_cm;
   beta = TMath::Sqrt(1 - 1 / gam / gam);

   double weight, theta_deg, ep;

   TH1F *h_PS_2n = new TH1F("h_PS_2n", "Phase Space:^{16}C + p #rightarrow t + ^{12}C + 2n", Nbin, MinLim, MaxLim);
   TH2F *h_kn_2n = new TH2F("h_kn_2n", "h_kn:^{16}C + p #rightarrow t + ^{12}C + 2n", 180, 0, 50, 600, 0, 60);

   TH1F *h_theta_debug2n = new TH1F("h_theta_debug2n", "#theta_{lab} de eventos generados", 180, 0, 180);
   TH1F *h_t3_debug2n = new TH1F("h_t3_debug2n", "kinetic energy deuteron", 600, 0, 100);

   TFile *EfficiencySimFile; // = new TFile("./16Cpp_phasespace/PhaseSpace_16C_pp_1n.root","recreate");
   TTree *out_tree;
   double ThetaCM_cal, Ex_cal, Weight_sim;

   if (SavingDataFile) {
      EfficiencySimFile = new TFile(OutputFileName_2n, "recreate");
      out_tree = new TTree("simulated_tree", " mysim tree ");
      TBranch *ThetaCM_cal_b = out_tree->Branch("ThetaCM_cal", &ThetaCM_cal, "ThetaCM_cal/D");
      TBranch *Ex_cal_b = out_tree->Branch("Ex_cal", &Ex_cal, "Ex_cal/D");
      TBranch *Weight_sim_b = out_tree->Branch("Weight_sim", &Weight_sim, "Weight_sim/D");
   }

   int I_percent = 0;
   int Nsim = 10000000;
   int factor_100 = Nsim / 100;
   for (Int_t n = 0; n < Nsim; n++) {
      if (n % factor_100 == 0) {
         cout << I_percent << " /100 progress " << endl;
         I_percent += 1;
      }
      // Generates a Random final state.
      weight = event.Generate();

      // Returns the Lorentz Vector for a final state.
      TLorentzVector *pTriton = event.GetDecay(0);
      TLorentzVector *p12C = event.GetDecay(1);

      TLorentzVector *pNeut = event.GetDecay(2);
      TLorentzVector *pNeut2 = event.GetDecay(3);

      e3n = pTriton->Energy(); // From Lorentz vector
      t3 = e3n - md;           // MeV resolution calculated
      theta_deg = pTriton->Theta() * 180. / TMath::Pi();
      Ex_and_ThetaCM_calculation_v2(t3, theta_deg, mass_in_MeV, bE, Ex_ThetaCM);

      if (theta_deg >= theta_min && theta_deg < theta_max && t3 > T_min && t3 < T_max && Ex_ThetaCM[1] < thetaCM_max &&
          Ex_ThetaCM[1] > thetaCM_min) {
         h_theta_debug2n->Fill(theta_deg);
         h_t3_debug2n->Fill(t3);

         double Ex_res = rGen.Gaus(Ex_ThetaCM[0], sigma_Ex);
         h_PS_2n->Fill(Ex_res, weight);
         h_kn_2n->Fill(theta_deg, t3, weight);

         if (SavingDataFile) {
            ThetaCM_cal = Ex_ThetaCM[1];
            Ex_cal = Ex_res;
            Weight_sim = weight;
            out_tree->Fill();
         }
      }
   }

   TCanvas *C1_2n = new TCanvas();
   C1_2n->cd();
   h_kn_2n->Draw("colz");

   TCanvas *C3_2n = new TCanvas();
   C3_2n->cd();
   h_PS_2n->Draw("hist");

   TCanvas *C_t3_debug2n = new TCanvas("C_t3_debug2n", "kinetic energy triton", 800, 600);
   h_t3_debug2n->GetXaxis()->SetTitle("T_{3} [MeV]");
   h_t3_debug2n->Draw("hist");

   TCanvas *C_theta_debug2n = new TCanvas("C_theta_debug2n", "theta", 800, 600);
   h_theta_debug2n->GetXaxis()->SetTitle("theta");
   h_theta_debug2n->Draw("hist");

   if (SavingDataFile) {
      cout << "writting the data file" << endl;
      out_tree->Write();
      EfficiencySimFile->Close();
   }
}

void PhaseSpace_16Cpt14C_test()
{
   double MinLim, MaxLim;
   int Nbin;

   MinLim = -4.;
   MaxLim = 14.;
   Nbin = 180;

   // Acceptance and cuts
   // theta_min=0, theta_max=41, T_min=0, T_max=480*2: from LISE++ calculation
   double theta_min = 10;
   double theta_max = 60;

   double thetaCM_min = 0.0;
   double thetaCM_max = 180.0;

   double T_min = 0.0;
   double T_max = 70 * 16; // 960

   double sigma_Ex = 0.200; // This can be modified in the function if the resolution dependes on the Ex energy (or
                            // directly include the E anf angle resolution if preferred)

   int SavingDataFile = 1;

   TString OutputFileName_1n = "PhaseSpace_16C_pt_1n.root";
   PhaseSpace_1n_16Cpt14C_binWidth(Nbin, MinLim, MaxLim, sigma_Ex, thetaCM_min, thetaCM_max, theta_min, theta_max,
                                   T_min, T_max, SavingDataFile, OutputFileName_1n);

   TString OutputFileName_2n = "PhaseSpace_16C_pt_2n.root";
   PhaseSpace_2n_16Cpt14C_binWidth(Nbin, MinLim, MaxLim, sigma_Ex, thetaCM_min, thetaCM_max, theta_min, theta_max,
                                   T_min, T_max, SavingDataFile, OutputFileName_2n);
}
