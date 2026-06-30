#include <TString.h>

#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TTree.h"

#include <cstdio>
#include <fstream>
#include <iostream>
#include <set>
#include <vector>

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

double Ebin_max = 9.;
double Ebin_min = -1.;
int NumberBins = 80; // best chi2/ndf for 80 bins
const double Sn = 1.2181;
// Cut bounds (define once, use everywhere — avoids filename/cut mismatches)
double z_min = 2.0, z_max = 60.0;   // cm
double z_min1 = 2.0, z_max1 = 30.0;//cm
double T_min_cut = 5.0, T_max_cut = 15.0; // MeV, E_ej cut

double theta_excl_min_deg = 88.0, theta_excl_max_deg = 92.0;

void C16_pd_histogram()
{
   // gROOT->ProcessLine(".X /home/georgina/fair_install/ATTPCROOTv2/macro/a1975/myStyle.C");
   gStyle->SetOptStat(0);

   TH2F *Ang_Ener_Corr = new TH2F("Ang_Ener_Corr", "z<60cm", 120, 10, 40, 120, 0, 60.0);
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
   TH2F *Ang_Ener_Corr1 = new TH2F("Ang_Ener_Corr1", "z<30cm", 120, 10, 40, 120, 0, 60.0);
  

   auto *AngDistrCM = new TH1F("Ang_Distr_CM", "Ang_Distr_CM", NumberBins, 0, 180);
   auto *ExvsZpos = new TH2F("ExvsZpos", "ExvsZpos", NumberBins, Ebin_min, Ebin_max, 100, -5, 65);
   auto *ExCorrvsZpos = new TH2F("ExCorrvsZpos", "z<60cm", NumberBins, Ebin_min, Ebin_max, 100, 0, 65);
   auto *ExCorrvsZpos1 = new TH2F("ExCorrvsZpos1", "z<30cm", NumberBins, Ebin_min, Ebin_max, 100, 0, 35);
   auto *KineticEnergy = new TH1F("KineticEnergy", "KineticEnergy", 100, 0, 75);

   TH1F *h_PS_1n = new TH1F("h_PS_1n", "h_PS_1n", NumberBins, Ebin_min, Ebin_max);
   // auto *ExvsTrackLength = new TH2F("ExvsTrackLength", "ExvsTrackLength", 1000, -5, 15, 200, -20, 150);
   /*auto *hredchi2 = new TH1F("redchi2", "redchi2", 1000, 0, 0.0001);
   auto *hbredchi2 = new TH1F("bredchi2", "bredchi2", 1000, 0, 5);*/

   auto *hexvstheta_CM = new TH2F("hexVStheta_CM", "z<60cm", 180, 0, 180,NumberBins, Ebin_min, 20);
   auto *hexvstheta_CM1 = new TH2F("hexVStheta_CM1", "z<30cm", 180, 0, 180, NumberBins, Ebin_min, 20);
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

   Double_t m_p = 938.272076 / 1.0;     // masa nuclear del protón
   Double_t m_d = 1875.612931 / 1.0;    // masa nuclear del deuterón
   Double_t m_C15 = 13979.218707 / 1.0; // masa nuclear del 15C
   Double_t m_C16 = 14914.533798 / 1.0; // masa nuclear del 16C

   // Correct nuclear masses (MeV/c^2)
   // = atomic mass (u) * 931.494 - Z * 0.511 (electron mass)

   // Verify Q-value
   double Q = m_C16 + m_p - m_d - m_C15;
   std::cout << "Q-value = " << Q << " MeV " << std::endl;

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

   double densityH2 = 3.553e-5; // g/cm³
   AtTools::AtELossCATIMA elossH2(densityH2);

   double mass{16.0147};                        // Mass of C16 in u
   elossH2.SetMaterial(catima::Material(1, 1)); // Set material to H2

   elossH2.SetProjectile(16, 6, mass); // Set projectile to proton

   std::set<int> excluded = {111, 121, 148, 149};

   for (int i = 104; i <= 189; i++) {
      if (excluded.count(i))
         continue;
      char name[64];
      std::snprintf(name, sizeof(name), "run_%04d_2H.root", i);
      filenames.push_back(name);
   }

   std::cout << "Ebeam at the end of the TPC = " << elossH2.GetEnergy(Ebeam_buff, 100.0 * 10.0) << " MeV" << std::endl;

   // for (auto filename : filenames) {
   for (const TString &filename : filenames) {

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

      // TTree *physics = (TTree *)runFile->Get("parquettree");
      TTree *physics = (TTree *)runFile->Get("parquettree");

      // 2. SAFETY CHECK: Does 'parquettree' exist inside this file?
      if (!physics) {
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

      physics->SetBranchAddress("polar", &theta);   // rad
      physics->SetBranchAddress("azimuthal", &phi); // rad
      physics->SetBranchAddress("brho", &Brho);     // T·m
      physics->SetBranchAddress("redchisq", &redchi);
      physics->SetBranchAddress("vertex_z",
                                &zPos);     // zPos is in meters, I will convert it to cm when filling the histograms
      physics->SetBranchAddress("ke", &ke); // ke is in MeV
      Double_t vx_pos{},
         vy_pos{}; // vx_pos and vy_pos are in meters, I will convert them to cm when filling the histograms
      physics->SetBranchAddress("vertex_x", &vx_pos);
      physics->SetBranchAddress("vertex_y", &vy_pos);

      for (int i = 0; i < physics->GetEntries(); i++) {
         physics->GetEntry(i);

         Double_t p_ej = Brho * Z_ej * 2.99792458 / 10 * 1000;
         Double_t E_ej = (TMath::Sqrt(p_ej * p_ej + m_ej * m_ej) - m_ej); // p_ej: MeV/c, m_ej: MeV/c^2, E_ej: MeV

         auto [ex_energy, theta_cm] = kine_2b(m_C16, m_p, m_b, m_B, Ebeam_buff, theta, E_ej); // ke

         double dist3D = TMath::Sqrt(vx_pos * vx_pos + vy_pos * vy_pos + zPos * zPos) * 100.0; // cm
         // Double_t Ebeam_at_z = elossH2.GetEnergy(Ebeam_buff, zPos * 100); //
         Double_t Ebeam_at_z =
            elossH2.GetEnergy(Ebeam_buff, dist3D * 10.0); // Usar la distancia 3D para la corrección de energía, MeV/mm

         // cout << "Ebeam at position = " << dist3D << " cm: " << Ebeam_at_z << " MeV" << endl;

         // Corrección cinemática
         /*double kethe = 13.;
         double theta_lab_corr_tiltCorr = theta;
         (theta -
          (2.0 * TMath::Pi() / 4000) * (E_ej - kethe)); // theta: rad; theta_lab_corr: rad; E_ej-kethe: MeV 29.5
          */

         double theta_lab_corr = theta; // Por ahora sin corrección de ángulo, solo para probar la implementación de la
                                        // corrección de energía en el cálculo de Ex y theta_cm
         auto [ex_energy_corr, theta_cm_corr] = kine_2b(m_C16, m_p, m_b, m_B, Ebeam_at_z, theta_lab_corr, E_ej);

         /*auto [ex_energy_corr_tiltCorr, theta_cm_corr_tiltCorr] =
            kine_2b(m_C16, m_p, m_b, m_B, Ebeam_at_z, theta_lab_corr_tiltCorr,
                    E_ej); // energies: MeV, angles: radians
                    */

         KineticEnergy->Fill(E_ej); // Usar energía calibrada para el histograma de energía cinética

         double theta_deg = theta * TMath::RadToDeg();

         //Ang_Ener_tiltCorr->Fill(theta_lab_corr_tiltCorr * TMath::RadToDeg(), E_ej); // calibrado
        


         // Corte angular cerca de 90 deg (igual que en SPYRAL): descarta eventos
         // con theta_lab dentro de [88, 92] deg, mala reconstruccion cerca del eje del haz
         bool pass_theta_cut = (theta_deg < theta_excl_min_deg) || (theta_deg > theta_excl_max_deg);
         if (!pass_theta_cut){
            std::cout << "does not pass the cut" << std::endl;
          continue; // salta este evento, no rellena ningun histograma
         }


         // Fill corrected histogram

         if (zPos * 100 > z_min1 && zPos * 100 < z_max1){
            if (E_ej > T_min_cut && E_ej < T_max_cut) {
               hexCorr1->Fill(ex_energy_corr);  // Llenar el histograma con corrección de eficiencia
               ExCorrvsZpos1->Fill(ex_energy_corr, zPos * 100.0);
               hexvstheta_CM1->Fill(theta_cm_corr, ex_energy_corr); // theta_cm_corr is already in degrees,
                                                            // ex_energy_corr is in MeV

            }
              Ang_Ener_Corr1->Fill(theta_lab_corr * TMath::RadToDeg(),
                             E_ej); // theta lab!! -> I still have to implement the correction of catima?
         }

         if (zPos * 100 > z_min && zPos * 100 < z_max){


            Ang_Ener_Corr->Fill(theta_lab_corr * TMath::RadToDeg(),
                             E_ej); // theta lab!! -> I still have to implement the correction of catima?
                             
          if(E_ej > T_min_cut && E_ej < T_max_cut) {
            KineticEnergy->Fill(E_ej); // Usar energía calibrada para el histograma de energía cinética

            KEvsEx->Fill(E_ej, ex_energy_corr); // MeV, MeV sin corrección tilt

            // cm y MeV (zPos is in meters, E_ej is in MeV)
            ExCorrvsZpos->Fill(ex_energy_corr, zPos * 100.0); // MeV, cm
            ExvsZpos->Fill(ex_energy_corr, zPos * 100.0);              // MeV, cm
            hex->Fill(ex_energy_corr);               // Llenar el histograma con corrección de eficiencia
           
            hexCorr2->Fill(ex_energy_corr); // Llenar el histograma con corrección de eficiencia

            // Histograms

            Double_t vx = TMath::Sin(theta) * TMath::Sqrt(ke);
            Double_t vy = TMath::Cos(theta) * TMath::Sqrt(ke);

            hVxVy->Fill(vx, vy);

            AngDistr->Fill(theta * TMath::RadToDeg());
            AngDistrCM->Fill(theta_cm_corr);
            hexvstheta_CM->Fill(theta_cm_corr, ex_energy_corr);
            hexvstheta_lab->Fill(
               ex_energy_corr,
               theta_lab_corr *
                  TMath::RadToDeg()); // theta_lab_corr is in radians, convert to degrees for
                                      // the histogram ExvsTrackLength->Fill(ex_energy_corr, arclength);

            }
         }

         // tEvents->Fill();
      } // events
      runFile->Close();
      delete runFile;
   } // Files

   TFile *outFile = new TFile(Form("hexCorr2_pd_z_%g-%gcm_ke_%g-%gMeV_%dbins.root",
                                 z_min, z_max, T_min_cut, T_max_cut, NumberBins), "RECREATE");
   hexCorr2->Write();
   outFile->Close();

   TFile *outFile1 = new TFile(Form("hexCorr1_pd_z_%g-%gcm_ke_%g-%gMeV_%dbins.root",
                                 z_min, z_max1, T_min_cut, T_max_cut, NumberBins), "RECREATE");
   hexCorr1->Write();
   outFile1->Close();

   TFile * correction = new TFile(Form("hexvstheta_CM_pd_z_%g-%gcm_ke_%g-%gMeV_%dbins.root",
                                 z_min, z_max1, T_min_cut, T_max_cut, NumberBins), "RECREATE");
   hexvstheta_CM1->Write();
   correction->Close();

   //---------------- Combined Canvas ----------------//
   TCanvas *c_combined = new TCanvas("c_combined", "Combined Ang Energy", 1200, 800);
   c_combined->Divide(2, 1); // 2 columns, 2 rows

   // --- PAD 1: Ang_Ener_Corr (uncorrected kinematics) ---
   c_combined->cd(1);
   Ang_Ener_Corr->Draw("colz");
   Ang_Ener_Corr->GetXaxis()->SetTitle("#theta_{lab} (deg)");
   Ang_Ener_Corr->GetYaxis()->SetTitle("Kinetic Energy (MeV)");

   c_combined->cd(2);
   Ang_Ener_Corr1->Draw("colz");
   Ang_Ener_Corr1->GetXaxis()->SetTitle("#theta_{lab} (deg)");
   Ang_Ener_Corr1->GetYaxis()->SetTitle("Kinetic Energy (MeV)");

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
   /*c_combined->cd(2);
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
   leg2->Draw();*/
   c_combined->Update();

   TCanvas *c_ExvsZpos = new TCanvas("ExvsZpos", "Excited Energy spectra: comparison", 1200, 800);
   /*c_ExvsZpos->Divide(2, 1);
   c_ExvsZpos->cd(1);
   ExvsZpos->Draw("zcol");
   ExvsZpos->GetYaxis()->SetTitle("Z position (mm)");
   ExvsZpos->GetXaxis()->SetTitle("Excitation Energy (MeV)");
   c_ExvsZpos->cd(2);*/
   c_ExvsZpos->Divide(2, 1);
   c_ExvsZpos->cd(1);
   hexvstheta_CM->GetXaxis()->SetTitle("#theta_CM (deg)");
   hexvstheta_CM->GetYaxis()->SetTitle("Excitation Energy (MeV)");
   hexvstheta_CM->Draw("colz");
    c_ExvsZpos->cd(2);
   ExCorrvsZpos->GetYaxis()->SetTitle("Z position (cm)");
   ExCorrvsZpos->GetXaxis()->SetTitle("Excitation Energy (MeV)");
   ExCorrvsZpos->Draw("zcol");

    TCanvas *c_ExvsZpos1 = new TCanvas("ExvsZpos1", "Excited Energy spectra: z<=30 cm", 1200, 800);
   /*c_ExvsZpos->Divide(2, 1);
   c_ExvsZpos->cd(1);
   ExvsZpos->Draw("zcol");
   ExvsZpos->GetYaxis()->SetTitle("Z position (mm)");
   ExvsZpos->GetXaxis()->SetTitle("Excitation Energy (MeV)");
   c_ExvsZpos->cd(2);*/
   c_ExvsZpos1->Divide(2, 1);
   c_ExvsZpos1->cd(1);
   hexvstheta_CM1->GetXaxis()->SetTitle("#theta_CM (deg)");
   hexvstheta_CM1->GetYaxis()->SetTitle("Excitation Energy (MeV)");
   hexvstheta_CM1->Draw("colz");
    c_ExvsZpos1->cd(2);
   ExCorrvsZpos1->GetYaxis()->SetTitle("Z position (cm)");
   ExCorrvsZpos1->GetXaxis()->SetTitle("Excitation Energy (MeV)");
   ExCorrvsZpos1->Draw("zcol");

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

   TCanvas *kin = new TCanvas("kin", "kin", 1200, 800);
   KineticEnergy->Sumw2();
   kin->cd();
   KineticEnergy->Draw("E1");

   TCanvas *c_hex = new TCanvas("c_hex", "ChexCorr2", 1200, 800);
   c_hex->cd();
   hex->Draw("E1");


   TCanvas *c_hexvsthetaCM = new TCanvas("c_hexvsthetaCM","correlations", 1200,800);
   c_hexvsthetaCM->Divide(2,1);
   c_hexvsthetaCM->cd(1);
   hexvstheta_CM->Draw("colz"); //z <=60 cm
   c_hexvsthetaCM->cd(2);
   hexvstheta_CM1->Draw("colz"); //z <=30 cm
   // ---------------------------------------------------------
   // GUARDAR TODOS LOS CANVAS EN UN ÚNICO PDF MULTIPÁGINA
   // ---------------------------------------------------------

   // 1. ABRIR EL PDF: Imprimimos el PRIMER canvas y abrimos el archivo con "("
  TString pdfName = Form("plots_C16_z_%g-%gcm_ke_%g-%gMeV_%dbins.pdf",
                        z_min, z_max, T_min_cut, T_max_cut, NumberBins);
   // 2. ABRIR EL PDF: Le sumamos el "(" al nombre
   c_combined->Print(pdfName + "(");

   // 3. PÁGINAS INTERMEDIAS: Usamos el nombre tal cual
   c_ExvsZpos->Print(pdfName);
   c_ExvsZpos1->Print(pdfName);
   c_AngDistr->Print(pdfName);
   kin->Print(pdfName);

   // 4. CERRAR EL PDF: Le sumamos el ")" al final
   c_hex->Print(pdfName + ")");
}
