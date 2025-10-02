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


void C16_pd_ana_v16_17Sep()
{

   bool guardar_en_pdf = true; // ← cambia a false si quieres solo verlos en pantalla
   // Activar modo batch si estás guardando en PDF
   if (guardar_en_pdf) {
      gROOT->SetBatch(kTRUE); // ← esto evita que se abran ventanas
   }
   // FairRunAna *run = new FairRunAna();

   TH2F *Ang_Ener = new TH2F("Ang_Ener", "Ang_Ener", 720, 10, 60, 1000, 0, 14.0);
   TH2F *Ang_Ener_Corr = new TH2F("Ang_Ener_Corr", "Ang_Ener_Corr",720, 10, 60, 1000, 0, 14.0);

   TH2F *ELossvsBrho = new TH2F("ELossvsBrho", "ELossvsBrho", 4000, 0, 25000, 1000, 0, 4);
   TH2F *dedxvsBrho = new TH2F("dedxvsBrho", "dedxvsBrho", 4000, 0, 10000, 1000, 0, 4);
   TH2F *hVxVy = new TH2F("hVxVy", "hVxVy", 1000, 0, 4, 1000, 0, 4);
   TH1F *henergyIC = new TH1F("henergyIC", "henergyIC", 2048, 0, 2047);

   auto *hex = new TH1F("hex", "hex", 90, -4, 14);
   auto *QvsEb = new TH2F("QvsEb", "QvsEb", 1000, -5, 15, 100, 0, 300);
   auto *QvsZpos = new TH2F("QvsZpos", "QvsZpos", 1000, -10, 50, 200, -100, 100);
   auto *hexCorr = new TH1F("hexCorr", "hexCorr", 90, -4, 14);

   auto *AngDistr = new TH1F("Ang_Distr", "Ang_Distr", 128, 0, 120);
   auto *AngDistrCM = new TH1F("Ang_Distr_CM", "Ang_Distr_CM", 90, 0, 180);
   auto *ExvsZpos = new TH2F("ExvsZpos", "ExvsZpos", 1000, -5, 15, 200, -20, 150);
   auto *ExvsTrackLength = new TH2F("ExvsTrackLength", "ExvsTrackLength", 1000, -5, 15, 200, -20, 150);
   auto *ExCorrvsZpos = new TH2F("ExCorrvsZpos", "ExCorrvsZpos", 1000, -10, 10, 200, -100, 100);
   auto *KineticEnergy = new TH1F("KineticEnergy", "KineticEnergy", 100, 0, 100);

   /*auto *hredchi2 = new TH1F("redchi2", "redchi2", 1000, 0, 0.0001);
   auto *hbredchi2 = new TH1F("bredchi2", "bredchi2", 1000, 0, 5);*/

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
   Double_t Ebeam_buff = 11.5 * 16; //11.5 
   cout<<" Beam energy in buffer gas : "<<Ebeam_buff<<"\n";
   Double_t m_b = m_d;
   Double_t m_B = m_C15;

   // Ejectile parameters:deuterium
   int A_ej = 2;
   int Z_ej = 1;
   Double_t m_ej = m_d;

   std::vector<TString> filenames;
   // Declare hHex at the beginning of your macro or function
   std::vector<TH1F*> hHex(18);
   // Assuming hHex is already defined and filled
   std::vector<TF1*> fExSpectra_vec(hHex.size(), nullptr);
   

   // Now you can safely initialize and use it
   for (int i = 0; i < hHex.size(); i++) {
      hHex[i] = new TH1F(Form("hex%d", i+1), Form("hex%d", i+1), 90, -4, 14);
   }

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

   filenames.push_back("run_0104_2H.root");
   filenames.push_back("run_0105_2H.root");
   filenames.push_back("run_0106_2H.root");
   filenames.push_back("run_0107_2H.root");
   filenames.push_back("run_0108_2H.root");
   filenames.push_back("run_0109_2H.root");
   filenames.push_back("run_0110_2H.root");
   //filenames.push_back("run_0111_2H.root");
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
   //filenames.push_back("run_0149_2H.root");
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
      TFile *runFile =
      new TFile("/home/georgina/C16_analysis/C16_H2/C16_pd_v16_root/" + filename, "R");
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

         Double_t p_ej = Brho * Z_ej * 2.99792458 / 10 * 1000;
         Double_t E_ej = TMath::Sqrt(p_ej * p_ej + m_ej * m_ej) - m_ej;

         auto [ex_energy, theta_cm] = kine_2b(m_C16, m_p, m_b, m_B, Ebeam_buff, theta, ke);
         
         /*if(ke < 0.0 || ke > 14.0)
                continue;

         if(theta*TMath::RadToDeg() > 60 || theta*TMath::RadToDeg() < 12.0)
                continue;
                */

                
         Double_t Ebeam_at_z = elossH2.GetEnergy(Ebeam_buff, zPos * 100); // 
        
        //Correccion JunRui

         double theta_lab_corr=(theta-(2.0*TMath::Pi()/4000) * (E_ej - kethe)); //theta: rad; theta_lab_corr: rad; E_ej-kethe: MeV
         //cout << "theta lab: " << theta << "  rad  " << " theta_corr_JR: " << theta_lab_corr<< "  rad" << "\n";
         
         auto [ex_energy_corr, theta_cm_corr]= kine_2b(m_C16, m_p, m_b, m_B, Ebeam_at_z, theta_lab_corr, ke);

         // Fill uncorrected histogram
         hex->Fill(ex_energy);      
         ExvsZpos->Fill(ex_energy, zPos*100.0);
         KineticEnergy->Fill(ke);
         
          // Fill corrected histogram
         if (zPos*100> 2.0 && zPos*100 < 60.0)  // only consider reactions occuring within the target region
         {
            //cout << " Ex corrected : " << ex_energy_corr << "  " << " zpos: " << zPos*100.0 << " Ebeam at z: " << Ebeam_at_z << "\n";
            ExCorrvsZpos->Fill(ex_energy_corr, zPos*100.0);
            hexCorr->Fill(ex_energy_corr);
         }

         // Histograms
         //hredchi2->Fill(redchi);

         Ang_Ener->Fill(theta* TMath::RadToDeg(), ke); //theta lab!! -> I still have to implement the correction of catima?
         Ang_Ener_Corr->Fill(theta_lab_corr* TMath::RadToDeg(), ke); //

         Double_t vx = TMath::Sin(theta) * TMath::Sqrt(ke);
         Double_t vy = TMath::Cos(theta) * TMath::Sqrt(ke);

         hVxVy->Fill(vx, vy);

         AngDistr->Fill(theta * TMath::RadToDeg());
         AngDistrCM->Fill(theta_cm);
         hexvstheta->Fill(ex_energy, theta * TMath::RadToDeg());
         ExvsTrackLength->Fill(ex_energy, arclength);

         // Plots of excitation energy in different angular ranges: initialization
         for (int i = 0; i < hHex.size(); i++) {
         double theta_min = 10 + i * 2.5;
         double theta_max = theta_min + 2.5;
         if (theta_cm > theta_min && theta_cm <= theta_max) {
            hHex[i]->Fill(ex_energy_corr);
            break; // Only fill one bin per event
         }
         }
      } // events
   } // Files

   AngDistrCM->Divide(new TF1("sin", "sin(x * TMath::DegToRad())", 0, 180));

//---------------- Fitting the experimental data ----------------//


   TF1 *fExSpectra = new TF1("fExSpectra", "gaus(0) + gaus(3) + [6] * TMath::BreitWigner(x, [7], [8]) + [9] * TMath::BreitWigner(x, [10], [11]) + [12] * TMath::BreitWigner(x, [13], [14])", -5, 14);
    double params[15] = {490, 0.57, 0.24, 
      560, 1.31, 0.218,
      200, 3.9, 0.2,
      50, 5.5, 0.2, 
      25, 6.9, 0.2}; //BW: amplitude mean width
   
   // fExSpectra->SetNpx(10000);
    fExSpectra->SetParameters(params);
    hex->Fit(fExSpectra);

    std::vector<double> p(15);
   for (int i = 0; i < 15; ++i) {
      p[i] = fExSpectra->GetParameter(i);
   }

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
   TString fileKine = Form("/home/georgina/fair_install/ATTPCROOTv2/macro/Kinematics/Decay_kinematics/%s", files[i].c_str());
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
  
TCanvas *c_AngEner = new TCanvas("AngEner", "Energy as a function of #theta", 1200, 800);
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
legend0->Draw();

// Excitation energy spectrum with fits --------------------------------------

   TCanvas *c_ExEner = new TCanvas("ExEner", "Corrected Excited Energy spectra", 1200, 800);
   hexCorr->Draw();
   hexCorr->GetXaxis()->SetTitle("Excitation Energy (MeV)");
   hexCorr->GetYaxis()->SetTitle("Counts");
   
   // Gaussian 1
   TF1 *gaus1 = new TF1("gaus1", "gaus(0)", -5, 14);
   gaus1->SetParameters(p[0], p[1], p[2]);
   gaus1->SetLineColor(kRed);
   gaus1->Draw("same");

   // Gaussian 2
   TF1 *gaus2 = new TF1("gaus2", "gaus(0)", -5, 14);
   gaus2->SetParameters(p[3], p[4], p[5]);
   gaus2->SetLineColor(kBlue);
   gaus2->Draw("same");

   // Breit-Wigner 1
   TF1 *bw1 = new TF1("bw1", "[0]*TMath::BreitWigner(x,[1],[2])", -5, 14);
   bw1->SetParameters(p[6], p[7], p[8]);
   bw1->SetLineColor(kGreen+2);
   bw1->Draw("same");

   // Breit-Wigner 2
   TF1 *bw2 = new TF1("bw2", "[0]*TMath::BreitWigner(x,[1],[2])", -5, 14);
   bw2->SetParameters(p[9], p[10], p[11]);
   bw2->SetLineColor(kMagenta);
   bw2->Draw("same");

   // Breit-Wigner 3
   TF1 *bw3 = new TF1("bw3", "[0]*TMath::BreitWigner(x,[1],[2])", -5, 14);
   bw3->SetParameters(p[12], p[13], p[14]);
   bw3->SetLineColor(kOrange+7);
   bw3->Draw("same");

   TLegend* legend2 = new TLegend(0.7, 0.15, 0.9, 0.3); // (x1, y1, x2, y2) en coordenadas del canvas
   //legend2->SetBorderSize(0); // sin borde
   legend2->SetFillStyle(0);  // fondo transparente
   legend2->AddEntry(hexCorr, "hexCorr", "l");
   legend2->AddEntry(gaus1, "gaus(0)", "l");
   legend2->AddEntry(gaus2, "gaus(1)", "l");
   legend2->AddEntry(bw1, "Breit-Wigner 1", "l");
   legend2->AddEntry(bw2, "Breit-Wigner 2", "l");
   legend2->AddEntry(bw3, "Breit-Wigner 3", "l");
   legend2->Draw();

   c_ExEner->Update();


   TCanvas *c_ExenerCorr = new TCanvas("ExenerCorr", "Excited Energy spectra corrected", 1200, 800);
   c_ExenerCorr->cd();
   c_ExenerCorr->Divide(2, 1);
   c_ExenerCorr->Draw();
   c_ExenerCorr->cd(1);
   hexCorr->Draw();
   c_ExenerCorr->cd(2);
   ExCorrvsZpos->Draw("zcol");

   TH1 *hexClone = (TH1*)hex->Clone("hexClone");
   hexClone->GetListOfFunctions()->Clear(); // This removes the fit line
   hexClone->SetLineColor(kBlue);
   hexCorr->SetLineColor(kRed);

   TCanvas *c_test = new TCanvas("test", "comparison correction hex", 1200, 800);
   c_test->cd();
   hexClone->Draw();
   hexCorr->Draw("same");

   TLegend *legendtest = new TLegend(0.7, 0.7, 0.9, 0.85); // x1, y1, x2, y2 (normalized coordinates)
   legendtest->AddEntry(hexClone, "Original Spectrum", "l");   // "l" for line style
   legendtest->AddEntry(hexCorr, "Corrected Spectrum", "l");
   legendtest->Draw();


   TCanvas *c_AngDistr = new TCanvas("AngDistr", "Angular Distribution", 1200, 600);
   c_AngDistr->cd();
   c_AngDistr->Divide(2, 1);
   c_AngDistr->cd(1);
   AngDistr->Draw();
   c_AngDistr->cd(2);
   AngDistrCM->Draw();
   AngDistrCM->GetXaxis()->SetTitle("Angle (deg)");
   AngDistrCM->GetYaxis()->SetTitle("#frac{d#sigma}{d#Omega} (a.u.)");

   TCanvas *c_ExvsZpos = new TCanvas( "ExvsZpos", "Excitation Energy vs z position and track length", 1200, 600);
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

   /*TCanvas *c_redchi2 = new TCanvas();
   c_redchi2->Divide(2, 1);
   c_redchi2->cd(1);
   hredchi2->Draw("zcol");
   c_redchi2->cd(2);
   hbredchi2->Draw("zcol");*/

   TCanvas *kin = new TCanvas( "kin", "kin", 1200, 600);
    KineticEnergy->Draw("");

   TCanvas *basura = new TCanvas( "basura", "basura", 1200, 600);
    ExvsTrackLength->Draw("zcol");
   
// Angular distribution 

std::vector<std::vector<double>> all_fit_params(hHex.size(), std::vector<double>(15, 0.0));

for (int i = 0; i < hHex.size(); i++) {
    fExSpectra_vec[i] = new TF1(Form("fExSpectra%d", i+1),
        "gaus(0) + gaus(3) + [6]*TMath::BreitWigner(x,[7],[8]) + [9]*TMath::BreitWigner(x,[10],[11]) + [12]*TMath::BreitWigner(x,[13],[14])", -5, 14);
        fExSpectra_vec[i]->SetParameters(params);
        // Set width constraints BEFORE fitting
         /*fExSpectra_vec[i]->SetParLimits(2, 0, 1.0);   // Width of first Gaussian
         fExSpectra_vec[i]->SetParLimits(5, 1.0, 2.0);   // Width of second Gaussian
         fExSpectra_vec[i]->SetParLimits(8, 3.0, 5.0);   // Width of first Breit-Wigner
         fExSpectra_vec[i]->SetParLimits(11, 5.0, 6.5);  // Width of second Breit-Wigner
         fExSpectra_vec[i]->SetParLimits(14, 6.5, 7.2);  // Width of third Breit-Wigner

         */

         /*double params[15] = {490, 0.57, 0.24, 
      560, 1.31, 0.218,
      200, 3.9, 0.2,
      50, 5.5, 0.2, 
      25, 6.9, 0.2}; //BW: amplitude mean width
      */
         
   if (hHex[i]->GetEntries() < 100) {
        fExSpectra_vec[i] = nullptr; // Optional: mark as skipped
        continue; // Skip fitting this histogram
    }
    hHex[i]->Fit(fExSpectra_vec[i]);

    for (int p = 0; p < 15; ++p) {
        all_fit_params[i][p] = fExSpectra_vec[i]->GetParameter(p);
        
    }
}
// Print fit parameters for each angular bin
TCanvas *c_hex_segmented1 = new TCanvas("c_hex_segmented1", "Hex Spectra 1 to 9", 1200, 800);
c_hex_segmented1->cd();
c_hex_segmented1->Divide(3, 3);

TCanvas *c_hex_segmented2 = new TCanvas("c_hex_segmented2", "Hex Spectra 10 to 18", 1200, 800);
c_hex_segmented2->cd();
c_hex_segmented2->Divide(3, 3);

// First canvas
for (int i = 0; i < 9; i++) {
   c_hex_segmented1->cd(i + 1);
   if (hHex[i]) hHex[i]->Draw();
   if (fExSpectra_vec[i]) fExSpectra_vec[i]->Draw("same");
}
// Second canvas

for (int i = 9; i < 18; i++) {
   c_hex_segmented2->cd(i - 8);
   if (hHex[i]) hHex[i]->Draw();
   if (fExSpectra_vec[i]) fExSpectra_vec[i]->Draw("same");
}

c_hex_segmented1->Update();
c_hex_segmented2->Update();

TCanvas *c1 = (TCanvas*)gROOT->GetListOfCanvases()->FindObject("c1");
if (c1) c1->Close();


gROOT->SetSelectedPad(nullptr);
gROOT->SetSelectedPrimitive(nullptr);

c_hex_segmented2->cd();
gPad->Modified();
gPad->Update();
gROOT->SetSelectedPad(nullptr);

if (gROOT->FindObject("basura")) {
    ((TCanvas*)gROOT->FindObject("basura"))->Close();
}

// last_20cm.close();


std::string nombre_pdf = "plots_C16_pd_C15.pdf";

if (guardar_en_pdf) {
    c_AngEner->Print((nombre_pdf + "(").c_str()); // abre el PDF multipágina
    c_AngEner_Corr->Print(nombre_pdf.c_str());
    c_test->Print(nombre_pdf.c_str());
    /*c_ExEner->Print(nombre_pdf.c_str());
    c_ExenerCorr->Print(nombre_pdf.c_str());
    c_AngDistr->Print(nombre_pdf.c_str());
    c_ExvsZpos->Print(nombre_pdf.c_str());
    c_hex_segmented1->Print(nombre_pdf.c_str());
    c_hex_segmented2->Print(nombre_pdf.c_str());*/
    kin->Print((nombre_pdf + ")").c_str()); // cierra el PDF multipágina
  


    gSystem->Exec(("xdg-open " + nombre_pdf).c_str()); // abre el PDF automáticamente
} else {
    c_AngEner->Draw();
    c_AngEner_Corr->Draw();
    c_ExEner->Draw();
    c_ExenerCorr->Draw();
    c_AngDistr->Draw();
    c_ExvsZpos->Draw();
    c_hex_segmented1->Draw();
    c_hex_segmented2->Draw();
}

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