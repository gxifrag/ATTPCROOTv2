#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TF1.h>
#include <TPaveText.h>
#include <iostream>
#include <fstream>


//void analyzeUKF_resolutions(TString fileName = "/home/georgina/fair_install/ATTPCROOTv2_KF_fork/ATTPCROOTv2/macro/tests/UKF/results_recoUKF/reco_ukf_output_pionssim_30MeV_Bfield_20kG_H300torr_theta30_catima_initialMom_10kEvt_hit1_Wrapping.root"){
//void analyzeUKF_resolutions(TString fileName = "/home/georgina/fair_install/ATTPCROOTv2_KF_fork/ATTPCROOTv2/macro/tests/UKF/reco_ukf_output_protonssim_40MeV_Bfield_20kG_H300torr_theta30_catima_initialMom_10kEvt_hit1_Wrapping.root"){
void analyzeUKF_resolutions(TString fileName = "reco_ukf_output_protonssim_40-80MeV_Bfield_20kG_H300torr_theta10-80_catima_initialMom_10kEvt_hit1_cluster10.root"){
//void analyzeUKF_resolutions(TString fileName = "/home/georgina/fair_install/ATTPCROOTv2_KF_fork/ATTPCROOTv2/macro/tests/UKF/reco_ukf_output_pionssim_30MeV_Bfield_20kG_H300torr_theta30_catima_initialMom_10kEvt_hit0_pointsToCluster1.root"){

    TFile *f = new TFile(fileName);
    if (!f || f->IsZombie()) {
        printf("Error: No se pudo abrir el archivo %s\n", fileName.Data());
        return;
    }

    TTree *tree = (TTree*)f->Get("UKFTree"); 
    if (!tree) {
        printf("Error: No se encontro el TTree 'UKFTree'\n");
        return;
    }
    double p_true, p_rec, theta_true, theta_rec, phi_true, phi_rec;
    tree->SetBranchAddress("p_true",     &p_true);
    tree->SetBranchAddress("p_rec",      &p_rec);
    tree->SetBranchAddress("theta_true", &theta_true);// de 0 a pi/2
    tree->SetBranchAddress("theta_rec",  &theta_rec);
    tree->SetBranchAddress("phi_true",   &phi_true); // de -pi a pi
    tree->SetBranchAddress("phi_rec",    &phi_rec);


    TH1F *hP = new TH1F("hP", "Resolution Momentum; (p_{rec} - p_{true})/p_{true} [%]; Counts", 100, -6, 4);
    TH1F *hTheta = new TH1F("hTheta", "Resolution #theta; #theta_{rec} - #theta_{true} [mrad]; Counts", 200, -50, 50);
    TH1F *hPhi   = new TH1F("hPhi",   "Resolution #phi; #phi_{rec} - #phi_{true} [mrad]; Counts", 200, -200, 40);


    TH1F *hP_sep = new TH1F("hP_sep", "Momentum; p_{rec} [MeV/c]; [rad]", 100, 35, 85);
    TH1F *hTheta_sep = new TH1F("hTheta_sep", " #theta_{rec} [rad]; [rad]", 100, -20, 20);
    TH1F *hPhi_sep   = new TH1F("hPhi_sep",   "#phi_{rec} [rad]; [rad]", 100, -4, 4); //-pi to pi

    TH1F *hP_true = new TH1F("hP_true", "Momentum; p_{true} [MeV/c]; [rad]", 100, 35, 85);
    TH1F *hTheta_true = new TH1F("hTheta_true", " #theta_{true} [rad]; [rad]", 100, -5, 5);
    TH1F *hPhi_true   = new TH1F("hPhi_true",   "#phi_{true} [rad]; [rad]", 100, -4, 4); //-pi to pi

    TH2F *hP_corr = new TH2F("hP_corr", "p_{rec} vs. p_{true}; p_{true} [MeV/c]]", 80, 35, 85, 80, 35, 85);
    TH2F *hTheta_corr = new TH2F("hTheta_corr", "#theta_{rec} vs. #theta_{true} [rad]; #theta_{true}", 80, -2, 2, 80, -2, 2);
    TH2F *hPhi_corr   = new TH2F("hPhi_corr",   "#phi_{rec} vs. #phi_{true} [rad]; #phi_{true}", 80, -5, 5, 80, -5, 5);

    /*TH2F *hP_corr = new TH2F("hP_corr", "p_{rec} vs. p_{true}; p_{true} [MeV/c]]", 60, 38, 42, 60, 35, 46);
    TH2F *hTheta_corr = new TH2F("hTheta_corr", "#theta_{rec} vs. #theta_{true} [mrad]; #theta_{true}", 50, -5, 20, 50, -5, 40);
    TH2F *hPhi_corr   = new TH2F("hPhi_corr",   "#phi_{rec} vs. #phi_{true} [mrad]; #phi_{true}", 50, -10, 10, 50, -10, 10);*/

    TH2F *hP_diff_vs_ptrue = new TH2F("hP_diff_vs_ptrue", 
                                     "Momentum Difference vs True Momentum; p_{true} [MeV/c]; p_{rec} - p_{true} [MeV/c]", 
                                     50, 35, 85,   // X-axis bins and range
                                     100, -2.0, 2.0); // Y-axis bins and range

    Long64_t nentries = tree->GetEntries();

    for (Long64_t i=0; i<nentries; i++) {

        tree->GetEntry(i);

        if (p_rec <= 0) continue;
     
        hP->Fill((p_rec - p_true) / p_true * 100.0); //%
        hTheta->Fill((theta_rec - theta_true) *1000); //  mrad
        double dphi = phi_rec - phi_true;
        //std::cout << "DEBUG: phi_rec = " << phi_rec << ", phi_true = " << phi_true << ", dphi (raw) = " << dphi << std::endl;
        double dphi_corr = TVector2::Phi_mpi_pi(dphi);
        if (dphi != dphi_corr) {
            std::cout << "DEBUG: Phi wrapping applied! Original dphi = " << dphi << ", Wrapped dphi = " << dphi_corr << std::endl;
        }
    
        hPhi->Fill(dphi_corr *1000); // mrad
        hP_sep->Fill(p_rec); //rad
        hTheta_sep->Fill(theta_rec); //rad
        hPhi_sep->Fill(phi_rec); //rad

        hP_true->Fill(p_true); //rad
        hTheta_true->Fill(theta_true); //rad
        hPhi_true->Fill(phi_true); //rad

        hP_corr->Fill(p_true, p_rec);//rad
        hTheta_corr->Fill(theta_true, theta_rec); //rad
        hPhi_corr->Fill(phi_true, phi_rec); //rad

        hP_diff_vs_ptrue->Fill(p_true, p_rec - p_true); //MeV/c

    }

    gStyle->SetOptStat(0); 
    gStyle->SetOptFit(1); 

    TCanvas *cRes = new TCanvas("cRes", "UKF Summary - H300 torr", 1800, 600);
    cRes->Divide(4, 1);

    double sigP, sigT, sigPh, meanP, meanT, meanPh;

    // --- Histogram Momentum ---
    cRes->cd(2);
    hP->SetFillColorAlpha(kAzure-9, 0.3);
    hP->SetLineColor(kAzure-9);
    hP->Draw();
    hP->Fit("gaus", "LQ");
    sigP  = hP->GetFunction("gaus")->GetParameter(2);
    meanP = hP->GetFunction("gaus")->GetParameter(1);

    // --- Histogram Theta ---
    cRes->cd(3);
    hTheta->SetFillColorAlpha(kOrange-9, 0.3);
    hTheta->SetLineColor(kOrange-9);
    hTheta->Draw();
    hTheta->Fit("gaus", "LQ");
    sigT  = hTheta->GetFunction("gaus")->GetParameter(2);
    meanT = hTheta->GetFunction("gaus")->GetParameter(1);

    // --- Histogram Phi ---
    cRes->cd(4);
    hPhi->SetFillColorAlpha(kGreen-10, 0.3);
    hPhi->SetLineColor(kGreen-10);
    hPhi->Draw();
    hPhi->Fit("gaus", "LQ");
    sigPh  = hPhi->GetFunction("gaus")->GetParameter(2);
    meanPh = hPhi->GetFunction("gaus")->GetParameter(1);


    // --- RESOLUTIONS SUMMARY BOX ---
    cRes->cd(1);
    TPaveText *pt = new TPaveText(0.1, 0.1, 0.9, 0.9, "NDC"); 
    pt->SetBorderSize(2);
    pt->SetFillColor(kYellow-10); 
    pt->SetTextAlign(22);        
    pt->SetTextFont(42);
    pt->SetTextSize(0.08);      

    pt->AddText("#bf{UKF RESOLUTIONS}");
    pt->AddText("#scale[0.8]{Summary (300 torr)}");
    pt->AddLine(0.1, 0.6, 0.9, 0.6);
    pt->AddText("");
    pt->AddText(Form("#sigma_{p} / p  :  %.3f %%", sigP));
    pt->AddText(Form("#sigma_{#theta}     :  %.3f mrad", sigT));
    pt->AddText(Form("#sigma_{#phi}     :  %.3f mrad", sigPh));
    pt->AddText("");
    pt->AddLine(0.1, 0.3, 0.9, 0.3);
    pt->AddText(Form("#scale[0.7]{Efficiency : %.2f %%}", (double)hP->GetEntries()/10000*100.0));
    pt->Draw();

    cRes->Update();

    //cRes->SaveAs("res_pionssim_30MeV_Bfield_20kG_H300torr_theta30_catima_initialMom_10kEvt_hit0_phiCorrections_PointsToCluster1.png");
    //cRes->SaveAs("test.png");
    cRes->SaveAs("/home/georgina/fair_install/ATTPCROOTv2_KF_fork/ATTPCROOTv2/macro/tests/UKF/results_recoUKF/plots/res_protonssim_40-80MeV_Bfield_20kG_H300torr_theta10-80_catima_initialMom_10kEvt_hit1_cluster10_phiCorrection.png");
    //cRes->SaveAs("res_pionssim_20-40MeV_Bfield_20kG_H300torr_theta0-90.png");

    TCanvas *correlationsCheck = new TCanvas("correlationsCheck", "UKF Correlations Check", 1200, 400);
    correlationsCheck->Divide(3, 2);
    correlationsCheck->cd(1);
    hP_corr->Draw("COLZ");

    correlationsCheck->cd(2);
    hTheta_corr->Draw("COLZ");

    correlationsCheck->cd(3);
    hPhi_corr->Draw("COLZ");

    correlationsCheck->cd(4);
    hP_sep->Draw(); 
    hP_true->SetLineColor(kRed);
    hP_true->SetLineStyle(2);
    hP_true->Draw("SAME");

    correlationsCheck->cd(5);
    hTheta_sep->Draw();
    hTheta_true->SetLineColor(kRed);
    hTheta_true->SetLineStyle(2);
    hTheta_true->Draw("SAME");

    correlationsCheck->cd(6);
    hPhi_sep->Draw();
    hPhi_true->SetLineColor(kRed);
    hPhi_true->SetLineStyle(2);
    hPhi_true->Draw("SAME");

    correlationsCheck->Update();
    //correlationsCheck->SaveAs("correlations_pions_H300torr_20kG_20-40MeV_theta0-90_hit1.png");
    //correlationsCheck->SaveAs("correlations_check.png");
    //correlationsCheck->SaveAs("correlations_pionssim_30MeV_Bfield_20kG_H300torr_theta30_catima_initialMom_10kEvt_hit0_phiCorrections:PointsToCLuster1.png");

    correlationsCheck->SaveAs("/home/georgina/fair_install/ATTPCROOTv2_KF_fork/ATTPCROOTv2/macro/tests/UKF/results_recoUKF/plots/correlations_reco_ukf_output_protonssim_40-80MeV_Bfield_20kG_H300torr_theta10-80_catima_initialMom_10kEvt_hit1_cluster10_phiCorrection.png");
   
   /*TCanvas *cMomDiff = new TCanvas("cMomDiff", "Momentum Difference vs True Momentum", 800, 600);
    cMomDiff->cd();
    
    // Draw with COLZ for the color palette
    hP_diff_vs_ptrue->Draw("COLZ"); 
    
    // Draw a dashed horizontal red line at Y = 0 to guide the eye
    TLine *line0 = new TLine(35, 0, 85, 0); 
    line0->SetLineStyle(2); // Dashed
    line0->SetLineColor(kRed);
    line0->SetLineWidth(2);
    line0->Draw("SAME");

    cMomDiff->Update();
    cMomDiff->SaveAs("mom_diff_vs_ptrue.png");*/
   
    //ofstream outfile("resolutions_H300torr_40MeV_theta30.txt");
    //outfile << "UKF ANALYSIS SUMMARY - File: " << fileName << endl;
    //outfile << "-----------------------------------------------" << endl;
    //outfile << "MOMENTUM (p):   Sigma = " << sigP << " %   | Bias = " << meanP << " %" << endl;
    //outfile << "THETA (theta):  Sigma = " << sigT << " mrad | Bias = " << meanT << " mrad" << endl;
    //outfile << "PHI (phi):      Sigma = " << sigPh << " mrad | Bias = " << meanPh << " mrad" << endl;
    //outfile << "-----------------------------------------------" << endl;
    //outfile << "Reconstruction Efficiency: " << (double)hP->GetEntries()/nentries*100.0 << " %" << endl;
    //outfile.close();

    //cout << "\n>>> Results saved in 'resolutions_H300torr_40MeV_theta30.txt'" << endl;
}