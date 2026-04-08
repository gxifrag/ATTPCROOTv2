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
//void analyzeUKF_resolutions(TString fileName = "reco_ukf_output_protonssim_40-80MeV_Bfield_20kG_H300torr_theta10-80_catima_initialMom_10kEvt_hit1_cluster10.root"){
void analyzeUKF_resolutions(TString fileName = "8April-reco_ukf_output_pionssim_20-40MeV_Bfield_20kG_H300torr_theta0-90_catima_initialMom_10kEvt_hit1_withStragg_mediaDensityH300_change.root"){

    TFile *f = new TFile(fileName);
    if (!f || f->IsZombie()) {
        printf("Error: Cannot open file %s\n", fileName.Data());
        return;
    }

    TTree *tree = (TTree*)f->Get("UKFTree"); 
    if (!tree) {
        printf("Error: Cannot find TTree 'UKFTree'\n");
        return;
    }
    double p_true, p_rec, theta_true, theta_rec, phi_true, phi_rec, res_theta, res_phi;
    std::vector<double> *hit_res_smooth = nullptr;
    std::vector<double> *hit_res_smooth_x = nullptr;
    std::vector<double> *hit_res_smooth_y = nullptr;
    std::vector<double> *hit_res_smooth_z = nullptr;

    tree->SetBranchAddress("p_true",     &p_true);
    tree->SetBranchAddress("p_rec",      &p_rec);
    tree->SetBranchAddress("theta_true", &theta_true);// de 0 a pi/2
    tree->SetBranchAddress("theta_rec",  &theta_rec);
    tree->SetBranchAddress("phi_true",   &phi_true); // de -pi a pi
    tree->SetBranchAddress("phi_rec",    &phi_rec);
    tree->SetBranchAddress("res_theta",  &res_theta);
    tree->SetBranchAddress("res_phi",    &res_phi);

    tree->SetBranchAddress("hit_res_smooth", &hit_res_smooth);
    tree->SetBranchAddress("hit_res_smooth_x", &hit_res_smooth_x);
    tree->SetBranchAddress("hit_res_smooth_y", &hit_res_smooth_y);
    tree->SetBranchAddress("hit_res_smooth_z", &hit_res_smooth_z);

    TH1F *hP = new TH1F("hP", "Resolution Momentum; (p_{rec} - p_{true})/p_{true} [%]; Counts", 80, -2, 2);
    TH1F *hTheta = new TH1F("hTheta", "Resolution #theta; #theta_{rec} - #theta_{true} [mrad]; Counts", 100, -10, 10);
    TH1F *hPhi   = new TH1F("hPhi",   "Resolution #phi; #phi_{rec} - #phi_{true} [mrad]; Counts", 100, -20, 20);

    TH1F *hP_sep = new TH1F("hP_sep", "Momentum; p_{rec} [MeV/c]; [rad]", 100, 25, 45);
    TH1F *hTheta_sep = new TH1F("hTheta_sep", " #theta_{rec} [rad]; [rad]", 100, -20, 20);
    TH1F *hPhi_sep   = new TH1F("hPhi_sep",   "#phi_{rec} [rad]; [rad]", 100, -4, 4); //-pi to pi

    TH1F *hP_true = new TH1F("hP_true", "Momentum; p_{true} [MeV/c]; [rad]", 100, 25, 45);
    TH1F *hTheta_true = new TH1F("hTheta_true", " #theta_{true} [rad]; [rad]", 100, -5, 5);
    TH1F *hPhi_true   = new TH1F("hPhi_true",   "#phi_{true} [rad]; [rad]", 100, -4, 4); //-pi to pi

    TH2F *hP_corr = new TH2F("hP_corr", "p_{rec} vs. p_{true}; p_{true} [MeV/c]]", 80, 25, 45, 80, 25, 45);
    TH2F *hTheta_corr = new TH2F("hTheta_corr", "#theta_{rec} vs. #theta_{true} [rad]; #theta_{true}", 80, -2, 2, 80, -2, 2);
    TH2F *hPhi_corr   = new TH2F("hPhi_corr",   "#phi_{rec} vs. #phi_{true} [rad]; #phi_{true}", 80, -5, 5, 80, -5, 5);

    /*TH2F *hP_corr = new TH2F("hP_corr", "p_{rec} vs. p_{true}; p_{true} [MeV/c]]", 60, 38, 42, 60, 35, 46);
    TH2F *hTheta_corr = new TH2F("hTheta_corr", "#theta_{rec} vs. #theta_{true} [mrad]; #theta_{true}", 50, -5, 20, 50, -5, 40);
    TH2F *hPhi_corr   = new TH2F("hPhi_corr",   "#phi_{rec} vs. #phi_{true} [mrad]; #phi_{true}", 50, -10, 10, 50, -10, 10);*/

    TH1F* hRes3D = new TH1F("hRes3D", "Smoother 3D Spatial Residuals;Distance to Geant4 Cluster [mm];Counts", 100, 0.0, 0.5);
    TH1F* hRes3D_x = new TH1F("hRes3D_x", "Smoother Residuals in X;Residual in X [mm];Counts", 100, -0.5, 0.5);
    TH1F* hRes3D_y = new TH1F("hRes3D_y", "Smoother Residuals in Y;Residual in Y [mm];Counts", 100, -0.5, 0.5);
    TH1F* hRes3D_z = new TH1F("hRes3D_z", "Smoother Residuals in Z;Residual in Z [mm];Counts", 100, -0.5, 0.5);

    Long64_t nentries = tree->GetEntries();

    for (Long64_t i=0; i<nentries; i++) {

        tree->GetEntry(i);

        if (p_rec <= 0) continue;
     
        hP->Fill((p_rec - p_true) / p_true * 100.0); //%  
        hTheta->Fill(res_theta *1000); // mrad
        hPhi->Fill(res_phi *1000); // mrad
        
        hP_sep->Fill(p_rec); //rad
        hTheta_sep->Fill(theta_rec); //rad
        hPhi_sep->Fill(phi_rec); //rad

        hP_true->Fill(p_true); //rad
        hTheta_true->Fill(theta_true); //rad
        hPhi_true->Fill(phi_true); //rad

        hP_corr->Fill(p_true, p_rec);//rad
        hTheta_corr->Fill(theta_true, theta_rec); //rad
        hPhi_corr->Fill(phi_true, phi_rec); //rad

        for (size_t p = 0; p < hit_res_smooth->size(); p++) {

            hRes3D->Fill(hit_res_smooth->at(p));
            hRes3D_x->Fill(hit_res_smooth_x->at(p));
            hRes3D_y->Fill(hit_res_smooth_y->at(p));
            hRes3D_z->Fill(hit_res_smooth_z->at(p));
        }

    }

    //gStyle->SetOptStat(0); 
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
    hTheta->Fit("gaus", "LQR", "", 0, 20); // Fit only within a reasonable range to avoid tails
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
    //cRes->SaveAs("/home/georgina/fair_install/ATTPCROOTv2_KF_fork/ATTPCROOTv2/macro/tests/UKF/results_recoUKF/plots/res_reco_ukf_output_pionssim_30MeV_Bfield_20kG_H300torr_theta30_catima_initialMom_10kEvt_hit1_withStragg_mediaDensityH300_change_covMatrixUnit.png");
    //cRes->SaveAs("res_pionssim_20-40MeV_Bfield_20kG_H300torr_theta0-90.png");

    TCanvas *correlationsCheck = new TCanvas("correlationsCheck", "UKF Correlations Check", 1800, 800);
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

   // correlationsCheck->SaveAs("/home/georgina/fair_install/ATTPCROOTv2_KF_fork/ATTPCROOTv2/macro/tests/UKF/results_recoUKF/plots/correlations_reco_ukf_output_pionssim_30MeV_Bfield_20kG_H300torr_theta30_catima_initialMom_10kEvt_hit1_withStragg_mediaDensityH300_change_covMatrixUnit.png");
   

    TCanvas* cRes3D = new TCanvas("cRes3D", "Residuals", 800, 600);
    hRes3D->Draw();
    hRes3D_x->SetLineColor(kRed);
    hRes3D_x->SetLineStyle(2);
    hRes3D_x->Draw("SAME");
    hRes3D_y->SetLineColor(kGreen+2);
    hRes3D_y->SetLineStyle(2);      
    hRes3D_y->Draw("SAME");
    hRes3D_z->SetLineColor(kViolet+2);
    hRes3D_z->SetLineStyle(2);
    hRes3D_z->Draw("SAME");
    //cRes3D->SaveAs("/home/georgina/fair_install/ATTPCROOTv2_KF_fork/ATTPCROOTv2/macro/tests/UKF/results_recoUKF/plots/residuals3D_reco_ukf_output_pionssim_30MeV_Bfield_20kG_H300torr_theta30_catima_initialMom_10kEvt_hit1_withStragg_mediaDensityH300_change_covMatrixUnit.png");


    TCut cut_ok = "status == 0";
    TCut cut_clean_theta = "abs(res_theta*1000) < 50"; 
    TCut cut_clean_phi   = "abs(res_phi*1000) < 100";  
    TCut cut_clean_p     = "abs((p_rec - p_true)/p_true * 100) < 10"; // Filtramos errores de momento > 10%

    TCut final_cut_theta = cut_ok && cut_clean_theta;
    TCut final_cut_phi   = cut_ok && cut_clean_phi;
    TCut final_cut_p     = cut_ok && cut_clean_p;

    TCanvas* c1 = new TCanvas("c1", "Resolucion Theta", 1200, 600);
    c1->Divide(2, 1);

    c1->cd(1);
    tree->Draw("res_theta*1000 >> hResTheta(100, -20, 20)", final_cut_theta);
    TH1D* hResTheta = (TH1D*)gDirectory->Get("hResTheta");
    hResTheta->SetTitle("Resolucion Polar (#theta);#theta_{rec} - #theta_{true} [mrad];Eventos");
    hResTheta->SetLineColor(kBlue+1);
    hResTheta->SetFillColor(kBlue-9);

    c1->cd(2);
    // Profile: Sesgo de theta en función del ángulo real
    tree->Draw("res_theta*1000 : theta_true >> hProfTheta(50, 0, 1.6)", final_cut_theta, "prof");
    TProfile* hProfTheta = (TProfile*)gDirectory->Get("hProfTheta");
    hProfTheta->SetTitle("Sesgo de #theta vs Angulo de disparo;#theta_{true} [rad];Sesgo medio [mrad]");
    hProfTheta->SetMarkerStyle(20);
    hProfTheta->SetMarkerColor(kRed);

    // ========================================================================
    // CANVAS 2: RESOLUCIÓN ANGULAR (PHI - Plano XY)
    // ========================================================================
    TCanvas* c2 = new TCanvas("c2", "Resolucion Phi", 1200, 600);
    c2->Divide(2, 1);

    c2->cd(1);
    tree->Draw("res_phi*1000 >> hResPhi(100, -50, 50)", final_cut_phi);
    TH1D* hResPhi = (TH1D*)gDirectory->Get("hResPhi");
    hResPhi->SetTitle("Resolucion Azimutal (#phi);#phi_{rec} - #phi_{true} [mrad];Eventos");
    hResPhi->SetLineColor(kGreen+2);
    hResPhi->SetFillColor(kGreen-9);

    c2->cd(2);
    tree->Draw("res_phi*1000 : theta_true >> hProfPhi(50, 0, 1.6)", final_cut_phi, "prof");
    TProfile* hProfPhi = (TProfile*)gDirectory->Get("hProfPhi");
    hProfPhi->SetTitle("Sesgo de #phi vs Angulo Polar;#theta_{true} [rad];Sesgo medio de #phi [mrad]");
    hProfPhi->SetMarkerStyle(20);
    hProfPhi->SetMarkerColor(kMagenta);

    // ========================================================================
    // CANVAS 3: RESOLUCIÓN DE MOMENTO (P)
    // ========================================================================
    TCanvas* c3 = new TCanvas("c3", "Resolucion Momento", 1200, 600);
    c3->Divide(2, 1);

    c3->cd(1);
    // Error relativo porcentual: (P_rec - P_true) / P_true * 100
    tree->Draw("(p_rec - p_true)/p_true * 100 >> hResP(100, -2, 2)", final_cut_p);
    TH1D* hResP = (TH1D*)gDirectory->Get("hResP");
    hResP->SetTitle("Resolucion Relativa de Momento;(p_{rec} - p_{true})/p_{true} [%];Eventos");
    hResP->SetLineColor(kRed+1);
    hResP->SetFillColor(kRed-9);

    c3->cd(2);
    // Vemos si la resolución de momento empeora a ciertos ángulos
    tree->Draw("(p_rec - p_true)/p_true * 100 : theta_true >> hProfP(50, 0, 1.6)", final_cut_p, "prof");
    TProfile* hProfP = (TProfile*)gDirectory->Get("hProfP");
    hProfP->SetTitle("Sesgo de Momento vs Angulo Polar;#theta_{true} [rad];Sesgo de p [%]");
    hProfP->SetMarkerStyle(20);
    hProfP->SetMarkerColor(kBlue);

    // ========================================================================
    // CANVAS 4: DIAGNÓSTICO DE PÉRDIDA DE ENERGÍA (ELOSS)
    // ========================================================================
    TCanvas* c4 = new TCanvas("c4", "Diagnostico Eloss", 1200, 600);
    c4->Divide(2, 1);

    c4->cd(1);
    // Dibujamos todos los puntos de los vectores (Fluctuaciones vs BetheBloch medio)
    tree->SetMarkerColor(kBlue);
    tree->Draw("sim_eloss_MC", cut_ok);
    tree->SetMarkerColor(kRed);
    tree->Draw("rec_eloss_Smooth", cut_ok, "SAME");
    // Truco para añadir un título rápido al gráfico generado automáticamente
    if (gPad->GetPrimitive("htemp")) {
        ((TH1F*)gPad->GetPrimitive("htemp"))->SetTitle("MC (Azul) vs Smoother (Rojo);Energia Perdida [MeV]");
    }

    c4->cd(2);
    // Distribución del error en la estimación de energía por hit
    tree->Draw("(rec_eloss_Smooth - sim_eloss_MC)*1e6 >> hDiffEloss(100, -5000, 5000)", cut_ok);
    TH1D* hDiffEloss = (TH1D*)gDirectory->Get("hDiffEloss");
    hDiffEloss->SetTitle("Error de dE/dx por hit (Smoother - MC);Diferencia [eV];Hits");
    hDiffEloss->SetLineColor(kOrange+7);
    hDiffEloss->SetFillColor(kOrange-9);

    // Actualizar todos los lienzos
    c1->Update();
    c2->Update();
    c3->Update();
    c4->Update();

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