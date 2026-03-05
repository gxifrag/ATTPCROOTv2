#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TF1.h>
#include <TPaveText.h>
#include <iostream>
#include <fstream>

double wrapAngle(double a){
    return TMath::ATan2(TMath::Sin(a), TMath::Cos(a));
}

std::vector<double> wrapAll(const std::vector<double>& delta){
    std::vector<double> out; out.reserve(delta.size());
    for(double d: delta) out.push_back(wrapAngle(d));
    return out;
}

td::pair<double,double> circularStats(const std::vector<double>& delta){
    double sx=0, cx=0;
    for(double d: delta){ sx += sin(d); cx += cos(d); }
    sx /= delta.size(); cx /= delta.size();
    double mean = atan2(sx, cx);                 // media circular
    double r = sqrt(sx*sx + cx*cx);              // resultant length
    double circ_std = sqrt(-2.0 * log(r));       // aproximación
    return {mean, circ_std};
}

//void analyzeUKF_resolutions(TString fileName = "reco_ukf_output_3T_H300torr_p40MeV_theta30.root") {
void analyzeUKF_resolutions(TString fileName = "reco_ukf_output_3T_H300torr_p40MeV_theta30.root") {

    // 1. Abrir el archivo y obtener el árbol
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

    // 2. Variables para las ramas
    double p_true, p_rec, theta_true, theta_rec, phi_true, phi_rec;
    tree->SetBranchAddress("p_true",     &p_true);
    tree->SetBranchAddress("p_rec",      &p_rec);
    tree->SetBranchAddress("theta_true", &theta_true);
    tree->SetBranchAddress("theta_rec",  &theta_rec);
    tree->SetBranchAddress("phi_true",   &phi_true);
    tree->SetBranchAddress("phi_rec",    &phi_rec);

    // 3. Definición de Histogramas (Rangos ajustados según tus datos de 300torr)
    TH1F *hP = new TH1F("hP", "Resolucion Momento; (p_{rec} - p_{true})/p_{true} [%]; Counts", 100, -5, 5);
    TH1F *hTheta = new TH1F("hTheta", "Resolucion #theta; #theta_{rec} - #theta_{true} [mrad]; Counts", 100, -100, 100);
    TH1F *hPhi   = new TH1F("hPhi",   "Resolucion #phi; #phi_{rec} - #phi_{true} [mrad]; Counts", 100, -50, 150);

    // 4. Llenar los histogramas
    Long64_t nentries = tree->GetEntries();
    for (Long64_t i=0; i<nentries; i++) {
        tree->GetEntry(i);
        if (p_rec > 0) {
            hP->Fill((p_rec - p_true) / p_true * 100.0);
            hTheta->Fill((theta_rec - theta_true)*1000*TMath::DegToRad());
            hPhi->Fill((phi_rec - phi_true)*1000*TMath::DegToRad());
        }
    }

    // 5. Estilo
    gStyle->SetOptStat(0); // Quitamos el stat box de los histogramas para que no tapen
    gStyle->SetOptFit(1);  // Pero dejamos el del fit si lo deseas (aparecerá arriba a la derecha)

    TCanvas *cRes = new TCanvas("cRes", "UKF Summary - H300 torr", 1800, 600);
    cRes->Divide(4, 1);

    // Procesamos cada panel y guardamos los resultados del Fit para el TXT y el cuadro
    double sigP, sigT, sigPh, meanP, meanT, meanPh;

    // --- Histograma Momento ---
    cRes->cd(2);
    hP->SetFillColorAlpha(kAzure-9, 0.3);
    hP->SetLineColor(kAzure-9);
    hP->Draw();
    hP->Fit("gaus", "LQ"); 
    sigP  = hP->GetFunction("gaus")->GetParameter(2);
    meanP = hP->GetFunction("gaus")->GetParameter(1);

    // --- Histograma Theta ---
    cRes->cd(3);
    hTheta->SetFillColorAlpha(kOrange-9, 0.3);
    hTheta->SetLineColor(kOrange-9);
    hTheta->Draw();
    hTheta->Fit("gaus", "LQ");
    sigT  = hTheta->GetFunction("gaus")->GetParameter(2);
    meanT = hTheta->GetFunction("gaus")->GetParameter(1);

    // --- Histograma Phi ---
    cRes->cd(4);
    hPhi->SetFillColorAlpha(kGreen-10, 0.3);
    hPhi->SetLineColor(kGreen-10);
    hPhi->Draw();
    hPhi->Fit("gaus", "LQ");
    sigPh  = hPhi->GetFunction("gaus")->GetParameter(2);
    meanPh = hPhi->GetFunction("gaus")->GetParameter(1);

    // --- Panel 1: CUADRO DE RESUMEN TIPO PÓSTER ---
    cRes->cd(1);
    TPaveText *pt = new TPaveText(0.1, 0.1, 0.9, 0.9, "NDC"); // Grande y centrado
    pt->SetBorderSize(2);
    pt->SetFillColor(kYellow-10); // Color crema/suave
    pt->SetTextAlign(22);        // Centrado total
    pt->SetTextFont(42);
    pt->SetTextSize(0.08);       // Texto bien grande

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

    // 6. Guardar imagen y TXT
    cRes->SaveAs("UKF_Summary_H300torr_40MeV_theta30.png");


    //coomprobaciones:

    double ssum=0, csum=0;
    for(int i=0;i<N;i++){
        ssum += w[i] * sin(sigmaPhi[i]);
        csum += w[i] * cos(sigmaPhi[i]);
    }
    double phi_mean = atan2(ssum, csum);

    cout << "\n>>> Circular Stats for Phi Residuals:" << endl;
    cout << "    - Mean (wrapped): " << phi_mean*1000 << " mrad" << endl;
    cout << "    - Sigma (approx): " << sqrt(-2.0 * log(sqrt(ssum*ssum + csum*csum))) * 1000 << " mrad" << endl;    

    ofstream outfile("resoluciones_H300torr_40MeV_theta30.txt");
    outfile << "RESUMEN ANALISIS UKF - Archivo: " << fileName << endl;
    outfile << "-----------------------------------------------" << endl;
    outfile << "MOMENTO (p):   Sigma = " << sigP << " %   | Bias = " << meanP << " %" << endl;
    outfile << "THETA (theta): Sigma = " << sigT << " mrad | Bias = " << meanT << " mrad" << endl;
    outfile << "PHI (phi):     Sigma = " << sigPh << " mrad | Bias = " << meanPh << " mrad" << endl;
    outfile << "-----------------------------------------------" << endl;
    outfile << "Eficiencia Reconstruccion: " << (double)hP->GetEntries()/nentries*100.0 << " %" << endl;
    outfile.close();

    cout << "\n>>> Resultados guardados en 'resoluciones_H300torr_40MeV_theta30.txt'" << endl;
}