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

double wrapDeg(double a){
    double r = fmod(a + 180.0, 360.0);
    if (r < 0) r += 360.0;
    return r - 180.0;
}

//void analyzeUKF_resolutions(TString fileName = "reco_ukf_output_3T_H300torr_p40MeV_theta30.root") {
void analyzeUKF_resolutions(TString fileName = "reco_ukf_output_hit1.root") {

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
    tree->SetBranchAddress("theta_true", &theta_true);// de -pi a pi
    tree->SetBranchAddress("theta_rec",  &theta_rec);
    tree->SetBranchAddress("phi_true",   &phi_true);
    tree->SetBranchAddress("phi_rec",    &phi_rec);

    // 3. Definición de Histogramas (Rangos ajustados según tus datos de 300torr)
    TH1F *hP = new TH1F("hP", "Resolucion Momento; (p_{rec} - p_{true})/p_{true} [%]; Counts", 160, -4, 4);
    TH1F *hTheta = new TH1F("hTheta", "Resolucion #theta; #theta_{rec} - #theta_{true} [mrad]; Counts", 160, -40, 40);
    TH1F *hPhi   = new TH1F("hPhi",   "Resolucion #phi; #phi_{rec} - #phi_{true} [mrad]; Counts", 160, -40, 40);

    // 4. Llenar los histogramas
    Long64_t nentries = tree->GetEntries();

    for (Long64_t i=0; i<nentries; i++) {
        tree->GetEntry(i);
        if (p_rec <= 0) continue; // Solo consideramos eventos con p_rec > 0 para evitar problemas de división y asegurar que son eventos reconstruidos
        
        hP->Fill((p_rec - p_true) / p_true * 100.0);
        double dtheta_mrad = (theta_rec - theta_true) * TMath::DegToRad() * 1000; // Convertimos a mrad
        hTheta->Fill(dtheta_mrad);

        //printf("entry %lld: p_true=%g  p_rec=%g  theta_true=%g  theta_rec=%g  dtheta_mrad=%g\n",
                //i, p_true, p_rec, theta_true, theta_rec, dtheta_mrad);
            
            hPhi->Fill((phi_rec - phi_true)*1000*TMath::DegToRad());
            //hPhi->Fill(phi_rec);

        double phi_t = wrapDeg(phi_true);
        double phi_r = wrapDeg(phi_rec);
        double ddeg = wrapDeg(phi_r - phi_t);

        double dphi_mrad = ddeg * TMath::DegToRad()*1000 ;              // a rad
        //printf("entry %lld: phi_true=%g  phi_t=%g  phi_rec=%g  phi_rec_wrapped=%g  dphi_deg=%g  dphi_mrad=%g\n",i, phi_true, phi_t, phi_rec, phi_r, ddeg, dphi_mrad);
        //hPhi->Fill(dphi_mrad);                 // mrad
        //deltasPhi_rad.push_back(dphi_rad);
                
        
    }

   /* for (int i=0;i<nentries;i++){
        tree->GetEntry(i);
         double phi_r_wrapped = wrapDeg(phi_rec);
       // comprobar phi_true raw (si esperas que siempre esté en [-180,180])
        if (phi_true > 180.0 || phi_true < -180.0) {
            double diff = phi_rec - phi_r_wrapped; // diferencia raw - wrapped (grados)
            printf("WARNING: Angulo phi_true fuera de rango en entry %d: phi_true=%g  phi_rec_raw=%g  phi_rec_wrapped=%g  diff_raw_wrapped=%g\n",
                i, phi_true, phi_rec, phi_r_wrapped, diff);
        }

        // comprobar el valor normalizado (phi_r_wrapped), no el raw
        if (phi_rec > 180.0 || phi_rec< -180.0) {
             double diff = phi_rec - phi_r_wrapped; 
            printf("WARNING: Angulo phi_rec_wrapped fuera de rango en entry %d: phi_true=%g  phi_rec_raw=%g  phi_rec_wrapped=%g  diff_raw_wrapped=%g\n",
                i, phi_true, phi_rec, phi_r_wrapped, diff);
        }
                    //printf("entry %d: phi_true=%g  phi_rec=%g  theta_true=%g  theta_rec=%g\n",
                //i, phi_true, phi_rec, theta_true, theta_rec);
    }*/

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