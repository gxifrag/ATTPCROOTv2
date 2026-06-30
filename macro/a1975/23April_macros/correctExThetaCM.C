// correctExThetaCM.C
//
// Objetivo: cuantificar y corregir el pequeño tilt de Ex con theta_CM
// para el estado fundamental de 15C (banda ~0-1 MeV).
//
// Uso:
//   root -l
//   .L correctExThetaCM.C+
//   correctExThetaCM("tu_archivo.root", "h2_Ex_vs_thetaCM")
//
// h2 debe tener Ex en X y theta_CM en Y (como en tu plot: Ex eje X, theta_CM eje Y)

#include "TFile.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <vector>
#include <iostream>

void correctExThetaCM(const char* filename, const char* histname,
                       double theta_min = 25., double theta_max = 55.,
                       double dtheta = 5.,
                       double fit_lo = -0.5, double fit_hi = 0.5)
{
    TFile* f = TFile::Open(filename, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "No puedo abrir " << filename << std::endl;
        return;
    }

    TH2D* h2 = (TH2D*)f->Get(histname);
    if (!h2) {
        std::cerr << "No encuentro el histograma " << histname << std::endl;
        return;
    }

    std::vector<double> vTheta, vThetaErr, vEx, vExErr;

    TCanvas* cSlices = new TCanvas("cSlices", "Slices en theta_CM", 1400, 900);
    int nSlices = (int)((theta_max - theta_min) / dtheta);
    int nx = (int)ceil(sqrt((double)nSlices));
    cSlices->Divide(nx, nx);

    int padIdx = 1;
    for (double th = theta_min; th < theta_max; th += dtheta) {
        double thLo = th, thHi = th + dtheta;

        int binYlo = h2->GetYaxis()->FindBin(thLo);
        int binYhi = h2->GetYaxis()->FindBin(thHi) - 1;
        if (binYhi < binYlo) continue;

        TH1D* proj = h2->ProjectionX(Form("px_%d", padIdx), binYlo, binYhi);
        if (proj->GetEntries() < 30) {
            // poca estadistica en esta franja, la saltamos
            continue;
        }

        // Hay dos picos muy juntos (g.s. y un estado vecino), asi que usamos
        // un DOBLE gaussiano y nos quedamos solo con el componente del g.s.
        double gs_seed = proj->GetBinCenter(proj->GetMaximumBin());
        double nb_seed = gs_seed + 0.8; // pico vecino, tipicamente ~0.8-1 MeV mas arriba

        TF1* fDouble = new TF1(Form("fDouble_%d", padIdx),
                                "gaus(0)+gaus(3)", fit_lo, fit_hi);
        fDouble->SetParameters(proj->GetMaximum(), gs_seed, 0.3,
                                proj->GetMaximum()*0.5, nb_seed, 0.3);
        fDouble->SetParLimits(1, gs_seed - 0.4, gs_seed + 0.4);
        fDouble->SetParLimits(4, nb_seed - 0.4, nb_seed + 0.6);
        fDouble->SetParLimits(2, 0.05, 0.7);
        fDouble->SetParLimits(5, 0.05, 0.7);
        proj->Fit(fDouble, "RQ0"); // Q=quiet, 0=no draw aun, R=usa rango

        double mean = fDouble->GetParameter(1);
        double meanErr = fDouble->GetParError(1);

        // Gaussiano individual del g.s. para dibujar superpuesto (solo visual)
        TF1* fGaus = new TF1(Form("fGaus_%d", padIdx), "gaus", fit_lo, fit_hi);
        fGaus->SetParameters(fDouble->GetParameter(0), fDouble->GetParameter(1), fDouble->GetParameter(2));
        fGaus->SetLineColor(kGreen+2);

        // theta_CM representativo del slice (centro) y su "anchura" como error en X
        double thetaMid = 0.5 * (thLo + thHi);
        double thetaErr = 0.5 * dtheta;

        vTheta.push_back(thetaMid);
        vThetaErr.push_back(thetaErr);
        vEx.push_back(mean);
        vExErr.push_back(meanErr);

        cSlices->cd(padIdx);
        proj->GetXaxis()->SetRangeUser(fit_lo - 1, fit_hi + 1);
        proj->SetTitle(Form("#theta_{CM} #in [%.0f, %.0f)", thLo, thHi));
        proj->Draw();
        fGaus->Draw("same");
        padIdx++;
    }

    if (vTheta.size() < 3) {
        std::cerr << "Muy pocos slices con estadistica suficiente. "
                   << "Prueba aumentando dtheta." << std::endl;
        return;
    }

    // --- TGraphErrors centroide vs theta_CM ---
    TGraphErrors* g = new TGraphErrors(vTheta.size(),
                                        &vTheta[0], &vEx[0],
                                        &vThetaErr[0], &vExErr[0]);
    g->SetTitle("Centroide Ex (g.s.) vs #theta_{CM};#theta_{CM} (deg);E_{x} centroid (MeV)");
    g->SetMarkerStyle(20);
    g->SetMarkerColor(kBlue+2);
    g->SetLineColor(kBlue+2);

    TF1* fLin = new TF1("fLin", "pol1", theta_min, theta_max);
    g->Fit(fLin, "Q");

    double beta = fLin->GetParameter(1);   // pendiente, MeV/deg
    double betaErr = fLin->GetParError(1);
    double offset = fLin->GetParameter(0);

    std::cout << "=====================================" << std::endl;
    std::cout << " Fit lineal: Ex(theta) = a + b*theta " << std::endl;
    std::cout << " a (offset) = " << offset << " MeV" << std::endl;
    std::cout << " b (beta)   = " << beta << " +/- " << betaErr << " MeV/deg" << std::endl;
    std::cout << " chi2/ndf   = " << fLin->GetChisquare() << " / " << fLin->GetNDF() << std::endl;
    std::cout << "=====================================" << std::endl;
    std::cout << " Para corregir, usa (con theta_ref a tu eleccion, ej. centro estadistico):" << std::endl;
    std::cout << "   Ex_corr = Ex_meas - beta * (theta_CM - theta_ref);" << std::endl;
    std::cout << " con beta = " << beta << " MeV/deg" << std::endl;
    std::cout << "=====================================" << std::endl;

    TCanvas* cFit = new TCanvas("cFit", "Tilt Ex vs theta_CM", 800, 600);
    g->Draw("AP");
    fLin->SetLineColor(kRed);
    fLin->Draw("same");

    // --- Histograma 2D corregido para verificacion visual ---
    double theta_ref = 0.5 * (theta_min + theta_max); // puedes cambiarlo
    TH2D* h2corr = (TH2D*)h2->Clone("h2_corrected");
    h2corr->Reset();
    h2corr->SetTitle("Ex corregido vs theta_CM");

    int nbinsX = h2->GetNbinsX();
    int nbinsY = h2->GetNbinsY();
    int nFilled = 0;
    for (int by = 1; by <= nbinsY; by++) {
        const double theta_by = h2->GetYaxis()->GetBinCenter(by); // distinto nombre que theta_ref a propósito
        const double shift = beta * (theta_by - theta_ref);
        for (int bx = 1; bx <= nbinsX; bx++) {
            double w = h2->GetBinContent(bx, by);
            if (w <= 0) continue;
            double ex = h2->GetXaxis()->GetBinCenter(bx);
            double exCorr = ex - shift;
            // usar theta_by (el theta real de este bin), NUNCA theta_ref
            h2corr->Fill(exCorr, theta_by, w);
            nFilled++;
        }
    }
    std::cout << "Bins con contenido procesados en el relleno: " << nFilled << std::endl;
    if (nFilled == 0) {
        std::cerr << "AVISO: h2corr no se ha rellenado con nada, revisa el histograma de entrada." << std::endl;
    }

    TCanvas* cComp = new TCanvas("cComp", "Antes vs despues", 1400, 700);
    cComp->Divide(2, 1);
    cComp->cd(1);
    h2->SetTitle("Original");
    h2->Draw("colz");
    cComp->cd(2);
    h2corr->SetTitle("Corregido");
    h2corr->Draw("colz");

    // Guardar resultados
    TFile* fout = new TFile("tilt_correction_output.root", "RECREATE");
    g->Write("g_centroid_vs_thetaCM");
    fLin->Write("fLin_tilt");
    h2corr->Write("h2_Ex_corrected");
    fout->Close();

    std::cout << "Guardado tilt_correction_output.root con el grafico, el fit y el h2 corregido." << std::endl;
}