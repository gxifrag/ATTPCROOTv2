// clang-format off
R__LOAD_LIBRARY(/home/georgina/fair_install/ATTPCROOTv2/macro/a1975/23April_macros/SpectralModel_cxx.so)
// clang-format on

// root -l
// .L SpectralModel.cxx++
#include <Math/MinimizerOptions.h>
#include <TLine.h>
#include <TPaveText.h>
#include <TVectorD.h>

#include "SpectralModel.h" // <-- la clase modular: SUSTITUYE a la vieja "class SpectralModel{...}" del macro
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TSpectrum.h"
#include "TSystem.h"
#include "TTree.h"

#include <iostream>

// Si en algún momento quieres usar AddBWTabulated() (fallback con tus tablas
// interpoladas) en vez de AddBW() (analítica), necesitarás estos includes y
// los TGraph gL0/gL1/gL2. Si no los usas, puedes comentarlos sin problema.
/*
#include
"/home/georgina/fair_install/ATTPCROOTv2/macro/a1975/23April_macros/penetrabilities/penetrabilities_neutron_15C_L_long0.C"
#include
"/home/georgina/fair_install/ATTPCROOTv2/macro/a1975/23April_macros/penetrabilities/penetrabilities_neutron_15C_L_long1.C"
#include
"/home/georgina/fair_install/ATTPCROOTv2/macro/a1975/23April_macros/penetrabilities/penetrabilities_neutron_15C_L_long2.C"

TGraph *gL0 = new TGraph(num_points_L_0, energies_neutron_15C_L_0, T0_neutron_15C_values);
TGraph *gL1 = new TGraph(num_points_L_1, energies_neutron_15C_L_1, T1_neutron_15C_values);
TGraph *gL2 = new TGraph(num_points_L_2, energies_neutron_15C_L_2, T2_neutron_15C_values);
*/

const double Sn = 1.2181;

TGraph *histoToTgraph(TH1F *h)
{
   auto g = new TGraph();
   for (int i = 1; i <= h->GetNbinsX(); i++) {
      g->SetPoint(i - 1, h->GetBinCenter(i), h->GetBinContent(i));
   }
   g->SetName("graphPS_1n");
   return g;
}

void C16_pd_fitter()
{
   // Int_t nbins;
   /*TFile *fIn = new TFile("hexCorr2_pd_z_2-60cm_ke_5-15MeV_80bins.root", "READ");

   TH1F *hexCorr2 = (TH1F *)fIn->Get("hexCorr2");
   hexCorr2->SetDirectory(nullptr);*/

   TFile *fIn = new TFile("hexCorr1_pd_z_2-30cm_ke_5-15MeV_80bins.root", "READ");

    TH1F *hexCorr2 = (TH1F *)fIn->Get("hexCorr1");
   hexCorr2->SetDirectory(nullptr);

   double PlotMax =9., PlotMin = -2.;

   // AUTOMÁTICO: El fitter se adapta a lo que tenga el histograma real
   int NumberBins = hexCorr2->GetNbinsX();
   double Ebin_min = hexCorr2->GetXaxis()->GetXmin();
   double Ebin_max = hexCorr2->GetXaxis()->GetXmax();
   std::cout << "NumberBins = " << NumberBins << std::endl;
   std::cout << "Ebin_max = " << Ebin_max << std::endl;
   std::cout << "Ebin_min = " << Ebin_min << std::endl;

   TH1F *h_PS_1n = new TH1F("h_PS_1n", "h_PS_1n", NumberBins, Ebin_min, Ebin_max);

   //-------Phase Space--------
   double ThetaCM_min = 0;
   double ThetaCM_max = 180;

   TFile *filePS =
      new TFile("/home/georgina/fair_install/ATTPCROOTv2/macro/a1975/PhaseSpace/PhaseSpace_16C_pd_1n_cuts_z_2-60cm_theta_10-60deg_ke_5-15MeV_sigmaEx_0.203342MeV_80bins.root", "READ");
   TTree *treePS = (TTree *)filePS->Get("simulated_tree");
   if (!treePS) {
      std::cout << "Error!";
      return;
   }

   double Weight_sim, Ex_cal, ThetaCM_cal;
   treePS->SetBranchAddress("Weight_sim", &Weight_sim);
   treePS->SetBranchAddress("Ex_cal", &Ex_cal);
   treePS->SetBranchAddress("ThetaCM_cal", &ThetaCM_cal);

   for (int i = 0; i < treePS->GetEntries(); i++) {
      treePS->GetEntry(i);
      if (ThetaCM_cal > ThetaCM_min && ThetaCM_cal < ThetaCM_max) {
         h_PS_1n->Fill(Ex_cal, Weight_sim);
      }else {
         std::cout << "no cumple la condición del pS" << std::endl;
      }
   }

   h_PS_1n->Smooth();
   TH1F *hPS_prefit = (TH1F *)h_PS_1n->Clone("hPS_prefit");
   TH1F *hPS_final = (TH1F *)h_PS_1n->Clone("hPS_final");
   TGraph *graphPS = histoToTgraph(h_PS_1n);

   //----------------------------minimization
   ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
   TCanvas *c_prefits = new TCanvas("prefits", "", 1200, 800);
   c_prefits->cd();

   TSpectrum *sp = new TSpectrum(6);
   int nfound = sp->Search(hexCorr2, 2, "", 0.02);
   Double_t *xpeaks = sp->GetPositionX();

   std::vector<double> sorted_peaks(xpeaks, xpeaks + nfound);
   std::sort(sorted_peaks.begin(), sorted_peaks.end());

   std::cout << "Picos encontrados por TSpectrum: " << nfound << std::endl;
   for (int i = 0; i < nfound; i++)
      std::cout << "Pico " << i + 1 << ": " << sorted_peaks[i] << " MeV" << std::endl;
   if (nfound < 4)
      std::cout << "[WARNING] TSpectrum no encontro los 4 picos." << std::endl;

   double m1 = sorted_peaks[0];
   double m2 = sorted_peaks[1];
   /*double m3 = sorted_peaks[2];
   double m4 = sorted_peaks[3];*/

   // ---- Prefit de las 2 gaussianas de fondo: SOLO para extraer sigma_mean ----
   TF1 *f2g = new TF1("f2g", "gaus(0) + gaus(3)", -1.0, 1.6);
   f2g->SetParameter(1, m1);
   f2g->SetParameter(2, 0.2);
   f2g->SetParameter(0, hexCorr2->GetBinContent(hexCorr2->FindBin(m1)));
   f2g->SetParameter(4, m2);
   f2g->SetParameter(5, 0.2);
   f2g->SetParameter(3, hexCorr2->GetBinContent(hexCorr2->FindBin(m2)));

   hexCorr2->Fit(f2g, "RQ");
   f2g->SetLineColor(kViolet + 2);
   f2g->SetNpx(500);
   hexCorr2->Draw("hist");
   f2g->Draw("same");

   double Ex1 = f2g->GetParameter(1);
   double sigma1 = f2g->GetParameter(2);
   double Ex2 = f2g->GetParameter(4);
   double sigma2 = f2g->GetParameter(5);
   double sigma_mean = 0.5 * (sigma1 + sigma2);

   std::cout << "sigma_mean = " << sigma_mean << " MeV" << std::endl;

   double phaseSpace = 0.0017;
   if (hPS_prefit->GetNbinsX() > 0)
      hPS_prefit->Scale(phaseSpace);
   hPS_prefit->SetLineColor(kGray + 2);
   hPS_prefit->SetLineWidth(2);
   hPS_prefit->SetLineStyle(2);
   hPS_prefit->SetMarkerSize(0);
   hPS_prefit->Draw("same HIST");
   
   //======================================================================
   // ---- MODELO ESPECTRAL: aqui esta el cambio principal ----
   //======================================================================

   // TODO: pon aqui los valores reales (los mismos que usa tu compañero)
   const double mass_C14_u = 14.003242;                                // 14C masa en amu
   const double mass_n_u = 1.008665;                                   // neutrón en amu
   const double R_ch = 1.2 * std::pow(14, 1.0 / 3.0);                  // Interaction radius in fm
   double mu_n14C = (mass_C14_u * mass_n_u) / (mass_C14_u + mass_n_u); // masa reducida del canal n+14C, en MeV/c^2
                                                                       // radio de canal, en fm

   // TODO: anchuras de literatura (o tus límites conocidos), en MeV.
   // Estas son Gamma0 = Gamma EN la resonancia, directamente -- ya no hace
   // falta convertir a gamma^2 como antes.
   double Gamma1_lit = 0.04;
   double Gamma2_lit = 0.014;
   double Gamma3_lit = 1.7;
   double Gamma4_lit = 0.02;

   // SpectralModel model;
   SpectralModel *model = new SpectralModel(); // USAR PUNTERO PARA EVITAR EL SEG FAULT

   model->AddGaussian("g0", f2g->GetParameter(0), f2g->GetParameter(1), f2g->GetParameter(2));
   model->AddGaussian("g1", f2g->GetParameter(3), f2g->GetParameter(4), f2g->GetParameter(5));

   // L=1, L=1, L=2 -- ajusta si tus 3 resonancias tienen otros L
   model->AddBW("bw0", /*L=*/1, Sn, mu_n14C, R_ch, /*ampInit=*/120, /*ER*/ 3.5, sigma_mean, Gamma1_lit);
   model->AddBW("bw1", /*L=*/2, Sn, mu_n14C, R_ch, /*ampInit=*/50, /*ER*/ 4.0, sigma_mean, Gamma2_lit);
   //model->AddBW("bw2", /*L=*/2, Sn, mu_n14C, R_ch, /*ampInit=*/70, /*ER*/ 4.2, sigma_mean, Gamma3_lit);
   // Add the new 5th peak from TSpectrum
   model->AddBW("bw2", /*L=*/2, Sn, mu_n14C, R_ch, /*ampInit=*/50, /*ER*/ 4.4, sigma_mean, Gamma4_lit);
   //model->AddBW("bw3", /*L=*/1, Sn, mu_n14C, R_ch, 30,  /*ER*/ 5.2, sigma_mean, Gamma4_lit);

   model->SetPhaseSpace(graphPS, phaseSpace);
   model->Finalize();

   TCanvas *c_ExEner = new TCanvas("ExEner", "Corrected Excited Energy spectra", 1200, 800);
   c_ExEner->cd();

   TF1 *fModel = model->Build("fModel_conv",Ebin_min-2., 8.2); //changing from Ebin_min-2. to Ebin_min, the FFT explodes
   fModel->SetNpx(1000);
   fModel->SetLineColor(kBlack);
   fModel->SetLineWidth(4);

   // Sigmas instrumentales fijas (las 3 comparten la misma resolución)0
   fModel->FixParameter(model->Idx("bw0_Sigma"), sigma_mean); // kGreen
   fModel->FixParameter(model->Idx("bw1_Sigma"), sigma_mean); // kMagenta
   fModel->FixParameter(model->Idx("bw2_Sigma"), sigma_mean); // kRed
   //fModel->FixParameter(model->Idx("bw3_Sigma"), sigma_mean); // kViolet

   // Limites de energia/anchura -- ahora con nombres, no con p[7], p[9]...
   // fModel->SetParLimits(model->Idx("bw0_ER"), 3.25, 3.50);
   fModel->SetParLimits(model->Idx("bw0_Gamma0"), 0.001, 1.0);

   fModel->SetParLimits(model->Idx("bw1_ER"), 3.9, 4.2);
   // fModel->SetParLimits(model->Idx("bw1_Gamma0"), 0.005, 0.08);
  
   fModel->SetParLimits(model->Idx("bw2_ER"), 4.8, 5.8); // Estreta el rang segons la teva esperança física
   fModel->SetParLimits(model->Idx("bw2_Gamma0"), 0.05, 1.5); // Limita l'amplada 
   //fModel->SetParLimits(model->Idx("bw3_ER"), 4.8, 5.3); // Estreta el rang segons la teva esperança física
   //fModel->SetParLimits(model->Idx("bw3_Gamma0"), 0.05, 1.0); // Limita l'amplada 
   //fModel->FixParameter(model->Idx("ps0_Amp"), 0.003);

   //..........................................................

   hexCorr2->Fit(fModel);
   int npar = fModel->GetNpar();
   for (int i = 0; i < npar; ++i) {
      double v = fModel->GetParameter(i);
      if (!std::isfinite(v)) {
         std::cerr << "Param[" << i << "] (" << model->ParName(i) << ") invalid: " << v << "\n";
         return;
      }
   }

  // Frame vacío con el rango y divisiones que quieres
TH1F *frame = (TH1F*)gPad->DrawFrame(PlotMin, 0.5, PlotMax, hexCorr2->GetMaximum()*1.1);
frame->GetXaxis()->SetTitle("Excitation Energy (MeV)");
frame->GetYaxis()->SetTitle("Counts / 125 keV");
frame->GetXaxis()->SetNdivisions(-11, kFALSE);   // ahora SÍ coincide: 12 sobre un rango de 12
gPad->Modified();
gPad->Update();

// Ahora dibuja encima
hexCorr2->Draw("same E1");
fModel->Draw("same L");

   // ---- Componentes individuales para el dibujo: ya no hay TF1 manuales ----
   std::vector<int> gaussColors = {kOrange + 7, kBlue};
   int colIdx = 0;
   for (const auto &label : model->GaussLabels()) {
      TF1 *fg = model->GetComponentTF1(label, fModel, PlotMin, 8.2);
      fg->SetLineColor(gaussColors[colIdx % gaussColors.size()]);
      fg->SetLineStyle(2);
      fg->SetNpx(500);
      fg->Draw("same L");
      colIdx++;
   }

   std::vector<int> bwColors = {kGreen + 2, kMagenta, kRed + 2};
   colIdx = 0;
   for (const auto &label : model->BWLabels()) {
      TF1 *fb = model->GetComponentTF1(label, fModel, PlotMin, 8.2);
      fb->SetLineColor(bwColors[colIdx % bwColors.size()]);
      fb->SetLineStyle(2);
      fb->SetNpx(500);
      fb->Draw("same L");
      colIdx++;
   }

   if (model->HasPhaseSpace()) {
      TF1 *fps = model->GetComponentTF1("ps0", fModel, PlotMin, 8.2);
      fps->SetLineColor(kGray + 2);
      fps->SetLineStyle(2);
      fps->SetFillColor(kGray + 2);   // color de las rayas
      fps->SetFillStyle(3004);        // rayas diagonales (prueba 3005 para el otro sentido)
      fps->SetLineColor(kGray + 2);   // contorno del mismo color, o distinto si prefieres
      fps->SetLineWidth(2);
      fps->Draw("same HIST");         // importante: HIST para que se vea el relleno, no solo el contorno
      fps->SetNpx(3000);
      fps->Draw("same L");

   }
  
   double chi2 = fModel->GetChisquare();
   int ndf = fModel->GetNDF();
   double chi2Ndf = chi2 / ndf;

   TPaveText *pt = new TPaveText(0.6, 0.63, 0.88, 0.85, "NDC");
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetBorderSize(0);
   pt->SetTextAlign(11);
   pt->SetTextSize(0.032);
   pt->SetTextFont(42);
   pt->AddText("^{16}C(p,d)^{15}C  E_{beam}/A = 11.5 MeV");
   pt->AddText(Form("#sigma_{det} = %.0f keV (from g.s. + 1st)", sigma_mean * 1000.0));
   pt->AddText(Form("#chi^{2}/NDF = %.2f", chi2Ndf));
   pt->Draw("same");

   gPad->Update();
   double ymax = gPad->GetUymax();
   TLine *vline0 = new TLine(Sn, 0, Sn, ymax);
   vline0->SetLineColor(kRed);
   vline0->SetLineWidth(3);
   vline0->Draw("SAME");

   c_ExEner->Update();

   //---------------------------------------------------
   // --- Guardar parametros del fit global en ROOT ---
   // Guardamos automáticamente con el número de bins correcto en el nombre
   TFile *fFitParams = new TFile(Form("fit_params_global_%dbins.root", NumberBins), "RECREATE");
   fModel->Write("fModel_global");
   if (graphPS)
      graphPS->Write("graphPS");

   TVectorD params(npar), errors(npar);
   for (int i = 0; i < npar; ++i) {
      params[i] = fModel->GetParameter(i);
      errors[i] = fModel->GetParError(i);
   }
   params.Write("fit_parameters");
   errors.Write("fit_errors");

   TNamed chi2_str("chi2_ndf", Form("%.4f / %d = %.4f", chi2, ndf, chi2Ndf));
   chi2_str.Write();
   fFitParams->Close();
   std::cout << "Fit params guardados en fit_params_global.root" << std::endl;

   // ---- Comprobacion del threshold para cada resonancia ----
   for (const auto &label : model->BWLabels()) {
      double ER = fModel->GetParameter(model->Idx(label + "_ER"));
      double sigma = fModel->GetParameter(model->Idx(label + "_Sigma"));
      double lower = ER - 5 * sigma;
      std::cout << label << " ER=" << ER << "  rango conv. inferior = " << lower << "  Sn=" << Sn
                << (lower < Sn ? "  --> CRUZA EL THRESHOLD" : "") << std::endl;
   }
}
