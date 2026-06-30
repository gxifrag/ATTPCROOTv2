#include "TCanvas.h"
#include "TF1.h"
#include "TF1Convolution.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
// #include "TH2.h"
#include <Math/MinimizerOptions.h>
#include <TLine.h>
#include <TPaveText.h>
#include <TVectorD.h>

#include "TProfile.h"
#include "TROOT.h"
#include "TSpectrum.h"
#include "TSystem.h"
#include "TTree.h"

// #include <fstream>
#include <iostream>

#include "/home/georgina/fair_install/ATTPCROOTv2/macro/a1975/23April_macros/penetrabilities/penetrabilities_neutron_15C_L_long0.C"
#include "/home/georgina/fair_install/ATTPCROOTv2/macro/a1975/23April_macros/penetrabilities/penetrabilities_neutron_15C_L_long1.C"
#include "/home/georgina/fair_install/ATTPCROOTv2/macro/a1975/23April_macros/penetrabilities/penetrabilities_neutron_15C_L_long2.C" // o el L que corresponda

// Crear los TGraph

TGraph *gL0 = new TGraph(num_points_L_0, energies_neutron_15C_L_0, T0_neutron_15C_values);
TGraph *gL1 = new TGraph(num_points_L_1, energies_neutron_15C_L_1, T1_neutron_15C_values);
TGraph *gL2 = new TGraph(num_points_L_2, energies_neutron_15C_L_2, T2_neutron_15C_values);

double Ebin_max = 9.;
double Ebin_min = -1.;
int NumberBins = 80; // best chi2/ndf for 80 bins, 0-8 MeV
int NumberBinsAux = 200;
const double Sn = 1.2181;

// double BreitWignerPenetrability(double E, double ER, double gamma2, TSpline3 *splinePen)
double BreitWignerPenetrability(double E, double ER, double gamma2, TGraph *gPen)
{

   double En = E - Sn;
   if (En < 0)
      return 0;

   double En_max = 7.78; // o el E_end que usaste al generar la tabla
   if (En > En_max)
      En = En_max; // clamp en vez de extrapolar sin control

   double P = gPen->Eval(En); // ← interpolación lineal
   // double P = gPen->Eval(En, nullptr, "S");

   double GammaE = 2.0 * P * gamma2;
   if (GammaE <= 0)
      return 0.0;
   double denom = (E - ER) * (E - ER) + 0.25 * GammaE * GammaE;

   if (denom < 1e-12)
      return 0.0;

   return GammaE / denom;
}

double BreitWignerPenetrabilityForTF1(double *x, double *p)
{
   double xx = x[0];
   double En = xx - Sn;
   double er = p[0];
   double gamma2 = p[1];
   if (En < 0)
      return 0;

   double En_max = 7.78; // o el E_end que usaste al generar la tabla
   if (En > En_max)
      En = En_max; // clamp en vez de extrapolar sin control

   double P = gL0->Eval(En); // ← interpolación lineal
   // double P = gPen->Eval(En, nullptr, "S");

   double GammaE = 2.0 * P * gamma2;
   if (GammaE <= 0)
      return 0.0;
   double denom = (xx - er) * (xx - er) + 0.25 * GammaE * GammaE;

   if (denom < 1e-12)
      return 0.0;

   return GammaE / denom;
}

double ConvolutedBW(double x, double ER, double sigma, double gamma2, TGraph *gPen)
{
   const int N = 1000;
   const double range = 5 * sigma;
   const double step = 2 * range / N;
   double sum = 0;

   for (int i = 0; i < N; i++) {
      double t = x - range + i * step;
      double bw = BreitWignerPenetrability(t, ER, gamma2, gPen);
      double ga = TMath::Gaus(x - t, 0, sigma, true);
      sum += bw * ga;
   }

   return sum * step;
}

class SpectralModel {
public:
   TGraph *graphPS;
   /* TSpline3 *splL0;
    TSpline3 *splL1;
    TSpline3 *splL2;*/

   TGraph *gL0;
   TGraph *gL1;
   TGraph *gL2;

   // Create pointer for convolution functions
   TF1 *fBW{};
   TF1 *fGauss{};
   TF1Convolution *fConv{};

   // SpectralModel(TGraph *ps, TSpline3 *s0, TSpline3 *s1, TSpline3 *s2) : graphPS(ps), splL0(s0), splL1(s1), splL2(s2)
   SpectralModel(TGraph *ps, TGraph *s0, TGraph *s1, TGraph *s2) : graphPS(ps), gL0(s0), gL1(s1), gL2(s2)
   {
      fBW = new TF1{"fBW", BreitWignerPenetrabilityForTF1, -10, 15, 2};
      fGauss = new TF1{"fGauss", "TMath::Gaus(x, 0, [1])", -10, 15};
      fConv = new TF1Convolution(fBW, fGauss);
   }

   double operator()(double *x, double *p)
   {
      double val = 0;
      val += p[0] * TMath::Gaus(x[0], p[1], p[2], false);
      val += p[3] * TMath::Gaus(x[0], p[4], p[5], false);

      // Set paramters for convolution of 1st unbound
      double er0 = p[7];
      double gamma0 = p[9];
      double sigma0 = p[8];
      // fConv->SetParameters(er0, gamma0, sigma0);

      // *** Componente 1: BW+penetrabilidad, convolucionada ***
      double amp1 = p[6];
      double ER1 = p[7];
      double sigma1 = p[8];
      double gamma1 = p[9];
      // val += amp1 * (*fConv)(&x[0], nullptr);
      //  val += amp1 * ConvolutedBW(x[0], ER1, sigma1, gamma1, gL1); // <- ajusta L
      val += amp1 * BreitWignerPenetrability(x[0], ER1, gamma1, gL1);

      // *** Componente 2: BW+penetrabilidad, convolucionada ***
      double amp2 = p[10];
      double ER2 = p[11];
      double sigma2 = p[12];
      double gamma2 = p[13];
      val += amp2 * BreitWignerPenetrability(x[0], ER2, gamma2, gL1);

      // *** Componente 3: BW+penetrabilidad, convolucionada ***
      double amp3 = p[14];
      double ER3 = p[15];
      double sigma3 = p[16];
      double gamma3 = p[17];
      val += amp3 * BreitWignerPenetrability(x[0], ER3, gamma3, gL2);

      val += p[18] * graphPS->Eval(x[0]);

      return val;
   }
};

TGraph *histoToTgraph(TH1F *h)
{

   auto g = new TGraph();
   for (int i = 1; i <= h->GetNbinsX(); i++) {
      g->SetPoint(i - 1, h->GetBinCenter(i), h->GetBinContent(i));
   }

   g->SetName("graphPS_1n");
   return g;
}

void C16_pd_fitter_old()
{
   Int_t nbins;
   TFile *fIn = new TFile("hexCorr2_pd_z_2-60cm_ke_5-15MeV.root", "READ");
   TH1F *hexCorr2 = (TH1F *)fIn->Get("hexCorr2");
   hexCorr2->SetDirectory(nullptr); // opcional, para que sobreviva si cierras fIn*/

   TH1F *h_PS_1n = new TH1F("h_PS_1n", "h_PS_1n", NumberBins, Ebin_min, Ebin_max);

   //-------Phase Space--------

   nbins = hexCorr2->GetNbinsX();
   int binmax = hexCorr2->GetMaximumBin();
   double ThetaCM_min = 0;
   double ThetaCM_max = 180;

   TFile *filePS =
      new TFile("/home/georgina/fair_install/ATTPCROOTv2/macro/a1975/PhaseSpace/PhaseSpace_16C_pd_1n.root", "READ");
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

   TSpectrum *sp = new TSpectrum(6);               // 5 maxima search
   int nfound = sp->Search(hexCorr2, 2, "", 0.02); // 2 = sigma of smoothing, last = threshold
   Double_t *xpeaks = sp->GetPositionX();

   // copy into a vector<double>
   std::vector<double> sorted_peaks(xpeaks, xpeaks + nfound);
   // sort them
   std::sort(sorted_peaks.begin(), sorted_peaks.end());

   double m1 = sorted_peaks[0];
   double m2 = sorted_peaks[1];
   double m3 = sorted_peaks[2];
   double m4 = sorted_peaks[3];

   std::cout << "Picos encontrados por TSpectrum: " << nfound << std::endl;
   std::cout << "Picos ordenados:" << std::endl;
   for (int i = 0; i < nfound; i++) {
      std::cout << "Pico " << i + 1 << ": " << sorted_peaks[i] << " MeV" << std::endl;
   }

   if (nfound < 4) {
      std::cout << "[WARNING] TSpectrum no encontró los 4 picos. Usando valores por defecto para los faltantes."
                << std::endl;
   }

   // Define a two-gaussian TF1 (ROOT built-in gaus uses amplitude = height)
   TF1 *f2g = new TF1("f2g", "gaus(0) + gaus(3)", -1.0, 1.6);

   // Peak 1
   f2g->SetParameter(1, m1);                                             // mean of gaus(0)
   f2g->SetParameter(2, 0.2);                                            // sigma guess
   f2g->SetParameter(0, hexCorr2->GetBinContent(hexCorr2->FindBin(m1))); // amplitude guess

   // Peak 2
   f2g->SetParameter(4, m2);                                             // mean of gaus(3)
   f2g->SetParameter(5, 0.2);                                            // sigma guess
   f2g->SetParameter(3, hexCorr2->GetBinContent(hexCorr2->FindBin(m2))); // amplitude guess

   hexCorr2->Fit(f2g, "RQ"); // R = use the range you specified, 0 no plot, M minuit
   f2g->SetLineColor(kViolet + 2);
   f2g->SetNpx(500); // aumenta el número de puntos para un ajuste más suave

   hexCorr2->Draw("hist");
   f2g->Draw("same"); // <--- REQUIRED so the fit curve is drawn

   // Extraer valores
   double Ex1 = f2g->GetParameter(1);
   double sigma1 = f2g->GetParameter(2);

   double Ex2 = f2g->GetParameter(4);
   double sigma2 = f2g->GetParameter(5);

   // Crear arrays AHORA (cuando ya existen los valores)
   double Ex_values[2] = {Ex1, Ex2};
   double sigma_values[2] = {sigma1, sigma2};

   /*TF1 *bwprefit1 = new TF1("bw1", "[0]*TMath::BreitWigner(x,[1],[2])", 2.5, 3.8);
   bwprefit1->SetParameters(150, m3, 0.01);
   bwprefit1->SetLineColor(kRed);
   bwprefit1->SetNpx(500); // aumenta el número de puntos para un ajuste más suave
   hexCorr2->Fit(bwprefit1, "RQ");
   bwprefit1->Draw("same");

   TF1 *bwprefit2 = new TF1("bw2", "[0]*TMath::BreitWigner(x,[1],[2])", 3.80, 4.2);
   bwprefit2->SetParameters(50, 4.2, 1.0);
   // bwprefit2->SetParLimits(2, 0.1, 1.0); // gamma acotado
   bwprefit2->SetLineColor(kBlue);
   bwprefit2->SetNpx(500); // aumenta el número de puntos para un ajuste más suave
   hexCorr2->Fit(bwprefit2, "RQ+");
   bwprefit2->Draw("same");

   TF1 *bwprefit3 = new TF1("bw3", "[0]*TMath::BreitWigner(x,[1],[2])", 4.5, 6.0);
   bwprefit3->SetParameters(50, m4, 1.5);
   // bwprefit3->SetParLimits(2, 0., 1.8); // gamma acotado
   bwprefit3->SetLineColor(kOrange);
   bwprefit3->SetNpx(500); // aumenta el número de puntos para un ajuste más suave
   hexCorr2->Fit(bwprefit3, "RQ+");
   bwprefit3->Draw("same");*/

   double phaseSpace = 0.0003; // valor inicial para el ajuste, se puede ajustar según la escala de los datos

   if (hPS_prefit->GetNbinsX() > 0) {
      hPS_prefit->Scale(phaseSpace); // Normaliza a la integral deseada
   }

   hPS_prefit->SetLineColor(kGray + 2);
   hPS_prefit->SetLineWidth(2);
   hPS_prefit->SetLineStyle(2);   // ← 2 = línea discontinua
   hPS_prefit->SetMarkerSize(0);  // ← elimina los puntos
   hPS_prefit->Draw("same HIST"); // ← HIST fuerza línea, sin marcadores

   /*double A1 = bwprefit1->GetParameter(0);
   double mean1 = bwprefit1->GetParameter(1);
   double gamma1 = bwprefit1->GetParameter(2);

   std::cout << "gamma1 =  " << gamma1 << std::endl;

   double A2 = bwprefit2->GetParameter(0);
   double mean2 = bwprefit2->GetParameter(1);
   double gamma2 = bwprefit2->GetParameter(2);

   double A3 = bwprefit3->GetParameter(0);
   double mean3 = bwprefit3->GetParameter(1);
   double gamma3 = bwprefit3->GetParameter(2);*/

   // Crear gráfico
   TGraph *gSigma = new TGraph(2, Ex_values, sigma_values);
   gSigma->SetTitle("Sigma vs Ex;E_{x} (MeV);#sigma (MeV)");
   gSigma->GetXaxis()->SetLimits(-1, 8);     // fija el rango en X
   gSigma->GetYaxis()->SetRangeUser(0, 0.5); // si quieres también rango en Y (opcional)
   gSigma->SetMarkerStyle(20);
   gSigma->SetMarkerColor(kRed + 1);
   gSigma->SetLineColor(kRed + 1);

   double Ex_mean = 0.5 * (Ex1 + Ex2);
   double sigma_mean = 0.5 * (sigma1 + sigma2);

   std::cout << "Ex_mean = " << Ex_mean << " MeV" << std::endl;
   std::cout << "sigma_mean = " << sigma_mean << " MeV" << std::endl;

   /*TGraph *gMean = new TGraph(1);
   gMean->SetPoint(0, Ex_mean, sigma_mean);
   gMean->SetMarkerStyle(29);
   gMean->SetMarkerSize(2.0);
   gMean->SetMarkerColor(kBlue + 2);*/

   TCanvas *c_ExEner = new TCanvas("ExEner", "Corrected Excited Energy spectra", 1200, 800);
   c_ExEner->cd();
   SpectralModel *model = new SpectralModel(graphPS, gL0, gL1, gL2);
   TF1 *fModel = new TF1("fModel_conv", model, -1.0, 8.5, 19, "SpectralModel");

   fModel->SetNpx(500);
   fModel->SetLineColor(kBlack);
   fModel->SetLineWidth(4);

   std::vector<double> globalParamsIni = {
      f2g->GetParameter(0),
      f2g->GetParameter(1),
      f2g->GetParameter(2), // p0-2  Gaus1
      f2g->GetParameter(3),
      f2g->GetParameter(4),
      f2g->GetParameter(5), // p3-5  Gaus2

      120,
      3.5,
      sigma_mean,
      1., // p6-9   Comp1 (conv)

      50,
      3.9,
      sigma_mean,
      1., // p10-13 Comp2 (conv)

      70,
      5.2,
      sigma_mean,
      1., // p14-17 Comp3 (conv)

      phaseSpace // p18
   };

   // Assign parameters
   for (size_t i = 0; i < globalParamsIni.size(); ++i) {
      fModel->SetParameter(i, globalParamsIni[i]);
   }

   // Sigmas instrumentales fijas (las 3 componentes comparten la misma resolución)
   /* fModel->FixParameter(8, sigma_mean);  // sigma Comp1
    fModel->FixParameter(12, sigma_mean); // sigma Comp2
    fModel->FixParameter(16, sigma_mean); // sigma Comp3
    */

   // Límites de energía/anchura
   fModel->SetParLimits(7, 2.8, 3.7);   // ER Comp1
   fModel->SetParLimits(9, 0.001, 1.0); // antes 0.001 — ahora evita el colapso a anchura ~0
   // fModel->SetParLimits(10, 15., 1e6);
   fModel->FixParameter(10, 0.0);
   fModel->FixParameter(11, 0.0);
   fModel->FixParameter(12, 0.0);
   fModel->FixParameter(13, 0.0);

   // fModel->SetParLimits(11, 3.9, 4.2); // ER Comp2
   // fModel->SetParLimits(13, 0.2, 3.0); // gamma Comp2 (opcional, antes estaba comentado)

   fModel->SetParLimits(15, 4.5, 6.0);   // ER Comp3
   fModel->SetParLimits(17, 0.001, 3.0); // gamma Comp3

   //..........................................................

   hexCorr2->Fit(fModel, "R");

   int npar = fModel->GetNpar();
   std::vector<double> globalParamsFinals(npar);
   for (int i = 0; i < npar; ++i) {
      globalParamsFinals[i] = fModel->GetParameter(i);
   }

   for (int i = 0; i < npar; ++i) {
      double v = fModel->GetParameter(i);
      if (!std::isfinite(v)) {
         std::cerr << "Param[" << i << "] invalid: " << v << "\n";
         return;
      }
   }

   hexCorr2->GetXaxis()->SetTitle("Excitation Energy (MeV)");
   hexCorr2->GetYaxis()->SetTitle("Counts");
   hexCorr2->Draw("E1"); // E1
   fModel->Draw("same L");
   // 1. Forzar la actualización del canvas actual
   gPad->Update();

   // 2. Pedirle al sistema que procese los eventos gráficos pendientes
   gSystem->ProcessEvents();

   // Gaussian 1
   TF1 *gaus1 = new TF1("gaus1", "gaus(0)", Ebin_min, Ebin_max);
   gaus1->SetParameters(globalParamsFinals[0], globalParamsFinals[1], globalParamsFinals[2]);
   gaus1->SetLineColor(kOrange + 7);
   gaus1->SetNpx(500);     // aumenta el número de puntos para un ajuste más suave
   gaus1->SetLineStyle(2); // línea discontinua
   gaus1->Draw("same L");

   // Gaussian 2
   TF1 *gaus2 = new TF1("gaus2", "gaus(0)", Ebin_min, Ebin_max);
   gaus2->SetParameters(globalParamsFinals[3], globalParamsFinals[4], globalParamsFinals[5]);
   gaus2->SetLineColor(kBlue);
   gaus2->SetNpx(500);     // aumenta el número de puntos para un ajuste más suave
   gaus2->SetLineStyle(2); // línea discontinua
   gaus2->Draw("same L");

   // Componente 1
   TF1 *comp1 = new TF1(
      "comp1",
      [=](double *x, double *p) {
         return p[0] * ConvolutedBW(x[0], p[1], p[2], p[3], model->gL1);
      }, // mismo L que arriba
      Ebin_min, Ebin_max, 4);
   comp1->SetParameters(globalParamsFinals[6], globalParamsFinals[7], globalParamsFinals[8], globalParamsFinals[9]);
   comp1->SetLineColor(kGreen + 2);
   comp1->SetLineStyle(2);
   comp1->SetNpx(500);
   comp1->Draw("same L");

   // Componente 2
   TF1 *comp2 = new TF1(
      "comp2",
      [=](double *x, double *p) {
         return p[0] * ConvolutedBW(x[0], p[1], p[2], p[3], model->gL1);
      }, // mismo L que arriba
      Ebin_min, Ebin_max, 4);
   comp2->SetParameters(globalParamsFinals[10], globalParamsFinals[11], globalParamsFinals[12], globalParamsFinals[13]);
   comp2->SetLineColor(kMagenta);
   comp2->SetLineStyle(2);
   comp2->SetNpx(500);
   comp2->Draw("same L");

   // Componente 3
   TF1 *comp3 = new TF1(
      "comp3",
      [=](double *x, double *p) {
         return p[0] * ConvolutedBW(x[0], p[1], p[2], p[3], model->gL2);
      }, // mismo L que arriba
      Ebin_min, Ebin_max, 4);
   comp3->SetParameters(globalParamsFinals[14], globalParamsFinals[15], globalParamsFinals[16], globalParamsFinals[17]);
   comp3->SetLineColor(kCyan + 2);
   comp3->SetLineStyle(2);
   comp3->SetNpx(500);
   comp3->Draw("same L");
   // Phase Space
   if (hPS_final->GetNbinsX() > 0) {
      hPS_final->Scale(globalParamsFinals[18]); // globalParamsFinals[14] es el factor de escala ajustado para el fondo
                                                // de fase espacio
   }
   hPS_final->SetLineColor(kGray + 2);
   hPS_final->SetLineWidth(2);
   hPS_final->SetMarkerStyle(24); // open circle
   hPS_final->SetMarkerSize(1);
   hPS_final->SetMarkerColor(kGray + 2);
   hPS_final->Draw("same P"); // P = solo marcadores, sin línea

   // Suppose you fitted with 'fitFcn' (could be gaus1, bw1, etc.)
   double chi2 = fModel->GetChisquare();
   int ndf = fModel->GetNDF();
   double chi2Ndf = chi2 / ndf;

   TPaveText *pt = new TPaveText(0.6, 0.63, 0.88, 0.85, "NDC");
   pt->SetFillColor(0);  // fondo blanco
   pt->SetFillStyle(0);  // sin relleno (transparente)
   pt->SetBorderSize(0); // sin frame
   pt->SetTextAlign(11); // left-center
   pt->SetTextSize(0.032);
   pt->SetTextFont(42); // misma fuente que ROOT por defecto
   pt->AddText("^{16}C(p,d)^{15}C  E_{beam}/A = 11.5 MeV");
   pt->AddText(Form("#sigma_{det} = %.0f keV (from g.s. + 1st)", sigma_mean * 1000.0));
   pt->AddText(Form("#chi^{2}/NDF = %.2f", chi2Ndf));
   pt->Draw("same");

   gPad->Update();                 // Asegura que el pad esté actualizado para obtener los límites correctos
   double ymax = gPad->GetUymax(); // máximo del eje Y en coordenadas del pa

   TLine *vline0 = new TLine(1.218, 0, 1.218, ymax); // línea vertical en x=1.218 MeV;
   vline0->SetLineColor(kRed);                       // opcional
   vline0->SetLineWidth(3);                          // opcional
   vline0->Draw("SAME");

   std::vector<double> Ex_theory = {0.0, 0.740, 3.103, 4.202, 4.780};
   std::vector<int> line_colors = {kOrange + 7, kBlue, kGreen + 2, kMagenta, kCyan + 1};

   c_ExEner->Update();

   /* for (int i = 0; i < Ex_theory.size(); ++i) {
      TLine *vline = new TLine(Ex_theory[i], 0, Ex_theory[i], ymax);
      vline->SetLineColor(line_colors[i]);
      // vline->SetLineStyle(2); // discontinua
      vline->SetLineWidth(2);
      vline->Draw("SAME");
   }*/

   //---------------------------------------------------
   // --- Guardar parámetros del fit global en ROOT ---
   /*TFile *fFitParams = new TFile("fit_params_global.root", "RECREATE");

   // Guardamos la TF1 completa (contiene parámetros + errores + chi2)
   fModel->Write("fModel_global");

   // Guardamos también graphPS si existe
   if (graphPS) {
      graphPS->Write("graphPS");
   } else {
      std::cout << "AVISO: graphPS es NULL, no se guarda." << std::endl;
   }

   // TVectorD para acceso rápido desde otros macros
   int npar_save = fModel->GetNpar(); // 13
   TVectorD params(npar_save), errors(npar_save);
   for (int i = 0; i < npar_save; ++i) {
      params[i] = fModel->GetParameter(i);
      errors[i] = fModel->GetParError(i);
   }
   params.Write("fit_parameters");
   errors.Write("fit_errors");

   // Guardamos también chi2 y NDF como TNamed para referencia
   TNamed chi2_str("chi2_ndf", Form("%.4f / %d = %.4f", fModel->GetChisquare(), fModel->GetNDF(),
                                    fModel->GetChisquare() / fModel->GetNDF()));
   chi2_str.Write();

   fFitParams->Close();
   std::cout << "Fit params guardados en fit_params_global.root" << std::endl;

   double Sn_val = 1.2181;
   for (auto [amp_idx, ER_idx, sigma_idx] :
        std::vector<std::tuple<int, int, int>>{{6, 7, 8}, {10, 11, 12}, {14, 15, 16}}) {
      double ER = fModel->GetParameter(ER_idx);
      double sigma = fModel->GetParameter(sigma_idx);
      double lower = ER - 5 * sigma;
      std::cout << "Comp ER=" << ER << "  rango conv. inferior = " << lower << "  Sn=" << Sn_val
                << (lower < Sn_val ? "  --> CRUZA EL THRESHOLD" : "") << std::endl;
   }*/
}