#include <Math/ProbFuncMathCore.h>

#include <cmath>
#include <iostream>
#include <vector>

// ============================================================
// angular_distribution_per_state.C
//
// Para cada bin de θ_CM:
//   1. Proyecta el espectro Ex desde el TChain de datos reales
//   2. Fita con SpectralModel (centroides/anchuras fijos del fit global,
//      solo amplitudes libres)
//   3. Extrae N_i(θ_CM) = integral de cada componente
//   4. Divide por eficiencia hEff_tCM
//   5. Produce dN/dθ_CM por estado (listo para convertir a dσ/dΩ)
// ============================================================

double Ebin_max = 8.;
double Ebin_min = -1.;
int NumberBins = 100;

// Rangos angulares definidos por ti
std::vector<std::pair<double, double>> thetaRanges = {{20, 30}, {30, 40}, {40, 50}, {50, 60}};
Double_t omega(Double_t x, Double_t y, Double_t z)
{
   return sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}

std::tuple<double, double>
kine_2b(Double_t m1, Double_t m2, Double_t m3, Double_t m4, Double_t K_proj, Double_t thetalab, Double_t K_eject)
{
   double Et1 = K_proj + m1;
   double Et2 = m2;
   double Et3 = K_eject + m3;
   double Et4 = Et1 + Et2 - Et3;
   double m4_ex, Ex, theta_cm;
   double s, t, u;
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
   theta_cm *= TMath::RadToDeg();
   return std::make_tuple(Ex, theta_cm);
}

// Modelo espectral (misma clase que en v16)
// Necesita el graphPS del phase space — se carga desde fit_params_global.root
class SpectralModel {
public:
   TGraph *graphPS;
   SpectralModel(TGraph *g) : graphPS(g) {}
   double operator()(double *x, double *p)
   {
      double val = 0;
      val += p[0] * TMath::Gaus(x[0], p[1], p[2], false);
      val += p[3] * TMath::Gaus(x[0], p[4], p[5], false);
      val += p[6] * TMath::BreitWigner(x[0], p[7], p[8]);
      val += p[9] * TMath::BreitWigner(x[0], p[10], p[11]);
      double ps_val = graphPS->Eval(x[0]);
      val += p[12] * ps_val;
      return val;
   }
};

class SpectralModelExtended {
public:
   TGraph *graphPS;
   SpectralModelExtended(TGraph *g) : graphPS(g) {}
   double operator()(double *x, double *p)
   {
      double val = 0;
      val += p[0] * TMath::Gaus(x[0], p[1], p[2], false);
      val += p[3] * TMath::Gaus(x[0], p[4], p[5], false);
      val += p[6] * TMath::BreitWigner(x[0], p[7], p[8]);
      val += p[9] * TMath::BreitWigner(x[0], p[10], p[11]);
      val += p[12] * TMath::BreitWigner(x[0], p[13], p[14]);
      double ps_val = graphPS->Eval(x[0]);
      val += p[15] * ps_val;
      return val;
   }
};

void angular_distribution_per_state()
{
   // gROOT->ProcessLine(".X /home/georgina/fair_install/ATTPCROOTv2/macro/a1975/myStyle.C");
   gStyle->SetTitleAlign(23);
   gStyle->SetTitleX(0.5);
   gStyle->SetOptFit(1111); // muestra parámetros, χ², ndf, prob, etc.
   gStyle->SetStatX(0.93);  // más a la izquierda
   gStyle->SetStatY(0.98);  // más arriba

   TFile *fParams = new TFile("fit_params_global.root", "READ");
   TVectorD *params = (TVectorD *)fParams->Get("fit_parameters");
   if (!params) {
      cerr << "ERROR: no se encuentran fit_parameters\n";
      return;
   }

   // Ahora sí puedes usar params[]
   double mean_GS = (*params)[1];
   double sigma_GS = (*params)[2];

   double mean_1st = (*params)[4];
   double sigma_1st = (*params)[5];

   double mean_2nd = (*params)[7];
   double gamma_2nd = (*params)[8];

   double mean_3rd = (*params)[10];
   double gamma_3rd = (*params)[11];

   cout << "Parámetros de forma cargados:" << endl;
   cout << Form("  GS:  mean=%.3f  sigma=%.3f", mean_GS, sigma_GS) << endl;
   cout << Form("  1st: mean=%.3f  sigma=%.3f", mean_1st, sigma_1st) << endl;
   cout << Form("  2nd: mean=%.3f  gamma=%.3f", mean_2nd, gamma_2nd) << endl;
   cout << Form("  3rd: mean=%.3f  gamma=%.3f", mean_3rd, gamma_3rd) << endl;

   TGraph *graphPS = (TGraph *)fParams->Get("graphPS");
   if (!graphPS) {
      cerr << "AVISO: graphPS no encontrado en fit_params_global.root\n"
           << "  → Añade graphPS->Write(\"graphPS\") al bloque de guardado en v16\n"
           << "  → O cámbia la ruta aquí:\n";
      // Alternativa: cargarlo de otro archivo
      // TFile *fPS = new TFile("ruta/a/phaseSpace.root","READ");
      // graphPS = (TGraph*)fPS->Get("graphPS");
      return;
   }

   // ============================================================
   // 1. Reconstruir 2D (Ex vs θ_CM) desde TChain de datos reales
   // ============================================================
   Double_t m_p = 938.272076;
   Double_t m_d = 1875.612931;
   Double_t m_C15 = 13979.218707;
   Double_t m_C16 = 14914.533798;
   Double_t Ebeam_buff = 11.5 * 16;
   int Z_ej = 1;

   double densityH2 = 3.553e-5;
   AtTools::AtELossCATIMA elossH2(densityH2);
   elossH2.SetMaterial(catima::Material(1, 1));
   elossH2.SetProjectile(16, 6, 16.0147);

   // Binning θ_CM — mismo que en efficiency macro
   const int nBinsCM = 90;
   const double tCM_min = 0., tCM_max = 180.;

   // 2D: eje X = Ex, eje Y = θ_CM
   TH2F *h2D = new TH2F("h2D", "Ex vs #theta_{CM};Ex (MeV);#theta_{CM} (#circ)", NumberBins, Ebin_min, Ebin_max,
                        nBinsCM, tCM_min, 80);

   TChain *chain = new TChain("parquettree");

   std::set<int> excluded = {111, 121, 148, 149};

   for (int i = 104; i <= 189; i++) {
      if (excluded.count(i))
         continue;

      char name[64];
      std::snprintf(name, sizeof(name), "run_%04d_2H.root", i);

      std::string fullPath =
         "/home/georgina/C16_analysis/C16_H2/output_a1975_22April_tb510_MMG20/InterpSolver/InterpSolver_pd_root/" +
         std::string(name);

      chain->Add(fullPath.c_str());
   }

   cout << "Total entries (data): " << chain->GetEntries() << endl;

   Double_t theta{}, Brho{}, zPos{}, vx_pos{}, vy_pos{};
   chain->SetBranchAddress("polar", &theta);
   chain->SetBranchAddress("brho", &Brho);
   chain->SetBranchAddress("vertex_z", &zPos);
   chain->SetBranchAddress("vertex_x", &vx_pos);
   chain->SetBranchAddress("vertex_y", &vy_pos);

   Double_t m_ej = m_d;
   for (Long64_t i = 0; i < chain->GetEntries(); i++) {
      chain->GetEntry(i);
      Double_t p_ej = Brho * Z_ej * 2.99792458 / 10 * 1000;
      Double_t E_ej = TMath::Sqrt(p_ej * p_ej + m_ej * m_ej) - m_ej;
      double dist3D = TMath::Sqrt(vx_pos * vx_pos + vy_pos * vy_pos + zPos * zPos) * 100.0;
      Double_t Ebeam_at_z = elossH2.GetEnergy(Ebeam_buff, dist3D * 10.0);
      double kethe = 13.;
      double theta_lab_corr = theta - (2.0 * TMath::Pi() / 4000) * (E_ej - kethe);
      auto [ex, theta_cm] = kine_2b(m_C16, m_p, m_d, m_C15, Ebeam_at_z, theta_lab_corr, E_ej);

      if (zPos * 100 > 2.0 && zPos * 100 < 60.0 && E_ej > 5.0 && E_ej < 15.0)
         h2D->Fill(ex, theta_cm);
   }
   cout << "2D llenado: " << h2D->GetEntries() << " eventos" << endl;

   // ============================================================
   // 2. Fit por slice de θ_CM
   // ============================================================

   // Histogramas de distribución angular por estado (sin corregir aún)
   TH1F *hAng_GS = new TH1F("hAng_GS", "GS;#theta_{CM} (#circ);dN/d#theta_{CM}", nBinsCM, tCM_min, tCM_max);
   TH1F *hAng_1st = new TH1F("hAng_1st", "1st;#theta_{CM} (#circ);dN/d#theta_{CM}", nBinsCM, tCM_min, tCM_max);
   TH1F *hAng_2nd = new TH1F("hAng_2nd", "2nd;#theta_{CM} (#circ);dN/d#theta_{CM}", nBinsCM, tCM_min, tCM_max);
   TH1F *hAng_3rd = new TH1F("hAng_3rd", "3rd;#theta_{CM} (#circ);dN/d#theta_{CM}", nBinsCM, tCM_min, tCM_max);

   // Los mismos, corregidos por eficiencia
   TH1F *hAng_GS_corr =
      new TH1F("hAng_GS_corr", "GS corr;#theta_{CM} (#circ);dN/d#theta_{CM}", nBinsCM, tCM_min, tCM_max);
   TH1F *hAng_1st_corr =
      new TH1F("hAng_1st_corr", "1st corr;#theta_{CM} (#circ);dN/d#theta_{CM}", nBinsCM, tCM_min, tCM_max);
   TH1F *hAng_2nd_corr =
      new TH1F("hAng_2nd_corr", "2nd corr;#theta_{CM} (#circ);dN/d#theta_{CM}", nBinsCM, tCM_min, tCM_max);
   TH1F *hAng_3rd_corr =
      new TH1F("hAng_3rd_corr", "3rd corr;#theta_{CM} (#circ);dN/d#theta_{CM}", nBinsCM, tCM_min, tCM_max);

   TCanvas *hexthetaCanvas = new TCanvas("hexvsthetaCanvas", "h2D", 800, 600);
   hexthetaCanvas->Divide(3, 2); // dos pads: arriba el 2D, abajo la proyección

   // Pad 1: el 2D completo
   hexthetaCanvas->cd();
   h2D->GetXaxis()->SetTitle("Excitation Energy (MeV)");
   h2D->GetYaxis()->SetTitle("#theta_{CM} (#circ)");

   h2D->Draw("colz");

   hexthetaCanvas->Update();
   //------------------------------------------------------

   // Modelo global
   SpectralModel *fmodel = new SpectralModel(graphPS);
   TF1 *fModel = new TF1("fModel", fmodel, -1.0, 8.0, 13, // number fitted parameters
                         "SpectralModel");

   // Inicializar con parámetros del fit global
   for (int i = 0; i < 13; i++) {
      fModel->SetParameter(i, (*params)[i]);
   }

   // Canvas
   TCanvas *cFits = new TCanvas("cFits", "Fits por slice", 1400, 1000);
   cFits->Divide(2, 2);

   int ipad = 1;

   for (auto &range : thetaRanges) {

      double tmin = range.first;
      double tmax = range.second;

      int bin_min = h2D->GetYaxis()->FindBin(tmin);
      int bin_max = h2D->GetYaxis()->FindBin(tmax);

      // Slice
      TH1D *hSlice = h2D->ProjectionX(Form("hEx_%d_%d", (int)tmin, (int)tmax), bin_min, bin_max);

      cFits->cd(ipad);

      // Dibujar histograma
      hSlice->SetTitle(Form("%.0f #circ < #theta_{CM} < %.0f #circ", tmin, tmax));
      hSlice->Draw("hist");

      // Crear el fit del slice
      TF1 *fSlice = new TF1(Form("fSlice_%d_%d", (int)tmin, (int)tmax), fmodel, -1.0, 8.0, 13, "SpectralModel");

      fSlice->SetNpx(1000);

      // Copiar parámetros del global como iniciales
      for (int i = 0; i < 13; i++)
         fSlice->SetParameter(i, fModel->GetParameter(i));

      // Fit libre
      hSlice->Fit(fSlice, "R0");
      gPad->Update(); // necesario para que el TPaveStats exista antes de modificarlo

      TPaveStats *st = (TPaveStats *)hSlice->FindObject("stats");
      if (st) {
         st->SetX1NDC(0.65); // esquina izquierda
         st->SetX2NDC(0.90); // esquina derecha
         st->SetY1NDC(0.35); // esquina inferior
         st->SetY2NDC(0.90); // esquina superior
         st->SetTextSize(0.028);
         gPad->Modified();
      }

      // Dibujar SOLO el fit del slice
      fSlice->SetLineColor(kRed);
      fSlice->Draw("same");

      // Componentes individuales con los parámetros fiteados del slice
      TF1 *fg1 = new TF1(Form("fg1_%d_%d", (int)tmin, (int)tmax), "gaus(0)", Ebin_min, Ebin_max);
      fg1->SetParameters(fSlice->GetParameter(0), fSlice->GetParameter(1), fSlice->GetParameter(2));
      fg1->SetLineColor(kOrange + 7);
      fg1->SetLineStyle(2);
      fg1->SetNpx(1000);
      fg1->Draw("same L");

      TF1 *fg2 = new TF1(Form("fg2_%d_%d", (int)tmin, (int)tmax), "gaus(0)", Ebin_min, Ebin_max);
      fg2->SetParameters(fSlice->GetParameter(3), fSlice->GetParameter(4), fSlice->GetParameter(5));
      fg2->SetLineColor(kBlue);
      fg2->SetLineStyle(2);
      fg2->SetNpx(1000);
      fg2->Draw("same L");

      TF1 *fbw1 =
         new TF1(Form("fbw1_%d_%d", (int)tmin, (int)tmax), "[0]*TMath::BreitWigner(x,[1],[2])", Ebin_min, Ebin_max);
      fbw1->SetParameters(fSlice->GetParameter(6), fSlice->GetParameter(7), fSlice->GetParameter(8));
      fbw1->SetLineColor(kGreen + 2);
      fbw1->SetLineStyle(2);
      fbw1->SetNpx(1000);
      fbw1->Draw("same L");

      TF1 *fbw2 =
         new TF1(Form("fbw2_%d_%d", (int)tmin, (int)tmax), "[0]*TMath::BreitWigner(x,[1],[2])", Ebin_min, Ebin_max);
      fbw2->SetParameters(fSlice->GetParameter(9), fSlice->GetParameter(10), fSlice->GetParameter(11));
      fbw2->SetLineColor(kMagenta);
      fbw2->SetLineStyle(2);
      fbw2->SetNpx(1000);
      fbw2->Draw("same L");

      // Phase space: p12 * graphPS evaluado como TF1
      TF1 *fPS = new TF1(
         Form("fPS_%d_%d", (int)tmin, (int)tmax),
         [graphPS](double *x, double *p) { return p[0] * graphPS->Eval(x[0]); }, Ebin_min, Ebin_max, 1);
      fPS->SetParameter(0, fSlice->GetParameter(12));
      fPS->SetLineColor(kGray + 2);
      fPS->SetLineStyle(2);
      fPS->SetNpx(1000);
      fPS->Draw("same L");

      ipad++;
   }

   cFits->Update();

   // ============================================================
   // 3. Guardar resultados
   // ============================================================
   /* TFile *fOut = new TFile("angular_distributions.root", "RECREATE");
    // Sin corregir
    hAng_GS->Write();
    hAng_1st->Write();
    hAng_2nd->Write();
    hAng_3rd->Write();
    // Corregidos por eficiencia
    hAng_GS_corr->Write();
    hAng_1st_corr->Write();
    hAng_2nd_corr->Write();
    hAng_3rd_corr->Write();
    // 2D para referencia
    h2D->Write();
    fOut->Close();
    cout << "Distribuciones angulares guardadas en angular_distributions.root" << endl;*/
}
