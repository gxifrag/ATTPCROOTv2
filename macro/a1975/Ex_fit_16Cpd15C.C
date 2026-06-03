#include "Rtypes.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TEventList.h"
#include "TFile.h"
#include "TGraph.h"
#include "TMath.h"
#include "TMinuit.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TString.h"
#include "TText.h"

#include <exception>

#include "iostream"

double ThetaCM_min = 5;
double ThetaCM_max = 150;

double step_ThetaCM = ThetaCM_max - ThetaCM_min;
double ThetaCM_center = 0.5 * (ThetaCM_min + ThetaCM_max);

void PlotHistoGaus(TF1 *func, TH1F *h_plot)
{

   double count_average, step_aux, xi_aux, xi_aux2;

   double binWidth = h_plot->GetBinWidth(1);
   int NumberBins = h_plot->GetNbinsX();

   for (int i = 0; i < NumberBins; i++) {
      xi_aux = h_plot->GetXaxis()->GetBinCenter(i + 1) - 0.5 * binWidth;
      count_average = 0;

      step_aux = binWidth / 20.;

      for (int j = 0; j < 20; j++) {
         xi_aux2 = xi_aux + (0.5 + j) * step_aux;
         count_average += func->Eval(xi_aux2);
      }

      count_average = count_average / 20.;

      // h_plot->Fill(xi_aux + 0.5 * binWidth ,count_average);
      if (count_average < 0.1) {
         h_plot->SetBinContent(i + 1, 0.);
         h_plot->SetBinError(i + 1, 0.001);
      } else {
         h_plot->SetBinContent(i + 1, count_average);
         h_plot->SetBinError(i + 1, sqrt(count_average));
      }
   }
}

// Definitions
TString foutput;
TString finput;
TString yaxis_name;
float yaxis_max;

std::vector<double> x_0;
std::vector<double> y_0;

const int npar = 30;
Double_t par[npar], parerr[npar];
Double_t nc_tot[200];
Double_t nc_Ex1[200];
Double_t nc_Ex2[200];
Double_t nc_Ex3[200];
Double_t nc_Ex4[200];
Double_t nc_Ex5[200];
Double_t nc_Ex6[200];
Double_t nc_Ex7[200];
Double_t nc_Ex8[200];
Double_t nc_PS_1n[200];
Double_t nc_cont_A4[200];
Double_t nc_PS_2n[200];
Double_t nc_GS[200];
Double_t nc_BKG_d[200];
Double_t x[200];
Int_t nbins;
Int_t npt; // number of bins with no zeros

double resolution_from_ThetaCM_GS()
{
   double sigma;
   double p0, p1;
   p0 = 0.0998;
   p1 = 0.004357;
   sigma = p0 + p1 * ThetaCM_center;

   return sigma;
}

double resolution_from_ThetaCM_2plus()
{
   double sigma;
   double p0, p1, p2;
   p0 = 0.183455;
   p1 = -0.00106317;
   p2 = 3.47021e-5;
   sigma = p0 + p1 * ThetaCM_center + p2 * ThetaCM_center * ThetaCM_center;
   return sigma;
}

double resolution_from_ThetaCM_Rest()
{
   double sigma;
   double p0, p1, p2;
   p0 = 0.19786;
   p1 = -0.00205544;
   p2 = 2.5888e-5;
   sigma = p0 + p1 * ThetaCM_center + p2 * ThetaCM_center * ThetaCM_center;

   return sigma;
}

double nc_sum(double *par, int j)
{

   double ncsum = 0;

   int Ninter = 20;

   double deltaX = x[1] - x[0];

   double stepX = deltaX / Ninter;

   double x_start = x[j] - 0.5 * deltaX;

   double xi;

   nc_GS[j] = 0;
   nc_Ex1[j] = 0;
   nc_Ex2[j] = 0;
   nc_Ex3[j] = 0;
   nc_Ex4[j] = 0;
   nc_Ex5[j] = 0;
   nc_Ex6[j] = 0;
   nc_Ex7[j] = 0;
   nc_Ex8[j] = 0;

   // Fixing resolution
   par[18] = resolution_from_ThetaCM_Rest();
   par[21] = resolution_from_ThetaCM_Rest();
   par[24] = resolution_from_ThetaCM_Rest();

   for (int i = 0; i < Ninter; i++) {

      xi = x_start + (i + 0.5) * stepX;

      nc_GS[j] += par[0] * TMath::Gaus(xi, par[1], par[2]);
      nc_Ex1[j] += par[3] * TMath::Gaus(xi, par[4], par[5]);
      nc_Ex2[j] += par[6] * TMath::Gaus(xi, par[7], par[8]);
      nc_Ex3[j] += par[9] * TMath::Gaus(xi, par[10], par[11]);

      nc_Ex4[j] += par[12] * TMath::Voigt(xi - par[13], par[14], par[15], 5);

      nc_Ex5[j] += par[16] * TMath::Gaus(xi, par[17], par[18]);
      nc_Ex6[j] += par[19] * TMath::Gaus(xi, par[20], par[21]);
      nc_Ex7[j] += par[22] * TMath::Gaus(xi, par[23], par[24]);
      nc_Ex8[j] += par[25] * TMath::Gaus(xi, par[26], par[27]);
   }

   nc_GS[j] = nc_GS[j] / Ninter;
   nc_Ex1[j] = nc_Ex1[j] / Ninter;
   nc_Ex2[j] = nc_Ex2[j] / Ninter;
   nc_Ex3[j] = nc_Ex3[j] / Ninter;
   nc_Ex4[j] = nc_Ex4[j] / Ninter;
   nc_Ex5[j] = nc_Ex5[j] / Ninter;
   nc_Ex6[j] = nc_Ex6[j] / Ninter;
   nc_Ex7[j] = nc_Ex7[j] / Ninter;
   nc_Ex8[j] = nc_Ex8[j] / Ninter;

   ncsum = nc_GS[j] + nc_Ex1[j] + nc_Ex2[j] + nc_Ex3[j] + nc_Ex4[j] + nc_Ex5[j] + nc_Ex6[j] + nc_Ex7[j] + nc_Ex8[j] +
           par[28] * nc_PS_1n[j] + par[29] * nc_PS_2n[j];

   return ncsum;
}

void chi2(Int_t &npar, Double_t *gin, Double_t &result, Double_t *par, Int_t flg)
{

   double chisq = 0.0;
   int entry = 0;
   double err_nc = 0.0;
   npt = 0;

   for (int ibin = 0; ibin < nbins; ibin++) {
      double err_nc = 0.0;
      double chi2_bin = 0.0;
      double fac = 0;
      int j = 0;

      if (x[ibin] < 3.5 && x[ibin] > -1.5) {

         if (nc_tot[ibin] == 0)
            err_nc = 1.84;
         if (nc_tot[ibin] == 1 && nc_sum(par, ibin) < 1)
            err_nc = 0.827;
         if (nc_tot[ibin] == 1 && nc_sum(par, ibin) > 1)
            err_nc = 2.300;
         if (nc_tot[ibin] == 2 && nc_sum(par, ibin) < 2)
            err_nc = 1.29;
         if (nc_tot[ibin] == 2 && nc_sum(par, ibin) > 2)
            err_nc = 2.64;

         if (nc_tot[ibin] > 2 && nc_sum(par, ibin) > nc_tot[ibin])
            err_nc = sqrt(nc_tot[ibin]) + 1;
         if (nc_tot[ibin] > 2 && nc_sum(par, ibin) < nc_tot[ibin])
            err_nc = sqrt(nc_tot[ibin]);

         fac = (nc_tot[ibin] - nc_sum(par, ibin));
         chisq += (fac * fac) / (err_nc * err_nc);
         npt++;
      }
   }

   result = chisq;
}

int Ex_fit_16Cpd15C()
{

   double Ebin_max = 15.0;
   double Ebin_min = -5.0;
   int NumberBins = 200;
   int NumberBinsAux = 200;

   TCanvas *c1 = new TCanvas("c1", "histograms", 400, 20, 1000, 1000);
   gROOT->SetStyle("Plain");
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(0);
   gStyle->SetCanvasColor(0);
   gStyle->SetStatColor(0);
   gStyle->SetPadColor(0);
   gStyle->SetTitleFillColor(0);
   gStyle->SetTitleColor(1);
   gStyle->SetPalette(1);
   gStyle->SetCanvasBorderMode(0);
   gStyle->SetPadBorderMode(0);
   gStyle->SetFrameBorderMode(0);

   TH1F *htot, *h_ContA4, *h_ContA4_aux, *htot_plot;

   TH1F *h_PS_1n_plot, *h_PS_2n_plot;

   TF1 *fres_GS, *fres_Ex1, *fres_Ex2, *fres_Ex3, *fres_Ex4, *fres_Ex5, *fres_Ex6, *fres_Ex7, *fres_Ex8;

   TH1F *h_GS_plot = new TH1F("h_GS_plot", "h_GS_plot", NumberBins, Ebin_min, Ebin_max);
   TH1F *h_Ex1_plot = new TH1F("h_Ex1_plot", "h_Ex1_plot", NumberBins, Ebin_min, Ebin_max);
   TH1F *h_Ex2_plot = new TH1F("h_Ex2_plot", "h_Ex2_plot", NumberBins, Ebin_min, Ebin_max);
   TH1F *h_Ex3_plot = new TH1F("h_Ex3_plot", "h_Ex3_plot", NumberBins, Ebin_min, Ebin_max);
   TH1F *h_Ex4_plot = new TH1F("h_Ex4_plot", "h_Ex4_plot", NumberBins, Ebin_min, Ebin_max);
   TH1F *h_Ex5_plot = new TH1F("h_Ex5_plot", "h_Ex5_plot", NumberBins, Ebin_min, Ebin_max);
   TH1F *h_Ex6_plot = new TH1F("h_Ex6_plot", "h_Ex6_plot", NumberBins, Ebin_min, Ebin_max);
   TH1F *h_Ex7_plot = new TH1F("h_Ex7_plot", "h_Ex7_plot", NumberBins, Ebin_min, Ebin_max);
   TH1F *h_Ex8_plot = new TH1F("h_Ex8_plot", "h_Ex8_plot", NumberBins, Ebin_min, Ebin_max);

   TH1F *vFrame;
   TLegend *leg;

   double fitpar[npar], fitparerr[npar];

   htot = new TH1F("htot", "htot", NumberBins, Ebin_min, Ebin_max);

   TString finput_unbinned = "FinalEx_Spyral_Feb13_2025_pd.root";

   TFile *file_data = new TFile(finput_unbinned, "READ");
   TTree *t1 = (TTree *)file_data->Get("yield_tree");

   double Ex_val, ThetaCM_val;

   t1->SetBranchAddress("Ex_val", &Ex_val);
   t1->SetBranchAddress("ThetaCM_val", &ThetaCM_val);

   for (int i = 0; i < t1->GetEntries(); i++) {
      t1->GetEntry(i);
      if (ThetaCM_val > ThetaCM_min && ThetaCM_val < ThetaCM_max) {
         htot->Fill(Ex_val);
      }
   }

   nbins = htot->GetNbinsX();

   int binmax = htot->GetMaximumBin();

   TString PhaseSpace_FileName = "PhaseSpace_16C_pp_1n.root";

   TFile *f_PS = new TFile(PhaseSpace_FileName, "READ");
   TTree *t_PS = (TTree *)f_PS->Get("simulated_tree");

   double Weight_sim, Ex_cal, ThetaCM_cal;

   t_PS->SetBranchAddress("Weight_sim", &Weight_sim);
   t_PS->SetBranchAddress("Ex_cal", &Ex_cal);
   t_PS->SetBranchAddress("ThetaCM_cal", &ThetaCM_cal);

   TH1F *h_PS_1n_aux =
      new TH1F("h_PS_1n_aux", "h_PS_1n_aux", NumberBins, Ebin_min, Ebin_max);      // This is for the 300 keV bin widt
   TH1F *h_PS_1n = new TH1F("h_PS_1n", "h_PS_1n", NumberBins, Ebin_min, Ebin_max); // This is for the 300 keV bin widt

   for (int i = 0; i < t_PS->GetEntries(); i++) {
      t_PS->GetEntry(i);
      if (ThetaCM_cal > ThetaCM_min && ThetaCM_cal < ThetaCM_max) {
         h_PS_1n->Fill(Ex_cal, Weight_sim);
      }
   }

   h_PS_1n->Smooth();

   // 2n phase space
   TString PhaseSpace_2n_FileName = "PhaseSpace_16C_pp_2n.root";

   TFile *f_PS_2n = new TFile(PhaseSpace_2n_FileName, "READ");
   TTree *t_PS_2n = (TTree *)f_PS_2n->Get("simulated_tree");

   double Weight_sim_2n, Ex_cal_2n, ThetaCM_cal_2n;

   t_PS_2n->SetBranchAddress("Weight_sim", &Weight_sim_2n);
   t_PS_2n->SetBranchAddress("Ex_cal", &Ex_cal_2n);
   t_PS_2n->SetBranchAddress("ThetaCM_cal", &ThetaCM_cal_2n);

   TH1F *h_PS_2n_aux =
      new TH1F("h_PS_2n_aux", "h_PS_2n_aux", NumberBins, Ebin_min, Ebin_max);      // This is for the 300 keV bin widt
   TH1F *h_PS_2n = new TH1F("h_PS_2n", "h_PS_2n", NumberBins, Ebin_min, Ebin_max); // This is for the 300 keV bin widt

   for (int i = 0; i < t_PS->GetEntries(); i++) {
      t_PS_2n->GetEntry(i);
      if (ThetaCM_cal_2n > ThetaCM_min && ThetaCM_cal_2n < ThetaCM_max) {
         h_PS_2n->Fill(Ex_cal_2n, Weight_sim_2n);
      }
   }

   h_PS_2n->Smooth();

   // Cast histo data into the arrays for the fitting

   for (Int_t i = 0; i < nbins; i++) {

      x[i] = htot->GetXaxis()->GetBinCenter(i + 1);

      nc_tot[i] = htot->GetBinContent(i + 1);
      nc_PS_1n[i] = h_PS_1n->GetBinContent(i + 1);
      nc_PS_2n[i] = h_PS_2n->GetBinContent(i + 1);
      nc_cont_A4[i] = 0.0; // h_ContA4 -> GetBinContent(i+1) ;

      // cout << nc_cont_A4[i] << endl ;
   }

   //////////////////////////////////////
   // Minimisation
   //////////////////////////////////////

   TStopwatch timer;
   printf("Starting timer\n");
   timer.Start();

   // Number of parameters to minimize
   TMinuit *min = new TMinuit(npar); //  New minimization class for up to 8 variables/parameters

   min->SetFCN(chi2); //  Tell Minuit about the function (above)
   double_t *arglist = new Double_t[npar * 2];
   int ierflg = 1; //  Needed if you want to check for errors

   // Set parameters
   // min->mnparm(0,"first",value,step,min,max,ierr);
   // min->mnparm(1,"second",value,step,min,max,ierr);

   // The mnparm methods define new parameters. The first argument numbers them starting from zero. The second argument
   // gives the parameter a name that you can recognize in the printout that Minuit generates. The other parameters are:

   // value - The initial value of the parameter
   // step - How far Minuit will move away from the initial value when calculating numerical derivatives.
   // min - If the parameter is bounded, this specifies the minimum possible value
   // max - If the parameter is bounded, this specifies the maximum possible value

   char pname[30];
   TString namepar;

   // initial Parameters Set up

   double ip_GS[3] = {1000, 0.0, 0.3};
   double ip_Ex1[3] = {500, 0.77, 0.3};
   double ip_Ex2[3] = {100, 3.0, 0.2};
   double ip_Ex3[3] = {0.1, 9.0, 0.3};
   double ip_Ex4[4] = {0.1, 9.7, 0.3, 0};
   double ip_Ex5[3] = {0.1, 9.5, 0.3};
   double ip_Ex6[3] = {0.1, 9.0, 0.3};
   double ip_Ex7[3] = {0.1, 9.0, 0.3};
   double ip_Ex8[3] = {0.1, 9, 0.3};

   // Fit initial parameters
   double ival[npar];
   double bmin[npar];
   double bmax[npar];
   double step[npar];

   int index_state;

   for (int i = 0; i < 3; i++) {

      index_state = 0;
      ival[i + 3 * index_state] = ip_GS[i];

      index_state = 1;
      ival[i + 3 * index_state] = ip_Ex1[i];

      index_state = 2;
      ival[i + 3 * index_state] = ip_Ex2[i];

      index_state = 3;
      ival[i + 3 * index_state] = ip_Ex3[i];

      index_state = 4;
      ival[i + 3 * index_state] = ip_Ex4[i];

      // We add the one to take into account that Ex4 has 4 index
      index_state = 5;
      ival[i + 3 * index_state + 1] = ip_Ex5[i];

      index_state = 6;
      ival[i + 3 * index_state + 1] = ip_Ex6[i];

      index_state = 7;
      ival[i + 3 * index_state + 1] = ip_Ex7[i];

      index_state = 8;
      ival[i + 3 * index_state + 1] = ip_Ex8[i];
   }

   int N_States = 9;

   for (int k = 0; k < N_States; k++) {
      for (int i = 0; i < 3; i++) {

         int index_mix;

         if (k < 5) {
            index_mix = k * 3 + i;
         } else {
            index_mix = k * 3 + i + 1;
         }

         if (i == 0) {
            bmin[index_mix] = 0;
            bmax[index_mix] = ival[index_mix] * 4;
            if (bmax[index_mix] < 10)
               bmax[index_mix] = 50;
            step[index_mix] = 0.01;
         }
         if (i == 1) {
            bmin[index_mix] = ival[index_mix] - 0.15;
            bmax[index_mix] = ival[index_mix] + 0.15;
            step[index_mix] = 0.001;
         }
         if (i == 2) {
            bmin[index_mix] = 0.15;
            bmax[index_mix] = 0.65;
            // if ( N_States == 7 ) bmax[index_mix] = 0.65 ;
            // else { bmax[index_mix] = 0.35 ; }
            step[index_mix] = 0.001;
         }
      }
   }

   index_state = 4;

   ival[3 * index_state + 3] = 0.00;
   bmin[3 * index_state + 3] = 0.0;
   bmax[3 * index_state + 3] = 0.5;
   step[3 * index_state + 3] = 0.001;

   // PS
   ival[28] = 0.0005;
   bmin[28] = 0.0001;
   bmax[28] = 0.001;
   step[28] = 0.0001;

   // PS
   ival[29] = 0.0; // 0.1 ;
   bmin[29] = 0.0001;
   bmax[29] = 1.5;
   step[29] = 0.0001;

   // ival[22] = 0 ;
   // bmin[22] = 0 ;

   for (int i = 0; i < npar; i++) {
      sprintf(pname, "p%d", i);
      min->mnparm(i, namepar, ival[i], step[i], bmin[i], bmax[i], ierflg);
   }

   arglist[0] = 1;
   min->mnexcm("SET ERR", arglist, 1, ierflg);

   arglist[0] = 50000;
   arglist[1] = 1.;

   min->FixParameter(3 * index_state + 3);
   // min->FixParameter( 5 ) ;

   // Fixing the amplitude of the voigt to 0
   // min->FixParameter( 12 ) ;
   min->FixParameter(15);
   min->FixParameter(23);
   // Fixing the 2n PS to 0
   min->FixParameter(29);
   // min->FixParameter( 28 ) ;

   for (int k = 0; k < N_States - 2; k++) {
      for (int i = 0; i < 3; i++) {

         int index_mix;

         if (k < 5) {
            index_mix = k * 3 + i;
         } else {
            index_mix = k * 3 + i + 1;
         }

         // if ( i == 2 ) min->FixParameter( index_mix ) ;
      }
   }

   // min->FixParameter( 22 ) ;
   // min->FixParameter( 25 ) ;

   cout << "    START: Minimisation procedure" << endl;
   min->mnexcm("MIGRAD", arglist, 2, ierflg);
   cout << "err= " << ierflg << endl;

   // min->Migrad();

   Double_t amin, edm, errdef, dof, prob;
   Int_t nvpar, nparx, icstat;
   min->mnstat(amin, edm, errdef, nvpar, nparx, icstat);
   min->mnprin(3, amin);

   cout << "chi2= " << amin << endl;
   dof = (npt - nvpar);
   cout << "npt= " << (npt) << endl;
   cout << "nvpar= " << (nvpar) << endl;
   cout << "ndof= " << (npt - nvpar) << endl;
   prob = TMath::Prob(amin, dof);
   cout << "prob= " << prob << endl;

   cout << "chi2 = " << amin << endl;
   cout << "n bins / nvpar = " << npt << " / " << nvpar << endl;
   cout << "chi2 / ndf= " << amin / (npt - nvpar) << endl;
   cout << "Npoints = " << npt << endl;
   cout << "Npar    = " << nvpar << endl;

   double covariance_matrix[20][20];
   min->mnemat(&covariance_matrix[0][0], 20);

   for (int i = 0; i < npar; i++) {
      double eplus, eminus, eparab, gcc;
      min->GetParameter(i, par[i], parerr[i]);
      min->mnerrs(i, eplus, eminus, eparab, gcc);
      fitpar[i] = par[i];
      fitparerr[i] = parerr[i];
      cout << "i= " << i << " par=  " << fitpar[i] << "/ parerr=  " << fitparerr[i] << endl;
      cout << "i= " << i << " eplus =  " << eplus << "/ eminus =  " << eminus << "/ eparab =  " << eparab
           << "/ gcc =  " << gcc << endl;
      // cout<<"i= "<< i << sqrt( covariance_matrix[i][i] ) << endl;
   }

   // Parameter print out

   cout << "====================================================================" << endl;
   cout << "===================================================================" << endl;

   cout << " " << endl;

   TString IndexName[20] = {"Ags", "Egs", "Ab1", "Eb1", "Ab2", "Eb2", "Ar6", "Er6", "Amr",   "Emr",
                            "Ar2", "Er2", "Ar3", "Er3", "Ar4", "Er4", "Ar5", "Er5", "A1nPS", "A2nPS"};

   for (int i = 0; i < 10; i++) {
      // cout << " i / j / cov = " << i << " / "  << sqrt( covariance_matrix[i][i] ) << endl;

      for (int j = 0; j < 10; j++) {
         double corr_coef =
            abs(covariance_matrix[i][j] / (sqrt(covariance_matrix[i][i]) * sqrt(covariance_matrix[j][j])));
         if (j > i && corr_coef > 0.2)
            cout << IndexName[i] << " / " << IndexName[j] << " / " << corr_coef << endl;
      }
   }

   cout << "====================================================================" << endl;
   cout << "===================================================================" << endl;

   cout << " " << endl;

   cout << "Parameters of the minimization which lead to the chi2 min" << endl;
   cout << "A GS  = " << fitpar[0] << " +/- " << fitparerr[0] << endl;
   cout << "E GS  = " << fitpar[1] << " +/- " << fitparerr[1] << endl;
   cout << "s GS  = " << fitpar[2] << " +/- " << fitparerr[2] << endl;
   cout << " " << endl;
   cout << "A Ex1 = " << fitpar[3] << " +/- " << fitparerr[3] << endl;
   cout << "E Ex1 = " << fitpar[4] << " +/- " << fitparerr[4] << endl;
   cout << "s Ex1 = " << fitpar[5] << " +/- " << fitparerr[5] << endl;
   cout << " " << endl;
   cout << "A Ex2 = " << fitpar[6] << " +/- " << fitparerr[6] << endl;
   cout << "E Ex2 = " << fitpar[7] << " +/- " << fitparerr[7] << endl;
   cout << "s Ex2 = " << fitpar[8] << " +/- " << fitparerr[8] << endl;
   cout << " " << endl;
   cout << "A Ex3 = " << fitpar[9] << " +/- " << fitparerr[9] << endl;
   cout << "E Ex3 = " << fitpar[10] << " +/- " << fitparerr[10] << endl;
   cout << "s Ex3 = " << fitpar[11] << " +/- " << fitparerr[11] << endl;
   cout << " " << endl;
   cout << "A Ex4 = " << fitpar[12] << " +/- " << fitparerr[12] << endl;
   cout << "E Ex4 = " << fitpar[13] << " +/- " << fitparerr[13] << endl;
   cout << "s Ex4 = " << fitpar[14] << " +/- " << fitparerr[14] << endl;
   cout << "W Ex4 = " << fitpar[15] << " +/- " << fitparerr[15] << endl;
   cout << " " << endl;
   cout << "A Ex5 = " << fitpar[16] << " +/- " << fitparerr[16] << endl;
   cout << "E Ex5 = " << fitpar[17] << " +/- " << fitparerr[17] << endl;
   cout << "s Ex5 = " << fitpar[18] << " +/- " << fitparerr[18] << endl;
   cout << " " << endl;
   cout << "A Ex6 = " << fitpar[19] << " +/- " << fitparerr[19] << endl;
   cout << "E Ex6 = " << fitpar[20] << " +/- " << fitparerr[20] << endl;
   cout << "s Ex6 = " << fitpar[21] << " +/- " << fitparerr[21] << endl;
   cout << " " << endl;
   cout << "A Ex7 = " << fitpar[22] << " +/- " << fitparerr[22] << endl;
   cout << "E Ex7 = " << fitpar[23] << " +/- " << fitparerr[23] << endl;
   cout << "s Ex7 = " << fitpar[24] << " +/- " << fitparerr[24] << endl;
   cout << " " << endl;
   cout << "A Ex8 = " << fitpar[25] << " +/- " << fitparerr[25] << endl;
   cout << "E Ex8 = " << fitpar[26] << " +/- " << fitparerr[26] << endl;
   cout << "s Ex8 = " << fitpar[27] << " +/- " << fitparerr[27] << endl;
   cout << " " << endl;
   cout << "A PS    = " << fitpar[28] << " +/- " << fitparerr[28] << endl;
   cout << " " << endl;
   cout << "A PS 2n = " << fitpar[29] << " +/- " << fitparerr[29] << endl;
   cout << " " << endl;

   cout << "====================================================================" << endl;
   cout << "====================================================================" << endl;

   h_PS_1n_plot = new TH1F("h_PS_1n_plot", "h_PS_1n_plot", NumberBins, Ebin_min, Ebin_max);
   h_PS_2n_plot = new TH1F("h_PS_2n_plot", "h_PS_2n_plot", NumberBins, Ebin_min, Ebin_max);

   fres_GS = new TF1("fres_GS", "[0] * TMath::Gaus( x, [1] , [2] )", -15., 50.);
   fres_Ex1 = new TF1("fres_Ex1", "[0] * TMath::Gaus( x, [1] , [2] )", -15., 50.);
   fres_Ex2 = new TF1("fres_Ex2", "[0] * TMath::Gaus( x, [1] , [2] )", -15., 50.);
   fres_Ex3 = new TF1("fres_Ex3", "[0] * TMath::Gaus( x, [1] , [2] )", -15., 50.);
   fres_Ex4 = new TF1("fres_Ex4", "[0] * TMath::Voigt( x - [1] , [2], [3] )", -15., 50.);
   fres_Ex5 = new TF1("fres_Ex5", "[0] * TMath::Gaus( x, [1] , [2] )", -15., 50.);
   fres_Ex6 = new TF1("fres_Ex6", "[0] * TMath::Gaus( x, [1] , [2] )", -15., 50.);
   fres_Ex7 = new TF1("fres_Ex7", "[0] * TMath::Gaus( x, [1] , [2] )", -15., 50.);
   fres_Ex8 = new TF1("fres_Ex8", "[0] * TMath::Gaus( x, [1] , [2] )", -15., 50.);

   fres_GS->SetNpx(10000);
   fres_Ex1->SetNpx(10000);
   fres_Ex2->SetNpx(10000);
   fres_Ex3->SetNpx(10000);
   fres_Ex4->SetNpx(10000);
   fres_Ex5->SetNpx(10000);
   fres_Ex6->SetNpx(10000);
   fres_Ex7->SetNpx(10000);
   fres_Ex8->SetNpx(10000);

   fres_GS->SetParameters(&fitpar[0]);
   fres_Ex1->SetParameters(&fitpar[3]);
   fres_Ex2->SetParameters(&fitpar[6]);
   fres_Ex3->SetParameters(&fitpar[9]);
   fres_Ex4->SetParameters(&fitpar[12]);
   fres_Ex5->SetParameters(&fitpar[16]);
   fres_Ex6->SetParameters(&fitpar[19]);
   fres_Ex7->SetParameters(&fitpar[22]);
   fres_Ex8->SetParameters(&fitpar[25]);

   PlotHistoGaus(fres_GS, h_GS_plot);
   PlotHistoGaus(fres_Ex1, h_Ex1_plot);
   PlotHistoGaus(fres_Ex2, h_Ex2_plot);
   PlotHistoGaus(fres_Ex3, h_Ex3_plot);
   PlotHistoGaus(fres_Ex4, h_Ex4_plot);
   PlotHistoGaus(fres_Ex5, h_Ex5_plot);
   PlotHistoGaus(fres_Ex6, h_Ex6_plot);
   PlotHistoGaus(fres_Ex7, h_Ex7_plot);
   PlotHistoGaus(fres_Ex8, h_Ex8_plot);

   double Residuo_array[200];
   double Xresiduo_array[200];
   double Residuo_Tot = 0;

   int TotalCounts_PS = 0;

   for (int i = 0; i < NumberBinsAux; i++) {

      double cont = nc_sum(par, i);
      double error_cont = TMath::Sqrt(cont);

      // cout << "x / cont = " << x[i] << " / " << cont << " / " << error_cont << endl ;

      Xresiduo_array[i] = x[i];
      if (nc_tot[i] > 0)
         Residuo_array[i] = 1.0 * (nc_tot[i] - cont) / (sqrt(nc_tot[i]));
      else {
         Residuo_array[i] = 0;
      }

      x_0.push_back(x[i]);
      y_0.push_back(cont);

      Residuo_Tot = Residuo_Tot / nbins;

      double cont_ps_1n = 0;
      cont_ps_1n = par[28] * nc_PS_1n[i];

      h_PS_1n_plot->SetBinContent(i + 1, cont_ps_1n);
      h_PS_1n_plot->SetBinError(i + 1, sqrt(cont_ps_1n));

      double cont_ps_2n = 0;
      cont_ps_2n = par[29] * nc_PS_2n[i];

      h_PS_2n_plot->SetBinContent(i + 1, cont_ps_2n);
      h_PS_2n_plot->SetBinError(i + 1, sqrt(cont_ps_2n));
   }

   TGraph *FitCurve = new TGraph(x_0.size(), &x_0[0], &y_0[0]);

   cout << "flag0" << endl;
   cout << "A GS  = " << fitpar[0] << " +/- " << fitparerr[0] << endl;
   cout << "E GS  = " << fitpar[1] << " +/- " << fitparerr[1] << endl;
   cout << "s GS  = " << fitpar[2] << " +/- " << fitparerr[2] << endl;
   cout << " " << endl;
   cout << "A Ex1 = " << fitpar[3] << " +/- " << fitparerr[3] << endl;
   cout << "E Ex1 = " << fitpar[4] << " +/- " << fitparerr[4] << endl;
   cout << "s Ex1 = " << fitpar[5] << " +/- " << fitparerr[5] << endl;
   cout << " " << endl;
   cout << "A Ex2 = " << fitpar[6] << " +/- " << fitparerr[6] << endl;
   cout << "E Ex2 = " << fitpar[7] << " +/- " << fitparerr[7] << endl;
   cout << "s Ex2 = " << fitpar[8] << " +/- " << fitparerr[8] << endl;
   cout << " " << endl;
   cout << "A Ex3 = " << fitpar[9] << " +/- " << fitparerr[9] << endl;
   cout << "E Ex3 = " << fitpar[10] << " +/- " << fitparerr[10] << endl;
   cout << "s Ex3 = " << fitpar[11] << " +/- " << fitparerr[11] << endl;
   cout << " " << endl;
   cout << "A Ex4 = " << fitpar[12] << " +/- " << fitparerr[12] << endl;
   cout << "E Ex4 = " << fitpar[13] << " +/- " << fitparerr[13] << endl;
   cout << "s Ex4 = " << fitpar[14] << " +/- " << fitparerr[14] << endl;
   cout << "W Ex4 = " << fitpar[15] << " +/- " << fitparerr[15] << endl;
   cout << " " << endl;
   cout << "A Ex5 = " << fitpar[16] << " +/- " << fitparerr[16] << endl;
   cout << "E Ex5 = " << fitpar[17] << " +/- " << fitparerr[17] << endl;
   cout << "s Ex5 = " << fitpar[18] << " +/- " << fitparerr[18] << endl;
   cout << " " << endl;
   cout << "A Ex6 = " << fitpar[19] << " +/- " << fitparerr[19] << endl;
   cout << "E Ex6 = " << fitpar[20] << " +/- " << fitparerr[20] << endl;
   cout << "s Ex6 = " << fitpar[21] << " +/- " << fitparerr[21] << endl;
   cout << " " << endl;
   cout << "A Ex7 = " << fitpar[22] << " +/- " << fitparerr[22] << endl;
   cout << "E Ex7 = " << fitpar[23] << " +/- " << fitparerr[23] << endl;
   cout << "s Ex7 = " << fitpar[24] << " +/- " << fitparerr[24] << endl;
   cout << " " << endl;
   cout << "A Ex8 = " << fitpar[25] << " +/- " << fitparerr[25] << endl;
   cout << "E Ex8 = " << fitpar[26] << " +/- " << fitparerr[26] << endl;
   cout << "s Ex8 = " << fitpar[27] << " +/- " << fitparerr[27] << endl;
   cout << " " << endl;
   cout << "A PS    = " << fitpar[28] << " +/- " << fitparerr[28] << endl;
   cout << " " << endl;
   cout << "A PS 2n = " << fitpar[29] << " +/- " << fitparerr[29] << endl;
   cout << " " << endl;
   /// PLOTTING-------------------/////

   // Plot Options
   gROOT->SetStyle("Plain");
   gStyle->SetOptStat(11111);
   gStyle->SetOptFit(111111);

   gStyle->SetOptStat(00000);
   gStyle->SetOptFit(000000);

   gStyle->SetOptStat(0);
   gStyle->SetOptFit(0);
   gStyle->SetCanvasColor(0);
   gStyle->SetStatColor(0);
   gStyle->SetPadColor(0);
   gStyle->SetTitleFillColor(0);
   gStyle->SetTitleColor(1);
   gStyle->SetPalette(1);
   gStyle->SetCanvasBorderMode(0);
   gStyle->SetPadBorderMode(0);
   gStyle->SetFrameBorderMode(0);

   yaxis_max = 1500.0;

   vFrame = c1->DrawFrame(-1.5, 0.1, 8.5, 500); // yaxis_max);
   vFrame->SetXTitle("E_{x} (MeV)");
   vFrame->GetXaxis()->SetLabelSize(0.035);
   vFrame->GetXaxis()->SetTitleSize(0.045);
   vFrame->GetXaxis()->SetTitleFont(22);

   yaxis_name = "Counts/200 keV";

   vFrame->SetYTitle(yaxis_name);
   // vFrame->SetTitleSize(0.1,axis= "y");
   vFrame->GetYaxis()->SetLabelSize(0.035);
   vFrame->GetYaxis()->SetTitleSize(0.05);
   vFrame->GetXaxis()->SetTitleSize(0.05);
   vFrame->GetYaxis()->SetTitleFont(22);

   gStyle->SetStatColor(0);
   /*
      htot->SetLineColor(1);
      htot->SetLineWidth(3);
      htot->Draw("same hist");
   */

   h_PS_1n_plot->SetLineColor(3);
   h_PS_1n_plot->SetFillColor(3);
   h_PS_1n_plot->SetFillStyle(3001);
   h_PS_1n_plot->SetLineWidth(2);
   h_PS_1n_plot->Draw("hist same");

   h_PS_2n_plot->SetLineColor(4);
   h_PS_2n_plot->SetFillColor(4);
   h_PS_2n_plot->SetFillStyle(3001);
   h_PS_2n_plot->SetLineWidth(2);
   h_PS_2n_plot->Draw("hist same");

   TH1F *h_FitCurve = new TH1F("h_FitCurve", "h_FitCurve", NumberBins, Ebin_min, Ebin_max);
   for (int i = 0; i < x_0.size(); ++i) {
      double x, y;
      FitCurve->GetPoint(i, x, y);
      h_FitCurve->Fill(x, y); // ?
   }

   cout << "flag1" << endl;

   h_GS_plot->SetLineColor(2);
   h_GS_plot->SetLineWidth(3);
   h_GS_plot->SetFillStyle(3001);
   h_GS_plot->SetFillColor(2);
   h_GS_plot->SetLineStyle(1);
   h_GS_plot->Draw("hist same");

   h_Ex1_plot->SetLineColor(3);
   h_Ex1_plot->SetFillStyle(3006);
   h_Ex1_plot->SetFillColor(3);
   h_Ex1_plot->SetLineWidth(2);
   // h_Ex1_plot->SetLineStyle(2);
   h_Ex1_plot->Draw("hist same");

   h_Ex2_plot->SetLineColor(4);
   h_Ex2_plot->SetFillStyle(3007);
   h_Ex2_plot->SetFillColor(4);
   h_Ex2_plot->SetLineWidth(2);
   // h_Ex2_plot->SetLineStyle(2);
   h_Ex2_plot->Draw("hist same");

   h_Ex3_plot->SetLineColor(5);
   h_Ex3_plot->SetFillStyle(3008);
   h_Ex3_plot->SetFillColor(5);
   // h_Ex3_plot->SetLineStyle(2);
   h_Ex3_plot->Draw("hist same");

   h_Ex4_plot->SetLineColor(6);
   h_Ex4_plot->SetFillStyle(3001);
   h_Ex4_plot->SetFillColor(6);
   h_Ex4_plot->SetLineWidth(3);
   h_Ex4_plot->SetLineStyle(1);
   h_Ex4_plot->Draw("hist same");

   h_Ex5_plot->SetLineColor(7);
   h_Ex5_plot->SetFillColor(7);
   h_Ex5_plot->SetFillStyle(3021);
   h_Ex5_plot->SetLineWidth(3);
   // h_Ex5_plot->SetLineStyle(2);
   h_Ex5_plot->Draw("hist same");

   h_Ex6_plot->SetLineColor(6664);
   h_Ex6_plot->SetFillColor(6664);
   h_Ex6_plot->SetFillStyle(3022);
   h_Ex6_plot->SetLineWidth(3);
   // h_Ex6_plot->SetLineStyle(2);
   h_Ex6_plot->Draw("hist same");

   h_Ex7_plot->SetLineColor(6664);
   h_Ex7_plot->SetFillStyle(3609);
   h_Ex7_plot->SetLineWidth(2);
   // h_Ex7_plot->SetLineStyle(2);
   h_Ex7_plot->Draw("hist same");

   h_Ex8_plot->SetLineColor(6664);
   h_Ex8_plot->SetLineWidth(2);
   // h_Ex8_plot->SetLineStyle(2);
   h_Ex8_plot->Draw("hist same");

   h_FitCurve->SetLineColor(1);
   h_FitCurve->SetMarkerColor(1);
   h_FitCurve->SetLineWidth(3);
   h_FitCurve->SetMarkerStyle(21);

   h_FitCurve->Draw("hist same");

   htot->SetLineColor(1);
   htot->SetLineWidth(3);
   htot->Draw("esame");

   cout << "flag2" << endl;

   // leg = new TLegend(0.6,0.65,0.95,0.90);
   leg = new TLegend(0.57, 0.65, 0.75, 0.95);
   leg->SetFillColor(0);
   leg->SetBorderSize(0);
   leg->SetTextSize(0.03);
   leg->SetTextFont(102);
   leg->AddEntry(htot, "Experimental Data", "le");
   leg->AddEntry(h_FitCurve, "Total fit", "l");
   leg->AddEntry(h_GS_plot, "G.S. 0p_{1/2}", "f");
   leg->AddEntry(h_Ex1_plot, "E_{x}=1.21 MeV, 0p_{3/2}", "f");
   leg->AddEntry(h_Ex2_plot, "E_{x}=2.58 MeV, 0p_{3/2}", "f");
   leg->AddEntry(h_Ex3_plot, "E_{x}=3.85 MeV, 0p_{3/2}", "f");
   leg->AddEntry(h_Ex4_plot, "E_{x}=5.30 MeV, 0p_{3/2}", "f");
   leg->AddEntry(h_Ex5_plot, "E_{x}=7.29 MeV, 0p_{3/2}", "f");
   leg->AddEntry(h_Ex6_plot, "E_{x}=9.39 MeV, 0p_{3/2}", "f");
   leg->AddEntry(h_PS_1n_plot, "Phase Space 1n", "f");
   // leg->Draw();

   c1->Update();

   gPad->SetTopMargin(0.03);
   gPad->SetBottomMargin(0.1);
   gPad->SetRightMargin(0.03);
   gPad->SetLeftMargin(0.1);
   printf("Time at the end of job = %f seconds\n", timer.CpuTime());

   TH1F *h_GS_plot_save = new TH1F("h_GS_plot_save", "h_GS_plot_save", NumberBins, Ebin_min, Ebin_max);
   TH1F *h_Ex1_plot_save = new TH1F("h_Ex1_plot_save", "h_Ex1_plot_save", NumberBins, Ebin_min, Ebin_max);
   TH1F *h_Ex2_plot_save = new TH1F("h_Ex2_plot_save", "h_Ex2_plot_save", NumberBins, Ebin_min, Ebin_max);
   TH1F *h_Ex3_plot_save = new TH1F("h_Ex3_plot_save", "h_Ex3_plot_save", NumberBins, Ebin_min, Ebin_max);
   TH1F *h_Ex4_plot_save = new TH1F("h_Ex4_plot_save", "h_Ex4_plot_save", NumberBins, Ebin_min, Ebin_max);
   TH1F *h_Ex5_plot_save = new TH1F("h_Ex5_plot_save", "h_Ex5_plot_save", NumberBins, Ebin_min, Ebin_max);
   TH1F *h_Ex6_plot_save = new TH1F("h_Ex6_plot_save", "h_Ex6_plot_save", NumberBins, Ebin_min, Ebin_max);

   TH1F *h_FitCurve_save = new TH1F("h_FitCurve_save", "h_FitCurve_save", NumberBins, Ebin_min, Ebin_max);
   TH1F *h_PS_1n_plot_save = new TH1F("h_PS_1n_plot_save", "h_PS_1n_plot_save", NumberBins, Ebin_min, Ebin_max);

   for (int i = 0; i < NumberBins; ++i) {
      if (h_GS_plot->GetBinContent(i + 1) > 0) {
         h_GS_plot_save->SetBinContent(i + 1, h_GS_plot->GetBinContent(i + 1));
      } else {
         h_GS_plot_save->SetBinContent(i + 1, 0);
      }

      if (h_Ex1_plot->GetBinContent(i + 1) > 0) {
         h_Ex1_plot_save->SetBinContent(i + 1, h_Ex1_plot->GetBinContent(i + 1));
      } else {
         h_Ex1_plot_save->SetBinContent(i + 1, 0);
      }

      if (h_Ex2_plot->GetBinContent(i + 1) > 0) {
         h_Ex2_plot_save->SetBinContent(i + 1, h_Ex2_plot->GetBinContent(i + 1));
      } else {
         h_Ex2_plot_save->SetBinContent(i + 1, 0);
      }

      if (h_Ex3_plot->GetBinContent(i + 1) > 0) {
         h_Ex3_plot_save->SetBinContent(i + 1, h_Ex3_plot->GetBinContent(i + 1));
      } else {
         h_Ex3_plot_save->SetBinContent(i + 1, 0);
      }

      if (h_Ex4_plot->GetBinContent(i + 1) > 0) {
         h_Ex4_plot_save->SetBinContent(i + 1, h_Ex4_plot->GetBinContent(i + 1));
      } else {
         h_Ex4_plot_save->SetBinContent(i + 1, 0);
      }

      if (h_Ex5_plot->GetBinContent(i + 1) > 0) {
         h_Ex5_plot_save->SetBinContent(i + 1, h_Ex5_plot->GetBinContent(i + 1));
      } else {
         h_Ex5_plot_save->SetBinContent(i + 1, 0);
      }

      if (h_Ex6_plot->GetBinContent(i + 1) > 0) {
         h_Ex6_plot_save->SetBinContent(i + 1, h_Ex6_plot->GetBinContent(i + 1));
      } else {
         h_Ex6_plot_save->SetBinContent(i + 1, 0);
      }

      if (h_FitCurve->GetBinContent(i + 1) > 0) {
         h_FitCurve_save->SetBinContent(i + 1, h_FitCurve->GetBinContent(i + 1));
      } else {
         h_FitCurve_save->SetBinContent(i + 1, 0);
      }

      if (h_PS_1n_plot->GetBinContent(i + 1) > 0) {
         h_PS_1n_plot_save->SetBinContent(i + 1, h_PS_1n_plot->GetBinContent(i + 1));
      } else {
         h_PS_1n_plot_save->SetBinContent(i + 1, 0);
      }
   }

   /*
      //Save histograms folder
      TFile *f_Histos19N = new TFile("/home/jlf-general/Desktop/19N_paper_files/Ex_and_fitted_states.root","recreate") ;
      htot -> Write() ;
      h_GS_plot_save    -> Write() ;
      h_Ex1_plot_save   -> Write() ;
      h_Ex2_plot_save   -> Write() ;
      h_Ex3_plot_save   -> Write() ;
      h_Ex4_plot_save   -> Write() ;
      h_Ex5_plot_save   -> Write() ;
      h_Ex6_plot_save   -> Write() ;
      h_FitCurve_save   -> Write() ;
      h_PS_1n_plot_save -> Write() ;
      fres_GS  -> Write() ;
      fres_Ex1 -> Write() ;
      fres_Ex2 -> Write() ;
      fres_Ex3 -> Write() ;
      fres_Ex4 -> Write() ;
      fres_Ex5 -> Write() ;
      fres_Ex6 -> Write() ;
      f_Histos19N -> Close() ;
   */

   cout << "chi2= " << amin << endl;

   double int_GS, int_Ex1, int_Ex2, int_Ex3, int_Ex4, int_Ex5, int_Ex6, int_Ex7, int_Ex8, int_Ex9;
   double u_int_GS, u_int_Ex1, u_int_Ex2, u_int_Ex3, u_int_Ex4, u_int_Ex5, u_int_Ex6, u_int_Ex7, u_int_Ex8, u_int_Ex9;

   int_GS = fres_GS->Integral(-20, 20) / htot->GetBinWidth(1);
   u_int_GS = (fitparerr[0] / fitpar[0]) * fres_GS->Integral(-20, 20) / htot->GetBinWidth(1);

   int_Ex1 = fres_Ex1->Integral(-20, 20) / htot->GetBinWidth(1);
   u_int_Ex1 = (fitparerr[3] / fitpar[3]) * fres_Ex1->Integral(-20, 20) / htot->GetBinWidth(1);

   int_Ex2 = fres_Ex2->Integral(-20, 20) / htot->GetBinWidth(1);
   u_int_Ex2 = (fitparerr[6] / fitpar[6]) * fres_Ex2->Integral(-20, 20) / htot->GetBinWidth(1);

   int_Ex3 = fres_Ex3->Integral(-20, 20) / htot->GetBinWidth(1);
   u_int_Ex3 = (fitparerr[9] / fitpar[9]) * fres_Ex3->Integral(-20, 20) / htot->GetBinWidth(1);

   int_Ex4 = fres_Ex4->Integral(-20, 20) / htot->GetBinWidth(1);
   u_int_Ex4 = (fitparerr[12] / fitpar[12]) * fres_Ex4->Integral(-20, 20) / htot->GetBinWidth(1);

   int_Ex5 = fres_Ex5->Integral(-20, 20) / htot->GetBinWidth(1);
   u_int_Ex5 = (fitparerr[16] / fitpar[16]) * fres_Ex5->Integral(-20, 20) / htot->GetBinWidth(1);

   int_Ex6 = fres_Ex6->Integral(-20, 20) / htot->GetBinWidth(1);
   u_int_Ex6 = (fitparerr[19] / fitpar[19]) * fres_Ex6->Integral(-20, 20) / htot->GetBinWidth(1);

   int_Ex7 = fres_Ex7->Integral(-20, 20) / htot->GetBinWidth(1);
   u_int_Ex7 = (fitparerr[22] / fitpar[22]) * fres_Ex7->Integral(-20, 20) / htot->GetBinWidth(1);

   int_Ex8 = fres_Ex8->Integral(-20, 20) / htot->GetBinWidth(1);
   u_int_Ex8 = (fitparerr[25] / fitpar[25]) * fres_Ex8->Integral(-20, 20) / htot->GetBinWidth(1);

   // int_Ex9   =  fres_Ex3 -> Integral( -20 , 20)/ htot -> GetBinWidth(1) ;
   // u_int_Ex9 = ( fitparerr[6] / fitpar[6] ) * fres_Ex3 -> Integral( -20, 20 ) / htot->GetBinWidth(1) ;

   // Alternative Ncounts integrated for gs

   // Find the corresponding bin numbers for the X limits
   int gs_min = htot->FindBin(fitpar[1] - 2 * fitpar[2]);
   int gs_max = htot->FindBin(fitpar[1] + 2 * fitpar[2]);

   int Ex1_min = htot->FindBin(fitpar[4] - 2 * fitpar[5]);
   int Ex1_max = htot->FindBin(fitpar[4] + 2 * fitpar[5]);

   int Ex2_min = htot->FindBin(fitpar[7] - 2 * fitpar[8]);
   int Ex2_max = htot->FindBin(fitpar[7] + 2 * fitpar[8]);

   int Ex3_min = htot->FindBin(fitpar[10] - 2 * fitpar[11]);
   int Ex3_max = htot->FindBin(fitpar[10] + 2 * fitpar[11]);

   // Integrate the histogram between bin_min and bin_max
   double integral_gs = htot->Integral(gs_min, gs_max);
   double integral_Ex1 = htot->Integral(Ex1_min, Ex1_max);
   double integral_Ex2 = htot->Integral(Ex2_min, Ex2_max);
   double integral_Ex3 = htot->Integral(Ex3_min, Ex3_max);

   // Amass selection correction

   int_GS = int_GS;
   u_int_GS = u_int_GS;

   int_Ex1 = int_Ex1;
   u_int_Ex1 = u_int_Ex1;

   int_Ex2 = int_Ex2;
   u_int_Ex2 = u_int_Ex2;

   int_Ex4 = int_Ex4;
   u_int_Ex4 = u_int_Ex4;

   int_Ex5 = int_Ex5;
   u_int_Ex5 = u_int_Ex5;

   int_Ex6 = int_Ex6;
   u_int_Ex6 = u_int_Ex6;

   cout << "-----------------------------------------------------------------------------------------------------------"
           "---------------------------------------------------------------------------------------------"
        << endl;

   cout << "   " << endl;
   cout << "   " << endl;
   /*
   cout << TString::Format("%.2f", ThetaCM_center)  << "   " << TString::Format("%.2f", step_ThetaCM ) << "   " <<
   TString::Format("%.4f", int_GS ) << "   " << TString::Format("%.4f", u_int_GS ) << "   " << TString::Format("%.4f",
   int_Ex1 ) << "   " << TString::Format("%.4f", u_int_Ex1 ) << "   " << TString::Format("%.4f", int_Ex2 ) << "   " <<
   TString::Format("%.4f", u_int_Ex2 ) << "   " << TString::Format("%.4f", int_Ex3) << "   " << TString::Format("%.4f",
   u_int_Ex3) << "   " << TString::Format("%.4f", int_Ex5 )<< "   " << TString::Format("%.4f", u_int_Ex5 ) << "   " <<
   TString::Format("%.4f", int_Ex6 ) << "   " << TString::Format("%.4f", u_int_Ex6 ) << "   " << TString::Format("%.4f",
   int_Ex7 ) << "   " << TString::Format("%.4f", u_int_Ex7 ) << endl ;
   */

   cout << TString::Format("%.2f", ThetaCM_center) << "   " << TString::Format("%.2f", step_ThetaCM) << "   "
        << TString::Format("%.4f", int_GS) << "   " << TString::Format("%.4f", u_int_GS) << "   "
        << TString::Format("%.4f", int_Ex1) << "   " << TString::Format("%.4f", u_int_Ex1) << "   "
        << TString::Format("%.4f", int_Ex2) << "   " << TString::Format("%.4f", u_int_Ex2) << "   "
        << TString::Format("%.4f", int_Ex3) << "   " << TString::Format("%.4f", u_int_Ex3) << "   " << endl;

   cout << "   " << endl;
   cout << "   " << endl;
   cout << "-----------------------------------------------------------------------------------------------------------"
           "---------------------------------------------------------------------------------------------"
        << endl;

   cout << " " << endl;
   cout << " " << endl;
   cout << " " << endl;

   cout << "-----------------------------------------------------------------------------------------------------------"
           "----------------------------------------------------------------------------------------"
        << endl;

   cout << TString::Format("%.2f", ThetaCM_center) << "   " << TString::Format("%.2f", step_ThetaCM) << "   "
        << TString::Format("%.2f", fitpar[1]) << "   " << TString::Format("%.2f", fitpar[2]) << "   "
        << TString::Format("%.2f", fitpar[4]) << "   " << TString::Format("%.2f", fitpar[5]) << "   "
        << TString::Format("%.2f", fitpar[7]) << "   " << TString::Format("%.2f", fitpar[8]) << "   "
        << TString::Format("%.2f", fitpar[10]) << "   " << TString::Format("%.2f", fitpar[11]) << "   "
        << TString::Format("%.2f", fitpar[17]) << "   " << TString::Format("%.2f", fitpar[18]) << "   "
        << TString::Format("%.2f", fitpar[20]) << "   " << TString::Format("%.2f", fitpar[21]) << "   "
        << TString::Format("%.2f", fitpar[23]) << "   " << TString::Format("%.2f", fitpar[24]) << endl;
   cout << "   " << endl;
   cout << "   " << endl;
   cout << "-----------------------------------------------------------------------------------------------------------"
           "----------------------------------------------------------------------------------------"
        << endl;

   cout << int_GS << "   " << integral_gs << endl;
   cout << int_Ex1 << "   " << integral_Ex1 << endl;
   cout << int_Ex2 << "   " << integral_Ex2 << endl;
   // cout << int_Ex3 << "   " << integral_Ex3 << endl ;

   return 0;
}
