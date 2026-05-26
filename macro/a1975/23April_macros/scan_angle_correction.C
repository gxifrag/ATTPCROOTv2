//==============================================================================
// scan_angle_correction.C
//
// Standalone macro to find the optimal parameters (alpha, kethe) for the
// empirical angle correction applied in the C16(p,d)C15 analysis:
//
//   theta_corr = theta - alpha * (E_ej - kethe)
//
// The optimal values are found by minimizing the ground state centroid
// offset as a function of (alpha, kethe) on a 2D grid scan.
//
// Usage: root -l scan_angle_correction.C
//
// Output:
//   - Terminal printout of GS centroid for each (alpha, kethe) pair
//   - 2D color map of |centroid| vs (alpha, kethe)
//   - Optimal values printed at the end
//
// Found optimal values (23 April dataset):
//   alpha = 0.00155 rad/MeV
//   kethe = 15.6    MeV
//   -> GS centroid = 0.000056 MeV (~0, essentially perfect)
//
// Author: based on C16_pd_ana_v16_23April.C analysis chain
//==============================================================================

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <tuple>
#include <cmath>

//------------------------------------------------------------------------------
// Helper: Källén function
//------------------------------------------------------------------------------
Double_t omega(Double_t x, Double_t y, Double_t z)
{
    return sqrt(x*x + y*y + z*z - 2*x*y - 2*y*z - 2*x*z);
}

//------------------------------------------------------------------------------
// 2-body kinematics
// m1=projectile, m2=target, m3=ejectile, m4=recoil (all in MeV/c^2)
// K_proj = beam kinetic energy (MeV)
// thetalab = ejectile lab angle (rad)
// K_eject  = ejectile kinetic energy (MeV)
// Returns: {excitation energy (MeV), theta_cm (deg)}
//------------------------------------------------------------------------------
std::tuple<double, double>
kine_2b(Double_t m1, Double_t m2, Double_t m3, Double_t m4,
        Double_t K_proj, Double_t thetalab, Double_t K_eject)
{
    double Et1 = K_proj + m1;
    double Et2 = m2;
    double Et3 = K_eject + m3;
    double Et4 = Et1 + Et2 - Et3;

    double s = pow(m1,2) + pow(m2,2) + 2*m2*Et1;
    double u = pow(m2,2) + pow(m3,2) - 2*m2*Et3;

    double m4_ex = sqrt((cos(thetalab) * omega(s, pow(m1,2), pow(m2,2))
                                       * omega(u, pow(m2,2), pow(m3,2))
                         - (s - pow(m1,2) - pow(m2,2)) * (pow(m2,2) + pow(m3,2) - u))
                        / (2*pow(m2,2))
                        + s + u - pow(m2,2));

    double Ex = m4_ex - m4;

    double t = pow(m2,2) + pow(m4_ex,2) - 2*m2*Et4;

    double theta_cm = TMath::Pi()
        - acos((pow(s,2) + s*(2*t - pow(m1,2) - pow(m2,2) - pow(m3,2) - pow(m4_ex,2))
                + (pow(m1,2) - pow(m2,2)) * (pow(m3,2) - pow(m4_ex,2)))
               / (omega(s, pow(m1,2), pow(m2,2)) * omega(s, pow(m3,2), pow(m4_ex,2))));

    theta_cm *= TMath::RadToDeg();
    return std::make_tuple(Ex, theta_cm);
}

//------------------------------------------------------------------------------
// Main scan function
//------------------------------------------------------------------------------
void scan_angle_correction()
{
    //--------------------------------------------------------------------------
    // USER PARAMETERS — adjust these to match your analysis
    //--------------------------------------------------------------------------

    // Data path
    TString dataPath = "/home/georgina/C16_analysis/C16_H2/"
                       "output_a1975_22April_tb510_MMG20/InterpSolver/"
                       "InterpSolver_pd_root/";

    // Beam energy at AT-TPC entrance (MeV) — from Spyral energy loss chain
    Double_t Ebeam_buff = 186.39;

    // Nuclear masses (MeV/c^2) — nuclear masses, no electron masses
    Double_t m_p   = 938.272076;
    Double_t m_d   = 1875.612931;
    Double_t m_C15 = 13979.218707;
    Double_t m_C16 = 14914.533798;

    // Reaction: 16C(p,d)15C in inverse kinematics
    // m1=16C (beam), m2=p (target), m3=d (ejectile), m4=15C (recoil)
    Double_t m_b  = m_d;    // ejectile
    Double_t m_B  = m_C15;  // recoil
    int      Z_ej = 1;
    Double_t m_ej = m_d;

    // CATIMA energy loss for 16C in H2 gas
    double densityH2 = 3.3084e-5; // g/cm³
    AtTools::AtELossCATIMA elossH2(densityH2);
    double massC16_u = 16.0147;   // u
    elossH2.SetMaterial(catima::Material(1, 1));
    elossH2.SetProjectile(16, 6, massC16_u);

    // Event selection cuts
    double zPos_min_cm = 2.0;   // cm
    double zPos_max_cm = 60.0;  // cm
    double E_ej_max    = 14.0;  // MeV

    // Excitation energy histogram range
    double Ebin_min = -1.0;
    double Ebin_max =  9.0;
    int    NumberBins = 100;

    // GS fit range (MeV) — adjust if needed
    double GS_fit_min = -0.5;
    double GS_fit_max =  0.5;

    //--------------------------------------------------------------------------
    // SCAN GRID — adjust step sizes for coarse/fine scan
    //--------------------------------------------------------------------------

    // Coarse scan (first pass):
    // for (double a = 0.0005; a <= 0.003; a += 0.0002)
    // for (double k = 10.0;   k <= 16.0;  k += 0.5)

    // Fine scan (second pass, around optimal region):
    std::vector<double> alpha_values, kethe_values;
    for (double a = 0.0014; a <= 0.0020; a += 0.00005)
        alpha_values.push_back(a);
    for (double k = 14.0; k <= 16.0; k += 0.2)
        kethe_values.push_back(k);

    //--------------------------------------------------------------------------
    // FILE LIST
    //--------------------------------------------------------------------------
    std::vector<TString> filenames;
    for (int run = 104; run <= 189; run++) {
        // Skip known bad runs
        if (run == 111 || run == 121 || run == 148 || run == 149) continue;
        filenames.push_back(Form("run_%04d_2H.root", run));
    }

    cout << "========================================\n";
    cout << " Angle correction parameter scan\n";
    cout << " Reaction: 16C(p,d)15C\n";
    cout << " Ebeam = " << Ebeam_buff << " MeV\n";
    cout << " Grid: " << alpha_values.size() << " x " << kethe_values.size()
         << " = " << alpha_values.size()*kethe_values.size() << " points\n";
    cout << " Files: " << filenames.size() << " runs\n";
    cout << "========================================\n\n";

    //--------------------------------------------------------------------------
    // 2D SCAN HISTOGRAM
    //--------------------------------------------------------------------------
    TH2F *h_scan = new TH2F("h_scan",
                             "|GS centroid| vs #alpha and kethe;"
                             "#alpha (rad/MeV);kethe (MeV);|centroid| (MeV)",
                             alpha_values.size(), alpha_values.front()-0.000025, alpha_values.back()+0.000025,
                             kethe_values.size(), kethe_values.front()-0.1,      kethe_values.back()+0.1);

    // Also store signed centroid for interpolation
    TH2F *h_scan_signed = new TH2F("h_scan_signed",
                                    "GS centroid (signed) vs #alpha and kethe;"
                                    "#alpha (rad/MeV);kethe (MeV);centroid (MeV)",
                                    alpha_values.size(), alpha_values.front()-0.000025, alpha_values.back()+0.000025,
                                    kethe_values.size(), kethe_values.front()-0.1,      kethe_values.back()+0.1);

    double best_centroid = 999.;
    double best_alpha    = 0.;
    double best_kethe    = 0.;

    //--------------------------------------------------------------------------
    // DOUBLE LOOP OVER (alpha, kethe)
    //--------------------------------------------------------------------------
    for (size_t ia = 0; ia < alpha_values.size(); ia++) {
        for (size_t ik = 0; ik < kethe_values.size(); ik++) {

            double alpha_test = alpha_values[ia];
            double kethe_test = kethe_values[ik];

            TH1F *h_test = new TH1F("h_test", "", NumberBins, Ebin_min, Ebin_max);

            //------------------------------------------------------------------
            // Inner file loop
            //------------------------------------------------------------------
            for (auto& filename : filenames) {
                TFile *f = new TFile(dataPath + filename, "R");
                if (!f || f->IsZombie()) { delete f; continue; }

                TTree *T = (TTree*)f->Get("parquettree");
                if (!T) { f->Close(); delete f; continue; }

                Double_t theta_s{}, Brho_s{}, zPos_s{}, vx_s{}, vy_s{};
                T->SetBranchAddress("polar",    &theta_s);
                T->SetBranchAddress("brho",     &Brho_s);
                T->SetBranchAddress("vertex_z", &zPos_s);
                T->SetBranchAddress("vertex_x", &vx_s);
                T->SetBranchAddress("vertex_y", &vy_s);

                for (int i = 0; i < T->GetEntries(); i++) {
                    T->GetEntry(i);

                    // Ejectile momentum and kinetic energy
                    Double_t p_ej  = Brho_s * Z_ej * 2.99792458 / 10 * 1000;
                    Double_t E_ej  = TMath::Sqrt(p_ej*p_ej + m_ej*m_ej) - m_ej;

                    // Basic cuts
                    if (E_ej  >= E_ej_max)    continue;
                    if (zPos_s*100 < zPos_min_cm || zPos_s*100 > zPos_max_cm) continue;

                    // Beam energy correction (CATIMA, distance in mm)
                    double dist3D  = TMath::Sqrt(vx_s*vx_s + vy_s*vy_s + zPos_s*zPos_s) * 100.0; // cm -> mm
                    Double_t Ebeam = elossH2.GetEnergy(Ebeam_buff, dist3D);

                    // Angle correction
                    double theta_corr = theta_s - alpha_test * (E_ej - kethe_test);

                    // Kinematics
                    auto [ex, tcm] = kine_2b(m_C16, m_p, m_b, m_B, Ebeam, theta_corr, E_ej);
                    h_test->Fill(ex);
                }

                f->Close();
                delete f;
            } // files

            //------------------------------------------------------------------
            // Fit GS peak
            //------------------------------------------------------------------
            TF1 *fg = new TF1("fg", "gaus", GS_fit_min, GS_fit_max);
            fg->SetParameters(h_test->GetMaximum(), 0.0, 0.2);
            int fitStatus = h_test->Fit(fg, "RQ0"); // Q=quiet, 0=don't draw

            double centroid = fg->GetParameter(1);
            double centroid_err = fg->GetParError(1);

            h_scan->SetBinContent(ia+1, ik+1, TMath::Abs(centroid));
            h_scan_signed->SetBinContent(ia+1, ik+1, centroid);

            if (TMath::Abs(centroid) < TMath::Abs(best_centroid)) {
                best_centroid = centroid;
                best_alpha    = alpha_test;
                best_kethe    = kethe_test;
            }

            cout << Form("alpha=%.5f  kethe=%.1f  ->  GS centroid = %+.6f MeV  (err=%.4f)",
                         alpha_test, kethe_test, centroid, centroid_err) << endl;

            delete fg;
            delete h_test;

        } // kethe
    } // alpha

    //--------------------------------------------------------------------------
    // RESULTS
    //--------------------------------------------------------------------------
    cout << "\n========================================\n";
    cout << " SCAN COMPLETE\n";
    cout << "========================================\n";
    cout << Form(" Optimal alpha  = %.5f rad/MeV\n", best_alpha);
    cout << Form(" Optimal kethe  = %.1f MeV\n",     best_kethe);
    cout << Form(" Min |centroid| = %.6f MeV\n",     TMath::Abs(best_centroid));
    cout << "\n Use in main macro:\n";
    cout << Form("   double kethe = %.1f;\n", best_kethe);
    cout << Form("   double theta_corr = theta - %.5f * (E_ej - kethe);\n", best_alpha);
    cout << "========================================\n";

    //--------------------------------------------------------------------------
    // PLOTS
    //--------------------------------------------------------------------------

    // --- 2D scan map ---
    TCanvas *c_scan = new TCanvas("c_scan", "Angle correction parameter scan", 1200, 500);
    c_scan->Divide(2, 1);

    c_scan->cd(1);
    gPad->SetRightMargin(0.18);
    h_scan->Draw("colz");
    h_scan->SetTitle("|GS centroid| vs (#alpha, kethe)");

    // Mark the best point
    TMarker *best_mark = new TMarker(best_alpha, best_kethe, 29); // star
    best_mark->SetMarkerColor(kRed);
    best_mark->SetMarkerSize(2.5);
    best_mark->Draw("same");

    TLatex lat;
    lat.SetNDC(); lat.SetTextSize(0.035);
    lat.DrawLatex(0.15, 0.92, Form("Best: #alpha=%.5f, kethe=%.1f MeV, |centroid|=%.4f MeV",
                                    best_alpha, best_kethe, TMath::Abs(best_centroid)));

    c_scan->cd(2);
    gPad->SetRightMargin(0.18);
    h_scan_signed->Draw("colz");
    h_scan_signed->SetTitle("GS centroid (signed) vs (#alpha, kethe)");

    // Draw zero contour
    h_scan_signed->SetContour(1);
    h_scan_signed->SetContourLevel(0, 0.0);
    TH2F *h_zero = (TH2F*)h_scan_signed->Clone("h_zero");
    h_zero->SetLineColor(kRed);
    h_zero->SetLineWidth(3);
    h_zero->Draw("cont3 same");

    best_mark->Draw("same");

    c_scan->Update();

    // --- 1D projections at optimal values ---
    TCanvas *c_proj = new TCanvas("c_proj", "1D projections at optimal", 1200, 500);
    c_proj->Divide(2, 1);

    // Profile along alpha at best kethe
    c_proj->cd(1);
    int best_ik = h_scan_signed->GetYaxis()->FindBin(best_kethe);
    TH1D *h_alpha_proj = h_scan_signed->ProjectionX("h_alpha_proj", best_ik, best_ik);
    h_alpha_proj->SetTitle(Form("GS centroid vs #alpha (kethe=%.1f MeV)", best_kethe));
    h_alpha_proj->GetXaxis()->SetTitle("#alpha (rad/MeV)");
    h_alpha_proj->GetYaxis()->SetTitle("GS centroid (MeV)");
    h_alpha_proj->SetLineColor(kBlue+1);
    h_alpha_proj->SetLineWidth(2);
    h_alpha_proj->Draw("hist");
    TLine *zero1 = new TLine(h_alpha_proj->GetXaxis()->GetXmin(), 0,
                              h_alpha_proj->GetXaxis()->GetXmax(), 0);
    zero1->SetLineColor(kRed); zero1->SetLineStyle(2); zero1->SetLineWidth(2);
    zero1->Draw("same");

    // Profile along kethe at best alpha
    c_proj->cd(2);
    int best_ia = h_scan_signed->GetXaxis()->FindBin(best_alpha);
    TH1D *h_kethe_proj = h_scan_signed->ProjectionY("h_kethe_proj", best_ia, best_ia);
    h_kethe_proj->SetTitle(Form("GS centroid vs kethe (#alpha=%.5f)", best_alpha));
    h_kethe_proj->GetXaxis()->SetTitle("kethe (MeV)");
    h_kethe_proj->GetYaxis()->SetTitle("GS centroid (MeV)");
    h_kethe_proj->SetLineColor(kBlue+1);
    h_kethe_proj->SetLineWidth(2);
    h_kethe_proj->Draw("hist");
    TLine *zero2 = new TLine(h_kethe_proj->GetXaxis()->GetXmin(), 0,
                              h_kethe_proj->GetXaxis()->GetXmax(), 0);
    zero2->SetLineColor(kRed); zero2->SetLineStyle(2); zero2->SetLineWidth(2);
    zero2->Draw("same");

    c_proj->Update();
}