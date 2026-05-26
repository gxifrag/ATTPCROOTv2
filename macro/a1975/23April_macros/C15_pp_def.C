#include <iostream>
#include <fstream>
//#include "frescostuff/fresco_handler.cpp"
#include "penetrabilities_neutron_15C_L_0.C"
#include "penetrabilities_neutron_15C_L_1.C"
#include "penetrabilities_neutron_15C_L_2.C"
#include "penetrabilities_neutron_15C_L_3.C"

// Breit-Wigner standard function
double BreitWigner(double E, double E0, double Gamma) {
    return Gamma / ((E - E0) * (E - E0) + pow(Gamma / 2, 2));
}

// Gaussian function for resolution
double Gaussian(double E, double mean, double sigma) {
    return (1 / (sigma * sqrt(2 * M_PI))) * exp(-0.5 * pow((E - mean) / sigma, 2));
}

// Convoluted Breit-Wigner function
double ConvolutedBW(double *x, double *par, double E_min, double E_max) {
    double E = x[0];
    double Amp = par[0];
    double E0 = par[1];
    double Gamma0 = par[2];
    double sigma = par[3];

    int n_points = 200; // Number of points for numerical integration
    // Limits for the integracion in each point of the function
    double lower = E - sigma; 
    double upper = E + sigma;

    if (lower < E_min) lower = E_min;
    if (upper > E_max) upper = E_max;

    double step = (upper - lower) / n_points; // Steps of each integration for each point
    double integral = 0.0;
// Calculation for each point
    for (int i = 0; i < n_points; ++i) {
        double t = lower + (i + 0.5) * step;
        double bw = BreitWigner(t, E0, Gamma0);
        double gauss = Gaussian(E, t, sigma);
        integral += bw * gauss * step;
    }

    return Amp * integral;
}
// BWModificada function to implement penetrability
double BWModificada(double *x, double *par, int l, int bw_index) {
    static const double E_min_values[4] = {1.22,1.22,1.22,1.22};
    static const double E_max_values[4] = {9,9,9,9};

    if (bw_index < 0 || bw_index >= 4) return 0.0;

    double E = x[0];
    double Amp = par[0];
    double E0 = par[1];
    double Gamma0 = par[2];
    double sigma = par[3];

    const int num_bins_bw = 10000;
    const int num_bins_pen = 10000;
    double bin_width_pen = (E_max_values[bw_index] - E_min_values[bw_index]) / num_bins_pen;

    int bin_index_bw = (E - E_min_values[bw_index]) / 0.000778; // This value is in order to keep 10000 points in the BW
    if (bin_index_bw < 0 || bin_index_bw >= num_bins_bw) return 0.0;

    int bin_index_E0 = static_cast<int>((E0 - E_min_values[bw_index]) / bin_width_pen);
    bin_index_E0 = std::clamp(bin_index_E0, 0, num_bins_pen - 1);

    int bin_index_pen = (E - E_min_values[bw_index]) / bin_width_pen;
    bin_index_pen = std::min(bin_index_pen, num_bins_pen - 1);

    double E_bin = E_min_values[bw_index] + (bin_index_bw + 0.5) * 0.000778;

    double* T[4] = {T0_neutron_15C_values, T1_neutron_15C_values, T2_neutron_15C_values, T3_neutron_15C_values};

    double Gamma_eff, E_eff;

         if (E_bin >= 1.22) {
        double Gamma_val = Gamma0 * T[l][bin_index_pen] / T[l][bin_index_E0];
        double E_val = -Gamma0 * (T[l][bin_index_pen] - T[l][bin_index_E0]) / (2 * T[l][bin_index_E0]);
        Gamma_eff = Gamma_val;
        E_eff = E0 - E_val;
    } 
    else {
        Gamma_eff = Gamma0;
        E_eff = E0;
    }

    double par_conv[4] = {Amp, E_eff, Gamma_eff, sigma};
    double E_min = E_min_values[bw_index];
    double E_max = E_max_values[bw_index];

    return ConvolutedBW(x, par_conv, E_min, E_max);
}


// Function to calculate FWHM of the fit, not implemented in the end
double calculateFWHM(TF1 *func, double E0, double E_min, double E_max) {
    double maxVal = func->Eval(E0);
    double halfMax = maxVal / 2.0;

    // Find the energy values where the function is at half its maximum value
    double E1 = E0, E2 = E0;
    double step = 0.0001; // Step size for searching

    // Search for E1 (left side)
    while (E1 > E_min && func->Eval(E1) > halfMax) {
        E1 -= step;
    }

    // Search for E2 (right side)
    while (E2 < E_max && func->Eval(E2) > halfMax) {
        E2 += step;
    }

    return E2 - E1;
}

Double_t omega(Double_t x, Double_t y, Double_t z)
{
   return sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}


 TGraph* histoToTgraph(TH1F* h) {

    auto g = new TGraph();
    for(int i = 1; i <= h->GetNbinsX(); i++ ) {
        g->SetPoint(i-1,
                    h->GetBinCenter(i),
                    h->GetBinContent(i));
    }

    g->SetName("g_PS_1n");
    return g;
}

   class SpectralModel {
   public:
      TGraph* g_PS;   // Phase-space TGraph
      TGraph* g_PS2;  // Second Phase-space TGraph

      SpectralModel(TGraph* g, TGraph* g2) : g_PS(g), g_PS2(g2) {}

      double operator()(double* x, double* p) {
    double val = 0;

    // Gaussians
    val += p[0] * TMath::Gaus(x[0], p[1], p[2], false);
    val += p[3] * TMath::Gaus(x[0], p[4], p[5], false);

    // Definition of l and bw_index BW
    int l_values[4] = {3, 1, 0, 0};         
    int bw_index_values[4] = {0, 1, 2, 3}; 

double par_bw1[4] = {p[6],  p[7],  p[8],  p[9]};
double par_bw2[4] = {p[10], p[11], p[12], p[13]};
double par_bw3[4] = {p[14], p[15], p[16], p[17]};
double par_bw4[4] = {p[18], p[19], p[20], p[21]};

val += BWModificada(x, par_bw1, l_values[0], bw_index_values[0]);
val += BWModificada(x, par_bw2, l_values[1], bw_index_values[1]);
val += BWModificada(x, par_bw3, l_values[2], bw_index_values[2]);
val += BWModificada(x, par_bw4, l_values[3], bw_index_values[3]);


    // Phase space from TGraph
    val += p[22] * g_PS->Eval(x[0]);
    val += p[23] * g_PS2->Eval(x[0]);

    return val;
}

   };

   std::tuple<double, double>
kine_2b(Double_t m1, Double_t m2, Double_t m3, Double_t m4, Double_t K_proj, Double_t thetalab, Double_t K_eject)
{

   // in this definition: m1(projectile); m2(target); m3(ejectile); and m4(recoil);
   double Et1 = K_proj + m1;
   double Et2 = m2;
   double Et3 = K_eject + m3;
   double Et4 = Et1 + Et2 - Et3;
   double m4_ex, Ex, theta_cm;
   double s, t, u; //---Mandelstam variables

   s = pow(m1, 2) + pow(m2, 2) + 2 * m2 * Et1;
   u = pow(m2, 2) + pow(m3, 2) - 2 * m2 * Et3;

   m4_ex = sqrt((cos(thetalab) * omega(s, pow(m1, 2), pow(m2, 2)) * omega(u, pow(m2, 2), pow(m3, 2)) -
                 (s - pow(m1, 2) - pow(m2, 2)) * (pow(m2, 2) + pow(m3, 2) - u)) /
                   (2 * pow(m2, 2)) +
                s + u - pow(m2, 2));
   Ex = m4_ex - m4;

   t = pow(m2, 2) + pow(m4_ex, 2) - 2 * m2 * Et4;

   // for inverse kinematics Note: this angle corresponds to the recoil
    theta_cm = TMath::Pi() - acos((pow(s, 2) + s * (2 * t - pow(m1, 2) - pow(m2, 2) - pow(m3, 2) - pow(m4_ex, 2)) +
                                  (pow(m1, 2) - pow(m2, 2)) * (pow(m3, 2) - pow(m4_ex, 2))) /
                                 (omega(s, pow(m1, 2), pow(m2, 2)) * omega(s, pow(m3, 2), pow(m4_ex, 2))));

   /*theta_cm = acos((pow(s, 2) + s * (2 * u - pow(m1, 2) - pow(m2, 2) - pow(m3, 2) - pow(m4_ex, 2)) +
                                  (pow(m1, 2) - pow(m2, 2)) * (pow(m4_ex, 2) - pow(m3, 2))) /
                                 (omega(s, pow(m1, 2), pow(m2, 2)) * omega(s, pow(m4_ex, 2), pow(m3, 2))));*/

   theta_cm = theta_cm * TMath::RadToDeg();
   return std::make_tuple(Ex, theta_cm);
}

double Rutherford(double theta, double Ekin, double z, double Z){
    double alpha = 1/137.035999177;
    double hbar = 6.582119569E-22;
    double c = 2.99792458E8;
    return pow(alpha * z * Z, 2) * pow(hbar * c / Ekin, 2) / (16 * pow(TMath::Sin(theta * TMath::DegToRad() / 2), 4)) * 1E28 * 1E3;
}

void GetEnergy(Double_t M, Double_t IZ, Double_t BRO, Double_t &E);
TGraph* ReadKinematics(std::string kineFile);
std::tuple<TGraph*, TGraph*, TGraph*> GetCalculatedXSections(std::string calcFile);


void C15_pp_def()
{

    double Ebin_max = 10.0 ;
	double Ebin_min = -2.0 ;
	int NumberBins = 140; //100
	int NumberBinsAux = 200 ;


    // We define the  histograms that we want to create.
    TH2F *Ang_Ener = new TH2F("Ang_Ener", "Ang_Ener", 1440, 0, 179, 2000, 0, 200.0);
    TH2F *Ang_Ener_PRAC = new TH2F("Ang_Ener_PRAC", "Ang_Ener_PRAC", 1000, 0, 100, 1000, 0, 200.0);
    TH2F *ELossvsBrho = new TH2F("ELossvsBrho", "ELossvsBrho", 4000, 0, 25000, 1000, 0, 4);
    TH2F *dedxvsBrho = new TH2F("dedxvsBrho", "dedxvsBrho", 4000, 0, 10000, 1000, 0, 4);
    TH2F *hVxVy = new TH2F("hVxVy", "hVxVy", 1000, 0, 20, 1000, 0, 20);
    TH1F *henergyIC = new TH1F("henergyIC", "henergyIC", 2048, 0, 2047);

    TH2F *Ang_EnerSim = new TH2F("Ang_EnerSim", "Ang_EnerSim", 1440, 0, 179, 2000, 0, 200.0);
    auto *hexSim = new TH1F("hexSim", "hexSim", 1000, -5, 15);


    auto *hex = new TH1F("hex", "hex", 250, -5, 15);
    auto *hexCoarse = new TH1F("hexCoarse", "hexCoarse", 250, -5, 15);
    auto *QvsEb = new TH2F("QvsEb", "QvsEb", 1000, -5, 15, 300, 0, 300);
    auto *QvsZpos = new TH2F("QvsZpos", "QvsZpos", 1000, -10, 50, 200, -100, 100);

    const double bin_deg = 1.;

    auto *AngDistr = new TH1F("AngDistr", "Ang_Distr", 2048, 0, 120);
    auto *AngDistrCM = new TH1F("AngDistrCM", "^{15}C(p,p)^{15}C elastic scattering", 180 / bin_deg, -0.5, 180-0.5);
    TH1F *AngDistrCMEfficiency = new TH1F("AngDistrCMEfficiency", "AngDistrCMEfficiency", 180 / bin_deg, -0.5, 180-0.5);
    auto *RatioToRutherford = new TH1F("RatioToRutherford", "^{15}C(p,p)^{15}C elastic scattering", 180 / bin_deg, -0.5, 180-0.5);
    auto *ExvsZpos = new TH2F("ExvsZpos", "ExvsZpos", 1000, -5, 15, 300, -0.1, 2);
    auto *thetavsZpos = new TH2F("thetavsZpos", "thetavsZpos", 180, 0, 180, 300, -0.1, 2);

    auto *hredchi2 = new TH1F("redchi2", "redchi2", 1000, 0, 0.0001);
    auto *hbredchi2 = new TH1F("bredchi2", "bredchi2", 1000, 0, 5);

    auto *hex11 = new TH1F("hex11", "hex (63-66deg)", 1000, -5, 55);
    auto *hex12 = new TH1F("hex12", "hex (66-69deg)", 1000, -5, 55);
    auto *hex13 = new TH1F("hex13", "hex (69-72deg)", 1000, -5, 55);
    auto *hex21 = new TH1F("hex21", "hex (72-75deg)", 1000, -5, 55);
    auto *hex22 = new TH1F("hex22", "hex (75-78deg)", 1000, -5, 55);
    auto *hex23 = new TH1F("hex23", "hex (78-81deg)", 1000, -5, 55);
    auto *hex31 = new TH1F("hex31", "hex (81-84deg)", 1000, -5, 55);
    auto *hex32 = new TH1F("hex32", "hex (84-87deg)", 1000, -5, 55);
    auto *hex33 = new TH1F("hex33", "hex (87-90deg)", 1000, -5, 55);

    auto *hexvstheta = new TH2F("hexVStheta", "hexVStheta", 1000, -5, 15, 90, 0, 90);
    auto *hcorrexvstheta = new TH2F("hcorrexVStheta", "hcorrexVStheta", 1000, -5, 15, 90, 0, 90);

    TH2F *hredchi2vsE = new TH2F("hredchi2vsE", "hredchi2vsE", 1000, 0, 0.0001, 1000, 0, 80);
    TH2F *hredchi2vsEx = new TH2F("hredchi2vsEx", "hredchi2vsEx", 1000, 0, 0.0001, 1000, -5, 15);
    TH2F *hredchi2vstheta = new TH2F("hredchi2vstheta", "hredchi2vstheta", 1000, 0, 0.0001, 180, 0, 180);

    std::vector<TH1F*> hHex(18);
    for (int i = 0; i < hHex.size(); i++) {
      hHex[i] = new TH1F(Form("hex%d", i+1), Form("hex%d", i+1), 90, -5, 14);
   }


    // Some useful transformation constants.
    const Double_t u_to_MeV = 931.49401;
    const Double_t Brho_to_p = 1.602176634E-19;
    const Double_t avogadro_number = 6.02214E23;

    // Some masses that may be useful for the experiment.
    Double_t m_p = 1.007825031898 * u_to_MeV;
    Double_t m_d = 2.0135532 * u_to_MeV;
    Double_t m_t = 3.016049281 * u_to_MeV;
    Double_t m_He3 = 3.016029 * u_to_MeV;
    Double_t m_a = 4.00260325415 * u_to_MeV;

    Double_t m_C12 = 12.00 * u_to_MeV;
    Double_t m_C13 = 13.00335484 * u_to_MeV;
    Double_t m_C14 = 14.003242 * u_to_MeV;
    Double_t m_C15 = 15.0105993 * u_to_MeV;
    Double_t m_C16 = 16.0147 * u_to_MeV;
    Double_t m_C17 = 17.0226 * u_to_MeV;

    // Beam and target parameters.
    Double_t Ebeam_buff = 194.5; // MeV
    const Double_t rho_H2 = 3.3084e-5 * avogadro_number / (2 * 1.007825031898); // Atoms/cm3

    // Ejectile parameters.
    int A_ej = 1;
    int Z_ej = 1;
    Double_t m_ej = m_p;
    Double_t m_b = m_p;
    Double_t m_B = m_C15;


    // ELoss tables.
    AtTools::AtELossTable *elossTableCF4 = new AtTools::AtELossTable();
    elossTableCF4->LoadSrimTable("ELossTables/15C_in_CF4_50Torr.txt");

    AtTools::AtELossTable *elossTableH2 = new AtTools::AtELossTable();
    elossTableH2->LoadSrimTable("ELossTables/15C_in_H2_300Torr.txt");

    // ELoss after IC.
    Ebeam_buff = elossTableCF4->GetEnergy(Ebeam_buff, 50);
    std::cout << "Ebeam = " << Ebeam_buff << " MeV.\n";

    // Loading files.
    std::vector<TString> filenames;
    /*for(int i = 138; i <= 182; i++){
     filenames.push_back("run_0" + std::to_string(i));
    }*/

    filenames.push_back("run_0138");
    filenames.push_back("run_0139");
    filenames.push_back("run_0141");
    filenames.push_back("run_0142");
    filenames.push_back("run_0143");
    filenames.push_back("run_0144"); //To check
    filenames.push_back("run_0145");
    filenames.push_back("run_0146");
    filenames.push_back("run_0147");
    filenames.push_back("run_0148");
    filenames.push_back("run_0149");
    filenames.push_back("run_0151");
    filenames.push_back("run_0152");
    filenames.push_back("run_0153");
    filenames.push_back("run_0157");
    filenames.push_back("run_0158");
    filenames.push_back("run_0159");
    filenames.push_back("run_0160");
    filenames.push_back("run_0161");
    filenames.push_back("run_0164"); //To check
    filenames.push_back("run_0165");
    filenames.push_back("run_0166"); //To check
    //filenames.push_back("run_0167");
    filenames.push_back("run_0169"); //To check
    //filenames.push_back("run_0170"); //To be recovered
    filenames.push_back("run_0171");
    filenames.push_back("run_0172");
    filenames.push_back("run_0173");
    filenames.push_back("run_0174");
    filenames.push_back("run_0175");
    filenames.push_back("run_0176");
    filenames.push_back("run_0177"); //To check
    filenames.push_back("run_0178");
    //filenames.push_back("run_0179");
    filenames.push_back("run_0181");
    filenames.push_back("run_0182");



    // Loading simulation files.
    std::vector<TString> filenamesSim;
    filenamesSim.push_back("run_0000");
    filenamesSim.push_back("run_0001");
    filenamesSim.push_back("run_0002");
    filenamesSim.push_back("run_0003");
    filenamesSim.push_back("run_0004");
    filenamesSim.push_back("run_0005");
    filenamesSim.push_back("run_0006");
    filenamesSim.push_back("run_0007");
    filenamesSim.push_back("run_0008");
    filenamesSim.push_back("run_0009");
    filenamesSim.push_back("run_0010");
    filenamesSim.push_back("run_0011");
    filenamesSim.push_back("run_0012");
    filenamesSim.push_back("run_0013");
    filenamesSim.push_back("run_0014");
    filenamesSim.push_back("run_0015");
    filenamesSim.push_back("run_0016");
    filenamesSim.push_back("run_0017");
    filenamesSim.push_back("run_0018");
    filenamesSim.push_back("run_0019");

    // zmin and zmax in meters for the cut on the vertex position.
    const Double_t zmin = 0.25;
    const Double_t zmax = 0.45;

    // Computation of the renormalization constant. Value of rho_H2 obtained from LISE++
    std::cout << "thickness = " << rho_H2 * (zmax - zmin) * 100 << " atoms/cm2.\n";

    // NI will count the total number of ions that entered the ATTPC. We iterate over the scaler files.
    Double_t NI {};
    
    std::cout << "N_I = " << NI << " total ions.\n";
    Double_t normalization_constant = 1 / (2 * TMath::Pi() * rho_H2 * 1E-27 * (zmax - zmin) * 100 * NI * bin_deg * TMath::Pi() / 180); // mbarn
    std::cout << "normalization_constant = "<< normalization_constant << " mbarn.\n";

    // Compute the reconstruction efficiency for the gs from the attpc_engine simulation.
    for(auto filename: filenamesSim){
        TFile *simFile = new TFile("./data/v0.13.0/analysis2/simulation_15C_pp_gs_new/" + filename + "_"  + "1H.root", "R");
        TTree *Tphysics = (TTree *) simFile->Get("parquettree");

        Double_t theta_rad {};
        Double_t phi_rad {};
        Double_t Brho {};
        Double_t redchi {};
        Double_t zPos {};
        Tphysics->SetBranchAddress("polar", &theta_rad);
        Tphysics->SetBranchAddress("azimuthal", &phi_rad);
        Tphysics->SetBranchAddress("brho", &Brho);
        Tphysics->SetBranchAddress("redchisq", &redchi);
        Tphysics->SetBranchAddress("vertex_z", &zPos);

        for(int i = 0; i < Tphysics->GetEntries(); i++){
            Tphysics->GetEntry(i);

            Double_t theta = theta_rad * TMath::RadToDeg();

            Double_t p_ej = Brho * Z_ej * 2.99792458E2; // MeV

            Double_t E_ej = TMath::Sqrt(p_ej * p_ej + m_ej * m_ej) - m_ej; // MeV

            //auto [ex_energy, theta_cm] = kine_2b(m_C15, m_p, m_b, m_B, Ebeam_buff, theta_rad, E_ej);
            auto [ex_energy, theta_cm] = kine_2b(m_C15, m_p, m_b, m_B, elossTableH2->GetEnergy(Ebeam_buff, zPos * 1000), theta_rad, E_ej);

            // Apply cuts.
            if(theta > 90 || theta < 10)
                continue;

            if(zPos < zmin || zPos > zmax)
                continue;

            if(E_ej > 16)
                continue;

            // Ground state cut.
            if(ex_energy < -0.4 || ex_energy > 0.475)
                continue;

            // Histograms.
            Ang_EnerSim->Fill(theta, E_ej);
            hexSim->Fill(ex_energy);
            AngDistrCMEfficiency->Fill(theta_cm, 180./100000. * 1. / (zmax - zmin));
        }
    }

    // Iterate over runs to load kinematic variables of interest.
    for(auto filename: filenames){
       // TFile *physicsFile = new TFile("./data/v0.13.0/analysis2/physics/" + filename + "_" + "1H.root", "R");
       // TFile *physicsFile = new TFile("/home/yassid/attpc_spyral_0.13.0/a2091_15C_pp/InterpSolver_root/" + filename + "_" + "1H.root", "R");
          TFile *physicsFile = new TFile("/home/daniel/Desktop/a2091/InterpSolver_root/" + filename + "_" + "1H.root", "R");
        //TFile *physicsFile = new TFile("/home/yassid/attpc_spyral_dev_test/a2091/C15_pp/InterpSolver_root/" + filename + "_" + "1H.root", "R");
        
        TTree *Tphysics = (TTree *) physicsFile->Get("parquettree");

        Double_t theta_rad {};
        Double_t phi_rad {};
        Double_t Brho {};
        Double_t redchi {};
        Double_t zPos {};
        Tphysics->SetBranchAddress("polar", &theta_rad);
        Tphysics->SetBranchAddress("azimuthal", &phi_rad);
        Tphysics->SetBranchAddress("brho", &Brho);
        Tphysics->SetBranchAddress("redchisq", &redchi);
        Tphysics->SetBranchAddress("vertex_z", &zPos);

        for(int i = 0; i < Tphysics->GetEntries(); i++){
            Tphysics->GetEntry(i);

            Double_t theta = theta_rad * TMath::RadToDeg();

            Double_t p_ej = Brho * Z_ej * 2.99792458E2; // MeV

            Double_t E_ej = TMath::Sqrt(p_ej * p_ej + m_ej * m_ej) - m_ej; // MeV

            //auto [ex_energy, theta_cm] = kine_2b(m_C15, m_p, m_b, m_B, Ebeam_buff, theta_rad, E_ej);
            auto [ex_energy, theta_cm] = kine_2b(m_C15, m_p, m_b, m_B, elossTableH2->GetEnergy(Ebeam_buff, zPos * 1000), theta_rad, E_ej);

            // Apply cuts.
            if(theta > 90 || theta < 10)
                continue;

            if(zPos < zmin || zPos > zmax)
                continue;

            if(E_ej > 16)
                continue;

            // Ground state cut.
           /* if(ex_energy < -0.4 || ex_energy > 0.475)
                continue;*/

            // Excitation energy vs Beam energy
            for (auto iEb = 0; iEb < 300; ++iEb) {
                auto [_ex_energy, _theta_cm] =  kine_2b(m_C15, m_p, m_b, m_B, iEb, theta_rad, E_ej);
                QvsEb->Fill(_ex_energy, iEb);
            }

            // Histograms
            hredchi2->Fill(redchi);
            Ang_Ener->Fill(theta, E_ej);
            hex->Fill(ex_energy);
            hexCoarse->Fill(ex_energy);

            Double_t vx = TMath::Sin(theta_rad) * TMath::Sqrt(E_ej);
            Double_t vy = TMath::Cos(theta_rad) * TMath::Sqrt(E_ej);
            hVxVy->Fill(vx, vy);

            AngDistr->Fill(theta);
            AngDistrCM->Fill(theta_cm, normalization_constant / TMath::Sin(theta_cm * TMath::DegToRad()));
            RatioToRutherford->Fill(theta_cm, normalization_constant / TMath::Sin(theta_cm * TMath::DegToRad()));

            ExvsZpos->Fill(ex_energy, zPos);
            thetavsZpos->Fill(theta, zPos);

            hredchi2vsE->Fill(redchi, E_ej);
            hredchi2vsEx->Fill(redchi, ex_energy);
            hredchi2vstheta->Fill(redchi, theta);

            hexvstheta->Fill(ex_energy, theta_cm);

            for (int i = 0; i < hHex.size(); i++) {
                double theta_min = 10 + i * 5.0;
                double theta_max = theta_min + 5.0;
            if (theta_cm > theta_min && theta_cm <= theta_max) {
                hHex[i]->Fill(ex_energy);
                break; // Only fill one bin per event
            }
         }  


        } // events
    } // Files

    //-------------------- PHASE SPACE ----------------------------------------------

   auto nbins = hex->GetNbinsX() ;
   int binmax = hex->GetMaximumBin() ;
   double ThetaCM_min = 0 ;
   double ThetaCM_max = 180 ;

   TString PhaseSpace_FileName = "./PhaseSpace/PhaseSpace_15C_pp_1n.root" ; 

	TFile *f_PS = new TFile( PhaseSpace_FileName, "READ" ) ; 
	TTree *t_PS = (TTree*)f_PS->Get("simulated_tree") ;

	double Weight_sim, Ex_cal, ThetaCM_cal ;

	t_PS->SetBranchAddress("Weight_sim",&Weight_sim);
	t_PS->SetBranchAddress("Ex_cal",&Ex_cal);
	t_PS->SetBranchAddress("ThetaCM_cal",&ThetaCM_cal);

	TH1F *h_PS_1n = new TH1F("h_PS_1n","h_PS_1n", NumberBins, Ebin_min, Ebin_max ) ; 

   for( int i = 0 ; i < t_PS->GetEntries() ; i++ ) {
		t_PS -> GetEntry(i) ;
		if ( ThetaCM_cal > ThetaCM_min && ThetaCM_cal < ThetaCM_max ) {
			h_PS_1n -> Fill( Ex_cal, Weight_sim ) ;	
		}	
	}
	h_PS_1n -> Smooth() ;


    TString PhaseSpace2n_FileName = "./PhaseSpace/PhaseSpace_15C_pp_2n.root" ; 
    TFile *f_PS2n = new TFile( PhaseSpace2n_FileName, "READ" ) ;
    TTree *t_PS2n = (TTree*)f_PS2n->Get("simulated_tree") ;
    double Weight2n_sim, Ex2n_cal, ThetaCM2n_cal ;
    t_PS2n->SetBranchAddress("Weight_sim",&Weight2n_sim);
    t_PS2n->SetBranchAddress("Ex_cal",&Ex2n_cal);
    t_PS2n->SetBranchAddress("ThetaCM_cal",&ThetaCM2n_cal);
    TH1F *h_PS_2n = new TH1F("h_PS_2n","h_PS_2n", NumberBins, Ebin_min, Ebin_max );
    for( int i = 0 ; i < t_PS2n->GetEntries() ; i++ ) {
          t_PS2n -> GetEntry(i) ;
          if ( ThetaCM2n_cal > ThetaCM_min && ThetaCM2n_cal < ThetaCM_max ) {
                h_PS_2n -> Fill( Ex2n_cal, Weight2n_sim ) ;
          }
     }
    h_PS_2n -> Smooth() ;

    

    //---------------- Fitting the experimental data ----------------//

    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
   
    TGraph* g_PS = histoToTgraph(h_PS_1n);
    TGraph* g_PS2 = histoToTgraph(h_PS_2n);

    // Define a two-gaussian TF1 (ROOT built-in gaus uses amplitude = height)
    TF1 *f2g = new TF1("f2g", "gaus(0) + gaus(3)",-1., 1.4);
    

    // Peak 1
    f2g->SetParameter(1, 0.15);   // mean of gaus(0)
    f2g->SetParameter(2, 0.2);  // sigma guess
    f2g->SetParameter(0, 5000 );  // amplitude guess (1400 - 500 bins)

    // Peak 2
    f2g->SetParameter(4, 0.9);   // mean of gaus(3)
    f2g->SetParameter(5, 0.2);  // sigma guess
    f2g->SetParameter(3, 100);   // amplitude guess (30 - 500 bins)

    hex->Fit(f2g, "R");   // R = use the range you specified, 0 no plot, M minuit
    f2g->SetLineColor(kViolet+2);

    // BW modified

    TF1 *bwprefit1 = new TF1("bw1", [=](double *x, double *par) { 
    return BWModificada(x, par, 3, 0); 
    }, 1.22, 9, 4);  

    TF1 *bwprefit2 = new TF1("bw2", [=](double *x, double *par) { 
    return BWModificada(x, par, 1, 1); 
    }, 1.22, 9, 4);  
    TF1 *bwprefit3 = new TF1("bw3", [=](double *x, double *par) { 
    return BWModificada(x, par, 0, 2); 
    }, 1.22, 9, 4);  

    TF1 *bwprefit4 = new TF1("bw4", [=](double *x, double *par) { 
    return BWModificada(x, par, 0, 3); 
    }, 1.22, 9, 4);  

    /*
    TF1 *bwprefit1 = new TF1("bw1", "[0]*TMath::BreitWigner(x,[1],[2])", 4.3,4.7);
    TF1 *bwprefit2 = new TF1("bw2", "[0]*TMath::BreitWigner(x,[1],[2])", 4.8,5.2);//4.7,5.5 for 500 bins
    TF1 *bwprefit3 = new TF1("bw3", "[0]*TMath::BreitWigner(x,[1],[2])", 6.5, 7.0);
    TF1 *bwprefit4 = new TF1("bw4", "[0]*TMath::BreitWigner(x,[1],[2])", 8.0, 8.8);
*/
    bwprefit1->SetParameters(50,4.5,0.1,0.175); 
    hex->Fit(bwprefit1, "R0");
    bwprefit2->SetParameters(50,5.1,0.02,0.175); 
    hex->Fit(bwprefit2, "R0+");
    bwprefit3->SetParameters(140,6.8,0.1,0.175); 
    hex->Fit(bwprefit3, "R0+");
    bwprefit4->SetParameters(50,8.4,0.1,0.175); 
    hex->Fit(bwprefit4, "R0+");

   double A1    = bwprefit1->GetParameter(0);
   double mean1 = bwprefit1->GetParameter(1);
   double gamma1 = bwprefit1->GetParameter(2);
   double resolution1 = bwprefit1->GetParameter(3);


   double A2    = bwprefit2->GetParameter(0);
   double mean2 = bwprefit2->GetParameter(1);
   double gamma2 = bwprefit2->GetParameter(2);
   double resolution2 = bwprefit2->GetParameter(3);

   double A3    = bwprefit3->GetParameter(0);
   double mean3 = bwprefit3->GetParameter(1);
   double gamma3 = bwprefit3->GetParameter(2);
   double resolution3 = bwprefit3->GetParameter(3);

   double A4    = bwprefit4->GetParameter(0);
   double mean4 = bwprefit4->GetParameter(1);
   double gamma4 = bwprefit4->GetParameter(2);
   double resolution4 = bwprefit4->GetParameter(3);

   SpectralModel* model = new SpectralModel(g_PS,g_PS2);
   TF1* fModel = new TF1("fModel", model, -1., 9.0, 24, "SpectralModel");

   std::vector<double> globalParamsIni = {
    f2g->GetParameter(0), // Amp1
    f2g->GetParameter(1), // Mean1
    f2g->GetParameter(2), // Sigma1
    f2g->GetParameter(3), // Amp2
    f2g->GetParameter(4), // Mean2
    f2g->GetParameter(5), // Sigma2
    A1, mean1, gamma1, resolution1,
    A2, mean2, gamma2, resolution2,
    A3, mean3, gamma3, resolution3,
    A4, mean4, gamma4, resolution4,
    0.001,
    //0.100
};

    // Assigns all the parameters at the same time
for (size_t i = 0; i < globalParamsIni.size(); ++i) {
    fModel->SetParameter(i, globalParamsIni[i]);
    std::cout << "i = " << i << " - globalParamsIni =   " << globalParamsIni[i] << endl;
}

    // Parameter limits and fixed resolution
    fModel->SetParLimits(6, 1, 100);
    fModel->SetParLimits(7, 4.52, 4.63);
    fModel->SetParLimits(8, 0.180, 0.300);
    fModel->FixParameter(9, 0.175);
    fModel->SetParLimits(10, 1, 100);
    fModel->SetParLimits(11, 5.00, 5.10);
    fModel->SetParLimits(12, 0.180, 0.250);
    fModel->FixParameter(13, 0.175);
    fModel->SetParLimits(14, 1, 100);
    fModel->SetParLimits(15, 6.75, 6.85);
    fModel->SetParLimits(16, 0.200, 0.450);
    fModel->FixParameter(17, 0.175);
    fModel->SetParLimits(18, 1, 100);
    fModel->SetParLimits(19, 8.37, 8.50);
    fModel->SetParLimits(20, 0.300, 0.500);
    fModel->FixParameter(21, 0.175);
   
    fModel->SetLineWidth(3);
    fModel->SetLineColor(kOrange+7);
    fModel->SetNpx(1000);
    hex->Fit(fModel, "R");

    int npar = fModel->GetNpar();
   std::vector<double> globalParamsFinals(npar);
   for (int i = 0; i < npar; ++i) {
      globalParamsFinals[i] = fModel->GetParameter(i);
        std::cout << "i = " << i << " - globalParamsFinals =   " << globalParamsFinals[i] << endl;
   }

   for (int i = 0; i < npar; ++i) {
      double v = fModel->GetParameter(i);
      if (!std::isfinite(v)) {
         std::cerr << " - Param[" << i << "] invalid: " << v << "\n";
         return;
      }
   }

   // Gaussian 1
   TF1 *gaus1 = new TF1("gaus1", "gaus(0)", Ebin_min, Ebin_max);
   gaus1->SetParameters(globalParamsFinals[0], globalParamsFinals[1], globalParamsFinals[2]);
   gaus1->SetNpx(1000);
   gaus1->SetLineColor(kOrange+7);
   gaus1->Draw("same");

   // Gaussian 2
   TF1 *gaus2 = new TF1("gaus2", "gaus(0)",Ebin_min, Ebin_max);
   gaus2->SetParameters(globalParamsFinals[3], globalParamsFinals[4], globalParamsFinals[5]);
   gaus2->SetNpx(1000);
   gaus2->SetLineColor(kBlue);
   gaus2->Draw("same");

// Breit-Wigner 1
TF1 *bw1 = new TF1("bw1", [=](double *x, double *par) { 
    return BWModificada(x, par, 0, 0); 
}, 1.22, 9, 4);   
bw1->SetParameters(globalParamsFinals[6], globalParamsFinals[7], globalParamsFinals[8], globalParamsFinals[9]);
bw1->SetNpx(1000);
bw1->SetLineColor(kGreen+2);
bw1->Draw("same");

// Breit-Wigner 2
TF1 *bw2 = new TF1("bw2", [=](double *x, double *par) { 
    return BWModificada(x, par, 1, 1); 
}, 1.22, 9, 4);
bw2->SetParameters(globalParamsFinals[10], globalParamsFinals[11], globalParamsFinals[12], globalParamsFinals[13]);
bw2->SetNpx(1000);
bw2->SetLineColor(kMagenta);
bw2->Draw("same");

// Breit-Wigner 3
TF1 *bw3 = new TF1("bw3", [=](double *x, double *par) { 
    return BWModificada(x, par, 1, 2); 
}, 1.22, 9, 4);
bw3->SetParameters(globalParamsFinals[14], globalParamsFinals[15], globalParamsFinals[16], globalParamsFinals[17]);
bw3->SetLineColor(kRed+2);
bw3->SetNpx(1000);
bw3->Draw("same");

// Breit-Wigner 4
TF1 *bw4 = new TF1("bw4", [=](double *x, double *par) { 
    return BWModificada(x, par, 1, 3); 
}, 1.22, 9, 4);
bw4->SetParameters(globalParamsFinals[18], globalParamsFinals[19], globalParamsFinals[20], globalParamsFinals[21]);
bw4->SetLineColor(kCyan+2);
bw4->SetNpx(1000);
bw4->Draw("same");


   // Phase Space
   if (h_PS_1n->GetNbinsX() > 0) {
      h_PS_1n->Scale(globalParamsFinals[22]);
   }
   h_PS_1n->SetLineColor(kGray+2);
   h_PS_1n->SetLineWidth(2);
   h_PS_1n->Draw("same"); 

    // Phase Space
   if (h_PS_2n->GetNbinsX() > 0) {
      h_PS_2n->Scale(globalParamsFinals[23]);
   }
   h_PS_2n->SetLineColor(kRed);
   h_PS_2n->SetLineWidth(2);
   h_PS_2n->Draw("same"); 

    ///////////// Plots

   TCanvas *c_ExEner = new TCanvas();

   // End Fits
   hex->Draw();
   hex->GetXaxis()->SetTitle("Excitation Energy (MeV)");
   hex->GetYaxis()->SetTitle("Counts");
}
