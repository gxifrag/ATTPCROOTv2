#include <iostream>
#include <fstream>
//#include "frescostuff/fresco_handler.cpp"

Bool_t compareEventName(std::string &getname, std::string &fribname)
{
   // Parsing FRIB event number
   std::regex fribregex("evt(\\d+)_\\d+");
   std::string result = std::regex_replace(fribname, fribregex, "$1\n");

   int fribnumber;
   std::istringstream iss(result);
   while (iss >> fribnumber) {
      // std::cout << fribnumber << std::endl;
   }

   // Parsing GET event name
   std::regex getregex("evt(\\d+)_data");
   result = std::regex_replace(getname, getregex, "$1\n");

   int getnumber;
   std::istringstream isss(result);
   while (isss >> getnumber) {
      // std::cout << getnumber << std::endl;
   }

   return (fribnumber == getnumber) ? 1 : 0;

   return 0;
}

Double_t omega(Double_t x, Double_t y, Double_t z)
{
   return sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}

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

   //Present case refers to a scattering process in the lab frame
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


void C15_dd_ana()
{
    TCutG *cutg = new TCutG("CUTG",18);
   cutg->SetVarX("Ang_Ener");
   cutg->SetVarY("");
   cutg->SetTitle("Graph");
   cutg->SetFillStyle(1000);
   cutg->SetPoint(0,43.0606,34.5727);
   cutg->SetPoint(1,41.724,33.5584);
   cutg->SetPoint(2,42.7742,29.248);
   cutg->SetPoint(3,44.461,27.283);
   cutg->SetPoint(4,49.4258,19.4863);
   cutg->SetPoint(5,53.0539,14.1617);
   cutg->SetPoint(6,55.2499,11.2458);
   cutg->SetPoint(7,55.8864,6.99881);
   cutg->SetPoint(8,57.4459,5.92121);
   cutg->SetPoint(9,57.9232,8.32996);
   cutg->SetPoint(10,56.4911,13.6546);
   cutg->SetPoint(11,53.5313,18.3453);
   cutg->SetPoint(12,49.5531,24.0502);
   cutg->SetPoint(13,46.4342,29.0579);
   cutg->SetPoint(14,43.9836,33.812);
   cutg->SetPoint(15,43.2198,34.6994);
   cutg->SetPoint(16,42.997,34.4459);
   cutg->SetPoint(17,43.0606,34.5727);

    std::cout << "C15_dd_ana started.\n";
    // We define the  histograms that we want to create.
    TH2F *Ang_Ener = new TH2F("Ang_Ener", "Ang_Ener", 720, 0, 90, 900, 0, 90.0);
    TH2F *Ang_Ener_PRAC = new TH2F("Ang_Ener_PRAC", "Ang_Ener_PRAC", 1000, 0, 100, 1000, 0, 200.0);
    TH2F *ELossvsBrho = new TH2F("ELossvsBrho", "ELossvsBrho", 4000, 0, 25000, 1000, 0, 4);
    TH2F *dedxvsBrho = new TH2F("dedxvsBrho", "dedxvsBrho", 4000, 0, 10000, 1000, 0, 4);
    TH2F *hVxVy = new TH2F("hVxVy", "hVxVy", 1000, 0, 20, 1000, 0, 20);
    TH1F *henergyIC = new TH1F("henergyIC", "henergyIC", 2048, 0, 2047);

    /*TH2F *Ang_EnerSim = new TH2F("Ang_EnerSim", "Ang_EnerSim", 1440, 0, 179, 2000, 0, 200.0);
    auto *hexSim = new TH1F("hexSim", "hexSim", 1000, -5, 15);*/


    auto *hex = new TH1F("hex", "hex", 1000, -5, 15);
    auto *QvsEb = new TH2F("QvsEb", "QvsEb", 1000, -5, 15, 300, 0, 300);
    auto *QvsZpos = new TH2F("QvsZpos", "QvsZpos", 1000, -10, 50, 200, -100, 100);

    const double bin_deg = 1.;

    auto *AngDistr = new TH1F("AngDistr", "Ang_Distr", 2048, 0, 120);
    auto *AngDistrCM = new TH1F("AngDistrCM", "^{15}C(d,d)^{15}C elastic scattering", 180 / bin_deg, -0.5, 180-0.5);
    TH1F *AngDistrCMEfficiency = new TH1F("AngDistrCMEfficiency", "AngDistrCMEfficiency", 180 / bin_deg, -0.5, 180-0.5);
   // auto *RatioToRutherford = new TH1F("RatioToRutherford", "^{15}C(d,d)^{15}C elastic scattering", 180 / bin_deg, -0.5, 180-0.5);
    auto *ExvsZpos = new TH2F("ExvsZpos", "ExvsZpos", 1000, -5, 15, 300, -0.1, 2);
    auto *thetavsZpos = new TH2F("thetavsZpos", "thetavsZpos", 100, 0, 100, 300, 0, 1);

    /*auto *hredchi2 = new TH1F("redchi2", "redchi2", 1000, 0, 0.0001);
    auto *hbredchi2 = new TH1F("bredchi2", "bredchi2", 1000, 0, 5);*/

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

    /*TH2F *hredchi2vsE = new TH2F("hredchi2vsE", "hredchi2vsE", 1000, 0, 0.0001, 1000, 0, 80);
    TH2F *hredchi2vsEx = new TH2F("hredchi2vsEx", "hredchi2vsEx", 1000, 0, 0.0001, 1000, -5, 15);
    TH2F *hredchi2vstheta = new TH2F("hredchi2vstheta", "hredchi2vstheta", 1000, 0, 0.0001, 180, 0, 180);*/


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
    const Double_t rho_D2 = 6.5643e-05 * avogadro_number / (2 * 2.0135532); // Atoms/cm3

    // Ejectile parameters: en el meu cas, deuteró
    int A_ej = 2;
    int Z_ej = 1;
    Double_t m_ej = m_d;
    Double_t m_b = m_d;
    Double_t m_B = m_C15;


    // ELoss tables.
    /*AtTools::AtELossTable *elossTableCF4 = new AtTools::AtELossTable();
    elossTableCF4->LoadSrimTable("ELossTables/15C_in_CF4_50Torr.txt");

    AtTools::AtELossTable *elossTableH2 = new AtTools::AtELossTable();
    elossTableH2->LoadSrimTable("ELossTables/15C_in_H2_300Torr.txt");

    // ELoss after IC.
    Ebeam_buff = elossTableCF4->GetEnergy(Ebeam_buff, 50);
    std::cout << "Ebeam = " << Ebeam_buff << " MeV.\n";*/

    // Loading files.
    std::vector<TString> filenames; 
std::cout << "Loading files...\n";
    //filenames.push_back("run_0013");
    filenames.push_back("run_0014");
    filenames.push_back("run_0015");
    filenames.push_back("run_0016");
    filenames.push_back("run_0017");
    //filenames.push_back("run_0018"); // excluded
    //filenames.push_back("run_0019"); // excluded
    filenames.push_back("run_0020");
    filenames.push_back("run_0021");
    filenames.push_back("run_0022");
    filenames.push_back("run_0023");
    filenames.push_back("run_0024");
    //filenames.push_back("run_0025"); // excluded
    filenames.push_back("run_0026");
    //filenames.push_back("run_0027"); // excluded
    filenames.push_back("run_0028");
    filenames.push_back("run_0029");
    filenames.push_back("run_0030");
    //filenames.push_back("run_0031"); // excluded
    //filenames.push_back("run_0032"); // excluded
    filenames.push_back("run_0033");
    filenames.push_back("run_0034");
    //filenames.push_back("run_0035"); // excluded
    filenames.push_back("run_0036");
    filenames.push_back("run_0037");
    //filenames.push_back("run_0038"); // excluded
    //filenames.push_back("run_0039"); // excluded
    //filenames.push_back("run_0040"); // excluded
    //filenames.push_back("run_0041"); // excluded
    filenames.push_back("run_0042");
    filenames.push_back("run_0043");
    filenames.push_back("run_0044");
    //filenames.push_back("run_0045"); // excluded
    //filenames.push_back("run_0046"); // excluded
    //filenames.push_back("run_0047"); // excluded
    //filenames.push_back("run_0048"); // excluded
    //filenames.push_back("run_0049"); // excluded
    filenames.push_back("run_0050");
    filenames.push_back("run_0051");
    filenames.push_back("run_0052");
    //filenames.push_back("run_0053"); // excluded
    filenames.push_back("run_0054");
    //filenames.push_back("run_0055"); // excluded
    filenames.push_back("run_0056");
    //filenames.push_back("run_0057"); // excluded
    //filenames.push_back("run_0058"); // excluded
    //filenames.push_back("run_0059"); // excluded
    filenames.push_back("run_0060");
    filenames.push_back("run_0061");
    filenames.push_back("run_0062");
    //filenames.push_back("run_0063"); // excluded
    //filenames.push_back("run_0064"); // excluded
    //filenames.push_back("run_0065"); // excluded
    filenames.push_back("run_0066");
    //filenames.push_back("run_0067"); // excluded
    filenames.push_back("run_0068");
    filenames.push_back("run_0069");
    filenames.push_back("run_0070");
    filenames.push_back("run_0071");
    filenames.push_back("run_0072");
    filenames.push_back("run_0073");
    filenames.push_back("run_0074");
    filenames.push_back("run_0075");
    //filenames.push_back("run_0076"); // excluded
    //filenames.push_back("run_0077"); // excluded
    //filenames.push_back("run_0078"); // excluded
    //filenames.push_back("run_0079"); // excluded
    //filenames.push_back("run_0080"); // excluded
    //filenames.push_back("run_0081"); // excluded
    //filenames.push_back("run_0082"); // excluded
    filenames.push_back("run_0083");
    filenames.push_back("run_0084");
    //filenames.push_back("run_0085"); // excluded
    //filenames.push_back("run_0086"); // excluded
    //filenames.push_back("run_0087"); // excluded
    filenames.push_back("run_0088");
    filenames.push_back("run_0089");
    //filenames.push_back("run_0090"); // excluded
    //filenames.push_back("run_0091"); // excluded
    //filenames.push_back("run_0092"); // excluded
    //filenames.push_back("run_0093"); // excluded
    //filenames.push_back("run_0094"); // excluded
    //filenames.push_back("run_0095"); // excluded
    filenames.push_back("run_0096");
    filenames.push_back("run_0097");
    //filenames.push_back("run_0098"); // excluded
    filenames.push_back("run_0099");
    filenames.push_back("run_0100");
    filenames.push_back("run_0101");
    filenames.push_back("run_0102");
    filenames.push_back("run_0103");
    filenames.push_back("run_0104");
    filenames.push_back("run_0105");
    filenames.push_back("run_0106");
    filenames.push_back("run_0107");
    //filenames.push_back("run_0108"); // excluded
    //filenames.push_back("run_0109"); // excluded
    filenames.push_back("run_0110");
    //filenames.push_back("run_0111"); // excluded
    //filenames.push_back("run_0112"); // excluded
    filenames.push_back("run_0113");
    filenames.push_back("run_0114");
    filenames.push_back("run_0115");
    //filenames.push_back("run_0116"); //there is no scaler
    //filenames.push_back("run_0117"); // excluded
    //filenames.push_back("run_0118");
    filenames.push_back("run_0119");
    filenames.push_back("run_0120");
    filenames.push_back("run_0121");
    filenames.push_back("run_0122");
    filenames.push_back("run_0123");
    filenames.push_back("run_0124");
    //filenames.push_back("run_0125");
    //filenames.push_back("run_0126"); // excluded
    filenames.push_back("run_0127");
    filenames.push_back("run_0128");
    filenames.push_back("run_0129");
    filenames.push_back("run_0130");
    filenames.push_back("run_0131");
    //filenames.push_back("run_0132"); // excluded
    filenames.push_back("run_0133");
  
    /*  for(int i = 138; i <= 182; i++){
        if(i == 138 || i == 139 || i == 141 || i == 142 || i == 143 || i == 144 || i == 145 || i == 146 || i == 147 || i == 150 || i == 151 || i == 152 || i == 153 || i == 157 || i == 158 || i == 159 || i == 160 || i == 161 || i == 163 || i == 165 || i == 166 || i == 169 || i == 172 || i == 173 || i == 175 || i == 179 || i == 181)
            filenames.push_back("run_0" + std::to_string(i));
    }*/

/*
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
*/
    // zmin and zmax in meters for the cut on the vertex position.
    const Double_t zmin = 0.23;
    const Double_t zmax = 0.90; //0.63

    // Computation of the renormalization constant. Value of rho_H2 obtained from LISE++
    std::cout << "thickness = " << rho_D2 * (zmax - zmin) * 100 << " atoms/cm2.\n";

    // NI will count the total number of ions that entered the ATTPC. We iterate over the scaler files.
    Double_t NI {};
    for(auto filename: filenames){
        TFile *scalersFile = new TFile("rootfilesScalers/" + filename + "_scaler.root", "R");
        TTree *Tscalers = (TTree *) scalersFile->Get("parquettree");

        Long64_t Ic_ds {};
        Long64_t tF {};
        Long64_t tL {};
        Tscalers->SetBranchAddress("ic_ds", &Ic_ds); //downscale
        Tscalers->SetBranchAddress("trigger_free", &tF); 
        Tscalers->SetBranchAddress("trigger_live", &tL);

        for(int i = 0; i < Tscalers->GetEntries(); i++){
            Tscalers->GetEntry(i);
            if(tF > 0)
                NI = NI + 1000 * Ic_ds * tL / tF;
        }
    }
    std::cout << "N_I = " << NI << " total ions.\n";
    Double_t normalization_constant = 1 / (2 * TMath::Pi() * rho_D2 * 1E-27 * (zmax - zmin) * 100 * NI * bin_deg * TMath::Pi() / 180); // mbarn.
    std::cout << "normalization_constant = "<< normalization_constant << " mbarn.\n";

    // Compute the reconstruction efficiency for the gs from the attpc_engine simulation.
    /*for(auto filename: filenamesSim){
        TFile *simFile = new TFile("./data/v0.13.0/analysis2/simulation_15C_pp_gs_new/" + filename + "_" + std::to_string(A_ej) + "H.root", "R");
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
*/
    // Iterate over runs to load kinematic variables of interest.
    for(auto filename: filenames){
        TFile *physicsFile = new TFile("rootfiles/" + filename + "_" + std::to_string(A_ej) + "H.root", "R");
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

            auto [ex_energy, theta_cm] = kine_2b(m_C15, m_d, m_b, m_B, Ebeam_buff, theta_rad, E_ej);
            //auto [ex_energy, theta_cm] = kine_2b(m_C15, m_d, m_b, m_B, elossTableH2->GetEnergy(Ebeam_buff, zPos * 1000), theta_rad, E_ej);

            // Apply cuts.
            if(theta > 90 || theta < 10)
                continue;

            if(zPos < zmin || zPos > zmax)
              continue;

            //if(E_ej > 16)
                //continue;

            // Ground state cut.
            //if(ex_energy < -0.4 || ex_energy > 0.475)
                //continue;

            // Excitation energy vs Beam energy
            for (auto iEb = 0; iEb < 300; ++iEb) {
                auto [_ex_energy, _theta_cm] =  kine_2b(m_C15, m_d, m_b, m_B, iEb, theta_rad, E_ej);
                QvsEb->Fill(_ex_energy, iEb);
            }

            // Histograms
           //hredchi2->Fill(redchi);
            //if (ex_energy<-2.0 && ex_energy<2.0){
                Ang_Ener->Fill(theta, E_ej);
                if (cutg->IsInside(theta, E_ej)) {
                hex->Fill(ex_energy);
                }
            //}


            Double_t vx = TMath::Sin(theta_rad) * TMath::Sqrt(E_ej);
            Double_t vy = TMath::Cos(theta_rad) * TMath::Sqrt(E_ej);
            hVxVy->Fill(vx, vy);

            AngDistr->Fill(theta);
            AngDistrCM->Fill(theta_cm, normalization_constant / TMath::Sin(theta_cm * TMath::DegToRad()));
            //RatioToRutherford->Fill(theta_cm, normalization_constant / TMath::Sin(theta_cm * TMath::DegToRad()));

            ExvsZpos->Fill(ex_energy, zPos);
            thetavsZpos->Fill(theta, zPos);

           /* hredchi2vsE->Fill(redchi, E_ej);
            hredchi2vsEx->Fill(redchi, ex_energy);
            hredchi2vstheta->Fill(redchi, theta);*/

            hexvstheta->Fill(ex_energy, theta_cm);


        } // events
    } // Files

    //AngDistrCM->Divide(AngDistrCMEfficiency);

    auto ruther = new TF1("ruther","Rutherford(x, [0], [1], [2])", 1, 180);
    ruther->SetParameters(Ebeam_buff * m_d/(m_C15 + m_d), 1, 6);

    //RatioToRutherford->Divide(ruther);

    // Do fresco stuff
    //generateSFrescoSearchFile("test.search", "testTemplate.temp", AngDistrCM);
    /*performFitting("frescostuff/in_files/C15_pp.in", "frescostuff/in_files/C15_pp_in_template.sh",
                   "frescostuff/search_files/C15_pp.search", "frescostuff/search_files/C15_pp_search_template.temp",
                   "frescostuff/sfresco_template.sh", "frescostuff/sfresco.in", "frescostuff/C15_pp_sfresco.out",
                   Ebeam_buff, AngDistrCM);
*/
   // Kinematics
    TGraph *kine_gs = ReadKinematics("/home/georgina/fair_install/ATTPCROOTv2/macro/Kinematics/Decay_kinematics/Kine15C_gs.txt");
    TGraph *kine_1st = ReadKinematics("/home/georgina/fair_install/ATTPCROOTv2/macro/Kinematics/Decay_kinematics/Kine15C_1stState.txt");
   /*TF1 *fElasticPeak = new TF1("fElasticPeak", "gaus(0)", -1., 1.);
   Double_t paramsElastic[3] = {1000, 00 , 0.5};
   fElasticPeak->SetNpx(1000);
   fElasticPeak->SetParameters(paramsElastic);
   fElasticPeak->SetParLimits(1, -1., 1.);
   fElasticPeak->SetParLimits(2, 0., 2.);
   hex->Fit(fElasticPeak, "R");*/

    TCanvas *c_ExEnerBeam = new TCanvas();
    QvsEb->Draw("zcol");


   TCanvas *c1 = new TCanvas();
   c1->Divide(2, 1);
   c1->Draw();
   c1->cd(1);
   Ang_Ener->SetMarkerStyle(20);
   Ang_Ener->SetMarkerSize(0.5);
   Ang_Ener->Draw("col");
   Ang_Ener->GetXaxis()->SetTitle("Angle (deg)");
   Ang_Ener->GetYaxis()->SetTitle("Energy (MeV)");
   kine_gs->Draw("SAME");
   kine_1st->Draw("SAME");
   c1->cd(2);
   hVxVy->Draw("zcol");


   TCanvas *c_ExEner = new TCanvas();

   // Fits
   /*TF1 *fitFunc = new TF1("gs&1st", "gaus(0)+gaus(3)", -1.5, 1.5);
   fitFunc->SetParameters(1500, 0.1, 1, 1000, 0.6, 0.2);
   hex->Fit("gs&1st");*/
   /*TF1 *fitFunc = new TF1("gs&1st", "gaus(0)+gaus(3)+gaus(6)", -3, 1.5);
   fitFunc->SetParameters(1500, 0, 1, 1400, 0, 1, 1600, 0, 1);
   hex->Fit("gs&1st");*/
   // End Fits
   hex->Draw();
   hex->GetXaxis()->SetTitle("Excitation Energy (MeV)");
   hex->GetYaxis()->SetTitle("Counts");


   //auto [KK_XSection, VV_XSection, PP_XSection] = GetCalculatedXSections("calc_C15_pp_elastic.txt");


   TCanvas *c_AngDistr = new TCanvas();
   c_AngDistr->Divide(2,1);
   c_AngDistr->cd(1);
   AngDistr->Draw();
   c_AngDistr->cd(2);
   ruther->SetLineColor(1);
   //ruther->Draw("l");
   AngDistrCM->Draw();
   //KK_XSection->SetLineColor(2);
   //KK_XSection->Draw("same l");
   //VV_XSection->SetLineColor(3);
   //VV_XSection->Draw("same l");
   //PP_XSection->SetLineColor(4);
   //PP_XSection->Draw("same l");
   AngDistrCM->GetXaxis()->SetTitle("#theta_{CM} (deg)");
   AngDistrCM->GetYaxis()->SetTitle("#frac{d#sigma}{d#Omega} (#frac{mb}{sr})");

   /*TCanvas *c_RatioToRutherford = new TCanvas();
    RatioToRutherford->Draw();*/

   TCanvas *c_ExvsZpos = new TCanvas();
   ExvsZpos->Draw("zcol");
   ExvsZpos->GetXaxis()->SetTitle("Excitation Energy (MeV)");
   ExvsZpos->GetYaxis()->SetTitle("z (m)");


   /*TCanvas *c_redchi2 = new TCanvas();
   hredchi2->Draw("zcol");*/

   TCanvas *c_hex_vs_theta = new TCanvas();
   hexvstheta->Draw("zcol");
   hexvstheta->GetXaxis()->SetTitle("Energy (MeV)");
   hexvstheta->GetYaxis()->SetTitle("Angle (deg)");


    TCanvas *c_theta_vs_zPos = new TCanvas();
    thetavsZpos->Draw("zcol");
    thetavsZpos->GetXaxis()->SetTitle("theta_{LAB} (deg)");
    thetavsZpos->GetYaxis()->SetTitle("z (m)");


    /*TCanvas *c_redchi2_vs_E = new TCanvas();
    c_redchi2_vs_E->Divide(2,2);
    c_redchi2_vs_E->cd(1);
    hredchi2vsE->Draw("colz");
    c_redchi2_vs_E->cd(2);
    hredchi2vsEx->Draw("colz");
    c_redchi2_vs_E->cd(3);
    hredchi2vstheta->Draw("colz");*/

  /* TCanvas *cAngEnerSim = new TCanvas();
   Ang_EnerSim->SetMarkerStyle(20);
   Ang_EnerSim->SetMarkerSize(0.5);
   Ang_EnerSim->Draw("col");
   Ang_EnerSim->GetXaxis()->SetTitle("Angle (deg)");
   Ang_EnerSim->GetYaxis()->SetTitle("Energy (MeV)");
   kine_gs->Draw("SAME");
   kine_1st->Draw("SAME");

   TCanvas *cExEnerSim = new TCanvas();
   hexSim->Draw();
   hexSim->GetXaxis()->SetTitle("Excitation Energy (MeV)");
   hexSim->GetYaxis()->SetTitle("Counts");*/

   TCanvas *cAngDistrEfficiency = new TCanvas();
   AngDistrCMEfficiency->Draw();
   AngDistrCMEfficiency->GetXaxis()->SetTitle("#theta_{CM} (deg)");
   AngDistrCMEfficiency->GetYaxis()->SetTitle("#epsilon");
}

void GetEnergy(Double_t M, Double_t IZ, Double_t BRO, Double_t &E)
{

   // Energy per nucleon
   Float_t AM = 931.5;
   Float_t X = BRO / 0.1439 * IZ / M;
   X = pow(X, 2);
   X = 2. * AM * X;
   X = X + pow(AM, 2);
   E = TMath::Sqrt(X) - M;
}

std::tuple<TGraph*, TGraph*, TGraph*> GetCalculatedXSections(std::string calcFile){
    Double_t *theta = new Double_t[900];
    Double_t *KK = new Double_t[900];
    Double_t *VV = new Double_t[900];
    Double_t *PP = new Double_t[900];

    std::ifstream *calcStr = new std::ifstream(calcFile);
    int i = 0;

    while(!calcStr->eof()){
        *calcStr >> theta[i] >> KK[i] >> VV[i] >> PP[i];
        i++;
    }
    TGraph *KK_XSection = new TGraph(i, theta, KK);
    TGraph *VV_XSection = new TGraph(i, theta, VV);
    TGraph *PP_XSection = new TGraph(i, theta, PP);
    return std::make_tuple(KK_XSection, VV_XSection, PP_XSection);
}

TGraph* ReadKinematics(std::string kineFile){
    Double_t *ThetaCMS = new Double_t[20000];
    Double_t *ThetaLabRec = new Double_t[20000];
    Double_t *EnerLabRec = new Double_t[20000];
    Double_t *ThetaLabSca = new Double_t[20000];
    Double_t *EnerLabSca = new Double_t[20000];
    Double_t *MomLabRec = new Double_t[20000];

    TString fileKine = kineFile;
    std::ifstream *kineStr = new std::ifstream(fileKine.Data());
    Int_t numKin = 0;

    if (!kineStr->fail()){
        while (!kineStr->eof()){
            *kineStr >> ThetaCMS[numKin] >> ThetaLabRec[numKin] >> EnerLabRec[numKin] >>
                        ThetaLabSca[numKin] >> EnerLabSca[numKin];
            // numKin++;

            // MomLabRec[numKin] =( pow(EnerLabRec[numKin] + M_Ener,2) - TMath::Power(M_Ener, 2))/1000.0;
            // std::cout<<" Momentum : " <<MomLabRec[numKin]<<"\n";
            // Double_t E = TMath::Sqrt(TMath::Power(p, 2) + TMath::Power(M_Ener, 2)) - M_Ener;
            numKin++;
       }
    } else if (kineStr->fail())
        std::cout << " Warning : No Kinematics file found for this reaction!" << std::endl;

    TGraph *kine = new TGraph(numKin, ThetaLabRec, EnerLabRec);
    return kine;
}
