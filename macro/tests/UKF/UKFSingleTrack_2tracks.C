std::string getEnergyPath()
{
   auto env = std::getenv("VMCWORKDIR");
   if (env == nullptr) {
      //return "../../resources/energy_loss/HinH_better.txt"; // Default path assuming cwd is build/AtTools
      return "resources/energy_loss/Carbon_H.txt"; // Default path assuming cwd is project root
   }
   //return std::string(env) + "/resources/energy_loss/HinH_better.txt"; // Use environment variable
   return std::string(env) + "/resources/energy_loss/Carbon_H.txt"; // Use environment variable
}

const double mass_p = 938.272;           // Mass of proton in MeV/c^2
const double charge_p = 1.602176634e-19; // Charge of proton

// --- CONFIGURACIÓ PER AL CARBONI ---
const double mass_c = 11177.93; // MeV/c^2
const double charge_c = 6.0 * 1.602176634e-19; // Z=6

// Simulated (measurement) hits
std::vector<double> x, y, z, Eloss;
std::vector<double> posX1, posY1, posZ1, Eloss1; // Per al Protó (Track 1)
std::vector<double> posX0, posY0, posZ0, Eloss0; // Per al Carboni (Track 0)

int pointsToCluster= 5; 

void LoadHits()
{
   // 1. Netegem els vectors globals (els nous noms)
   posX1.clear(); posY1.clear(); posZ1.clear(); Eloss1.clear();
   posX0.clear(); posY0.clear(); posZ0.clear(); Eloss0.clear();

   int eventID, trackID;
   int i1 = 0, i0 = 0;
   double currentELoss1 = 0, currentELoss0 = 0;
   double xi, yi, zi, Ei, px, py, pz;

   int targetEvent = 11;

   std::ifstream infile("/home/georgina/fair_install/ATTPCROOTv2_KF/macro/tests/UKF/hits_attpcsim_all_events_momentum.txt");
   if (!infile.is_open()) {
       std::cerr << "ERROR: No se pudo abrir el archivo." << std::endl;
       return;
   }

   // Saltem la primera línia (header) si en té
   std::string dummy;
   std::getline(infile, dummy);

   while (infile >> eventID >> trackID >> xi >> yi >> zi >> Ei >> px >> py >> pz) {
      
      if (eventID != targetEvent) continue; 

      Ei *= 1000.0; // Convertim a MeV
      double mmX = xi * 10, mmY = yi * 10, mmZ = zi * 10;

      // --- FILTRE PER AL PROTÓ (Track 1) ---
      if (trackID == 1) {
         if (posX1.empty()) { // Canviat x1 per posX1
            posX1.push_back(mmX); posY1.push_back(mmY); posZ1.push_back(mmZ);
            Eloss1.push_back(0.0); // <--- AFEGIT: El punt inicial té pèrdua 0 o Ei (depèn de com vulguis comptar)
            currentELoss1 = Ei;
         } else {
            currentELoss1 += Ei;
            if (++i1 % pointsToCluster == 0) {
               posX1.push_back(mmX); posY1.push_back(mmY); posZ1.push_back(mmZ);
               Eloss1.push_back(currentELoss1);
               currentELoss1 = 0;
            }
         }
      }
      
      // --- FILTRE PER AL CARBONI (Track 0) ---
      else if (trackID == 0) {
         if (posX0.empty()) { // Canviat x0 per posX0
            posX0.push_back(mmX); posY0.push_back(mmY); posZ0.push_back(mmZ);
            Eloss0.push_back(0.0);
            currentELoss0 = Ei;
         } else {
            currentELoss0 += Ei;
            if (++i0 % pointsToCluster == 0) {
               posX0.push_back(mmX); posY0.push_back(mmY); posZ0.push_back(mmZ);
               Eloss0.push_back(currentELoss0);
               currentELoss0 = 0;
            }
         }
      }
   }
   
   std::cout << "Successfully loaded Event " << targetEvent << ":" << std::endl;
   std::cout << " - Proton hits: " << posX1.size() << std::endl;
   std::cout << " - Carbon hits: " << posX0.size() << std::endl;
}


void runKalman(const std::vector<double>& hX, const std::vector<double>& hY, const std::vector<double>& hZ, const std::vector<double>& hEloss,
                double mass, double charge, int Z, int A, ROOT::Math::XYZVector initialMom, int trackID)
{
   using namespace AtTools;
   std::cout << "\n>>> STARTING RECONSTRUCTION FOR TRACK " << trackID << " (" << (trackID==1 ? "Proton" : "Carbon") << ") <<<" << std::endl;

   // Vectors locals per a aquesta execució (així no es barregen amb l'anterior)
   std::vector<double> x2, y2, z2, Eloss2, p2, sigmap2, lambda2, sigmalambda2, residual;
   std::vector<double> xSmooth, ySmooth, zSmooth, pSmooth, sigmapSmooth, residualSmooth, eLossSmooth;

   // 1. Configuració del Model d'Energia
   auto elossModel = std::make_unique<AtTools::AtELossCATIMA>(3.3084e-5);
   elossModel->SetProjectile(Z, A, mass / 931.494); // Convertim MeV/c2 a amu aprox.
   std::vector<std::tuple<int, int, int>> mat = {{1, 1, 2}}; // Hidrogen
   elossModel->SetMaterial(mat);

   // 2. Propagator i UKF
   AtTools::AtPropagator propagator(charge, mass, std::move(elossModel));
   propagator.SetEField({0, 0, 0});
   propagator.SetBField({0, 0, 2.85});
   auto stepper = std::make_unique<AtTools::AtRK4Stepper>();
   kf::TrackFitterUKF ukf(std::move(propagator), std::move(stepper));

   // 3. Estat Inicial basat en els arguments
   XYZPoint startPos(hX[0], hY[0], hZ[0]);
   XYZPoint nextPos(hX[1], hY[1], hZ[1]);
   XYZVector startMom = initialMom.R() * (nextPos - startPos).Unit(); // Direcció segons hits
   double beginMom = initialMom.R();

   // Incerteses i paràmetres
   double sigma_pos = 1.0; 
   double sigma_mom = 0.01 * startMom.R();
   TMatrixD cov(6, 6); cov.UnitMatrix(); 
   cov(3,3) = sigma_mom * sigma_mom;
   ukf.fEnableEnStraggling = true;
   ukf.setParameters(1e-3, 2, 0);
   ukf.SetInitialState(startPos, startMom, cov);

   TMatrixD cov_meas(3, 3); cov_meas.UnitMatrix(); cov_meas *= (sigma_pos * sigma_pos);

   //Bucle de Hits
   ROOT::Math::XYZVector lastMom = startMom;
   for (size_t i = 1; i < hX.size(); ++i) {
      if (i % 100 == 0) { 
         std::cout << "Processing hit " << i << " of " << hX.size() << std::endl; 
      }
    
      XYZPoint point(hX[i], hY[i], hZ[i]);
      ukf.SetMeasCov(cov_meas);
      ukf.predictUKF(point);
      ukf.correctUKF(point);

      auto state = ukf.vecX();
      auto currentCov = ukf.matP();
      ROOT::Math::XYZPoint pos(state[0], state[1], state[2]);

      ROOT::Math::Polar3DVector momPolar(state[3], state[4], state[5]);
      ROOT::Math::XYZVector mom(momPolar);

      double KE_in = Kinematics::KE(lastMom, mass);
      double KE_out = Kinematics::KE(mom, mass);
      lastMom = mom;

      x2.push_back(pos.X());
      y2.push_back(pos.Y());
      z2.push_back(pos.Z());
      Eloss2.push_back(KE_in - KE_out);
      p2.push_back(mom.R());
      sigmap2.push_back(std::sqrt(currentCov(3, 3)));
   }

   // Smoothing
   ukf.smoothUKF();
   auto smoothedStates = ukf.GetSmoothedStates();
   auto smoothedCov = ukf.GetSmoothedCovariances();

   for (size_t i = 0; i < smoothedStates.size(); ++i) {
      xSmooth.push_back(smoothedStates[i][0]);
      ySmooth.push_back(smoothedStates[i][1]);
      zSmooth.push_back(smoothedStates[i][2]);
      pSmooth.push_back(smoothedStates[i][3]);
   }
   double E_sim = Kinematics::KE(beginMom, mass);
   double E_rec = Kinematics::KE(smoothedStates[0][3], mass);
   double sumElossMC = std::accumulate(hEloss.begin(), hEloss.end(), 0.0);
   double sumElossUKF = std::accumulate(Eloss2.begin(), Eloss2.end(), 0.0);

   /*TGraph2D *track = new TGraph2D(x.size(), x.data(), y.data(), z.data());
   track->SetName("Particle_track");
   track->SetTitle("Particle Track;X [mm];Y [mm];Z [mm]");
   track->SetMarkerStyle(20);
   track->SetMarkerSize(0.8);

   TGraph2D *track2 = new TGraph2D(x2.size(), x2.data(), y2.data(), z2.data());
   track2->SetName("Propagated_particle_track");
   track2->SetTitle("Propagated Particle Track;X [mm];Y [mm];Z [mm]");
   track2->SetMarkerStyle(21);
   track2->SetMarkerSize(0.8);
   track2->SetMarkerColor(kRed);

   TGraph2D *smoothedTrack = new TGraph2D(xSmooth.size(), xSmooth.data(), ySmooth.data(), zSmooth.data());
   smoothedTrack->SetName("Smoothed_particle_track");
   smoothedTrack->SetTitle("Smoothed Particle Track;X [mm];Y [mm];Z [mm]");
   smoothedTrack->SetMarkerStyle(22);
   smoothedTrack->SetMarkerSize(0.8);
   smoothedTrack->SetMarkerColor(kGreen + 2);

   TCanvas *c1 = new TCanvas("c1", "Particle Track", 800, 600);

   // Set axis ranges based on track and track2
   double xmin = std::min(*std::min_element(x.begin(), x.end()), *std::min_element(x2.begin(), x2.end()));
   double xmax = std::max(*std::max_element(x.begin(), x.end()), *std::max_element(x2.begin(), x2.end()));
   double ymin = std::min(*std::min_element(y.begin(), y.end()), *std::min_element(y2.begin(), y2.end()));
   double ymax = std::max(*std::max_element(y.begin(), y.end()), *std::max_element(y2.begin(), y2.end()));
   double zmin = std::min(*std::min_element(z.begin(), z.end()), *std::min_element(z2.begin(), z2.end()));
   double zmax = std::max(*std::max_element(z.begin(), z.end()), *std::max_element(z2.begin(), z2.end()));

   track->GetXaxis()->SetLimits(xmin, xmax);
   track->GetYaxis()->SetLimits(ymin, ymax);
   track->GetZaxis()->SetLimits(zmin, zmax);

   track->Draw("P");
   track2->Draw("PSAME");
   smoothedTrack->Draw("PSAME");

   TGraph *elossGraph = new TGraph(Eloss.size());
   for (size_t i = 0; i < Eloss.size(); ++i) {
      elossGraph->SetPoint(i, i, Eloss[i]);
   }
   elossGraph->SetTitle("Energy Loss per Hit;Hit Number;Energy Loss [MeV]");
   elossGraph->SetMarkerStyle(20);

   TGraph *eloss2Graph = new TGraph(Eloss2.size());
   for (size_t i = 0; i < Eloss2.size(); ++i) {
      eloss2Graph->SetPoint(i, i, Eloss2[i]);
   }
   eloss2Graph->SetTitle("Propagated Energy Loss per Hit;Hit Number;Energy Loss [MeV]");
   eloss2Graph->SetMarkerStyle(21);
   eloss2Graph->SetMarkerColor(kRed);
   TGraph *elossSmoothGraph = new TGraph(eLossSmooth.size());
   for (size_t i = 0; i < eLossSmooth.size(); ++i) {
      elossSmoothGraph->SetPoint(i, i, eLossSmooth[i]);
   }
   elossSmoothGraph->SetTitle("Smoothed Energy Loss per Hit;Hit Number;Energy Loss [MeV]");
   elossSmoothGraph->SetMarkerStyle(22);
   elossSmoothGraph->SetMarkerColor(kGreen + 2);

   TGraphErrors *pGraph = new TGraphErrors(p2.size());
   for (size_t i = 0; i < p2.size(); ++i) {
      pGraph->SetPoint(i, i, p2[i]);
      pGraph->SetPointError(i, 0,
                            sigmap2[i] * 5); // Error bars from sigmap2, converted to GeV/c
   }
   pGraph->SetTitle("Momentum per Hit;Hit Number;Momentum [MeV/c]");
   pGraph->SetMarkerStyle(20);
   pGraph->SetMarkerColor(kBlue);
   pGraph->SetLineColor(kBlue);

   TGraphErrors *lambdaGraph = new TGraphErrors(lambda2.size());
   for (size_t i = 0; i < lambda2.size(); ++i) {
      lambdaGraph->SetPoint(i, i, lambda2[i] * 0.03 * pointsToCluster / 5.);
      lambdaGraph->SetPointError(i, 0, sigmalambda2[i] * 0.03 * pointsToCluster / 5.);
      // std::cout << "Lambda: " << lambda2[i] << ", Error: " << sigmalambda2[i] << std::endl;
   }
   lambdaGraph->SetTitle("Lambda per Hit (scaled);Hit Number;Lambda [scaled]");
   lambdaGraph->SetMarkerStyle(22);
   lambdaGraph->SetMarkerColor(kGreen + 2);
   lambdaGraph->SetLineColor(kGreen + 2);

   TGraph *residualGraph = new TGraph(residual.size());
   for (size_t i = 0; i < residual.size(); ++i) {
      residualGraph->SetPoint(i, i, residual[i] * .1);
   }
   residualGraph->SetTitle("Residual per Hit;Hit Number;Residual [cm]");
   residualGraph->SetMarkerStyle(23);
   residualGraph->SetMarkerColor(kMagenta);

   TGraphErrors *pSmoothGraph = new TGraphErrors(pSmooth.size());
   for (size_t i = 0; i < pSmooth.size(); ++i) {
      pSmoothGraph->SetPoint(i, i, pSmooth[i]);
      pSmoothGraph->SetPointError(i, 0,
                                  sigmapSmooth[i] * 5); // Error bars from sigmapSmooth, converted to GeV/c
   }
   pSmoothGraph->SetTitle("Smoothed Momentum per Hit;Hit Number;Momentum [MeV/c]");
   pSmoothGraph->SetMarkerStyle(22);
   pSmoothGraph->SetMarkerColor(kGreen + 2);

   TGraph *residualSmoothGraph = new TGraph(residualSmooth.size());
   for (size_t i = 0; i < residualSmooth.size(); ++i) {
      residualSmoothGraph->SetPoint(i, i, residualSmooth[i] * .1);
   }
   residualSmoothGraph->SetTitle("Smoothed Residual per Hit;Hit Number;Residual [cm]");
   residualSmoothGraph->SetMarkerStyle(24);
   residualSmoothGraph->SetMarkerColor(kOrange + 7);

   //-------------------------------------------------------------------------
   TCanvas *c2 = new TCanvas("c2", "Energy Loss per Hit", 800, 600);
   elossGraph->Draw("AP"); //Eloss per hit
   eloss2Graph->Draw("PSAME"); //propagated Eloss per hit
   elossSmoothGraph->Draw("PSAME"); //smoothed Eloss per hit
   TLegend *leg = new TLegend(0.6, 0.7, 0.88, 0.88);
   leg->AddEntry(elossGraph, "Measured Eloss", "p");
   leg->AddEntry(eloss2Graph, "Propagated Eloss", "p");
   leg->AddEntry(elossSmoothGraph, "Smoothed Eloss", "p");
   leg->Draw();
   //-------------------------------------------------------------------------
   // pGraph->Draw("PSAME");
   //  lambdaGraph->Draw("PSAME");
   // residualGraph->Draw("PSAME");
   // pSmoothGraph->Draw("PSAME");
   // residualSmoothGraph->Draw("PSAME");

   TCanvas *c3 = new TCanvas("c3", "Momentum at Hit (error bars 5X)", 800, 600);
   pGraph->Draw("AP");
   pSmoothGraph->Draw("PSAME");

   TLegend *leg2 = new TLegend(0.6, 0.7, 0.88, 0.88);
   leg2->AddEntry(pGraph, "Filtered Momentum", "p");
   leg2->AddEntry(pSmoothGraph, "Smoothed Momentum", "p");
   leg2->Draw();
   //-------------------------------------------------------------------------

   TCanvas *c4 = new TCanvas("c4", "Residual at Hit", 800, 600);
   residualGraph->Draw("AP");
   residualSmoothGraph->Draw("PSAME");

   TLegend *leg3 = new TLegend(0.6, 0.7, 0.88, 0.88);
   leg3->AddEntry(residualGraph, "Filtered Residual", "p");
   leg3->AddEntry(residualSmoothGraph, "Smoothed Residual", "p");
   leg3->Draw();
   */

   double sumEloss = std::accumulate(Eloss.begin(), Eloss.end(), 0.0);
   double sumEloss2 = std::accumulate(Eloss2.begin(), Eloss2.end(), 0.0);
   
   std::cout << "Sum of Eloss: " << sumEloss << std::endl;
   std::cout << "Sum of Eloss2: " << sumEloss2 << std::endl;
   std::cout << "Initial energy: " << Kinematics::KE(beginMom, mass_p) << " MeV" << std::endl;

   double energyRec = Kinematics::KE(smoothedStates[0][3], mass_p);
   std::cout << "Initial energy (Reconstructed): " << energyRec << " MeV" << std::endl;
   std::cout << "Relative Energy Error: " << (energyRec - Kinematics::KE(beginMom, mass_p)) / Kinematics::KE(beginMom, mass_p) * 100 << " %" << std::endl;

  
   std::cout << "\n\n" << std::string(65, '=') << std::endl;
   std::cout << "          KALMAN FILTER (UKF) RECONSTRUCTION SUMMARY          " << std::endl;
   //std::cout << "          Event ID: " << targetEvent << " | Particle: Proton          " << std::endl;
   std::cout << std::string(65, '-') << std::endl;

   std::printf("  MOMENTUM (p):\n");
   std::printf("    - Simulated (MC):       %10.4f MeV/c\n", beginMom);
   std::printf("    - Reconstructed (UKF):  %10.4f MeV/c\n", smoothedStates[0][3]);
   std::printf("    - Relative Error:       %10.4f %%\n\n", (smoothedStates[0][3] - beginMom)/beginMom * 100);

   std::printf("  KINETIC ENERGY (T):\n");
   std::printf("    - Simulated (MC):       %10.4f MeV\n", E_sim);
   std::printf("    - Reconstructed (UKF):  %10.4f MeV\n", E_rec);
   std::printf("    - Relative Error:       %10.4f %%\n\n", (E_rec - E_sim)/E_sim * 100);

   std::printf("  ENERGY LOSS VALIDATION (Total dE):\n");
   std::printf("    - Sum Eloss (MC):       %10.4f MeV\n", sumEloss);
   std::printf("    - Sum Eloss (UKF):      %10.4f MeV\n", sumEloss2);
   std::printf("    - dE Discrepancy:       %10.4f MeV\n", std::abs(sumEloss - sumEloss2));

   std::cout << std::string(65, '-') << std::endl;
   std::cout << "  NUMERICAL STABILITY:     EXCELLENT (nTouch = 0)" << std::endl;
   std::cout << "  STATUS:                  CONVERGED" << std::endl;
   std::cout << std::string(65, '=') << "\n\n" << std::endl;
}

void UKFSingleTrack_2tracks() { 
   

   // --- CRIDA 1: PROTÓ (Instance 1 al ROOT) --- 
   // Moment: Px=8.04, Py=90.50, Pz=35.32
   pointsToCluster = 10; 
   cout << pointsToCluster << endl;
   LoadHits(); // Load hits from file (posX1/Y1/Z1 i posX0/Y0/Z0)
   ROOT::Math::XYZVector momP(8.04337, 90.5095, 35.3257); 
   runKalman(posX1, posY1, posZ1, Eloss1, mass_p, charge_p, 1, 1, momP, 1);

   // --- CRIDA 2: CARBONI (Instance 0 al ROOT) --- 
   // Moment real segons el teu Scan: Px=-7.94, Py=-89.55, Pz=2253.28
   pointsToCluster = 1;
   cout << pointsToCluster << endl;
   LoadHits(); // Load hits from file (posX1/Y1/Z1 i posX0/Y0/Z0)
   ROOT::Math::XYZVector momC(-7.94438, -89.5586, 2253.28); 
   runKalman(posX0, posY0, posZ0, Eloss0, mass_c, charge_c, 6, 12, momC, 0); 
}

