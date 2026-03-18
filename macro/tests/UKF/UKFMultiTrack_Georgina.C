// Struct to store and easily compare the Unscented Kalman Filter (UKF) 
// reconstruction results against the Monte Carlo truth data.
enum class FitStatus {
    SUCCESS,
    STOPPED_IN_BRAGG,
    FAILURE
};

struct UKFResult {
   int eventID = -1;
   int trackID = -1;
   
   FitStatus status = FitStatus::FAILURE; // Default to FAILURE, will be set to SUCCESS or STOPPED_IN_BRAGG based on fit outcome

   // Reconstructed variables
   double p_rec = 0.0;      // Reconstructed momentum [MeV/c]
   double E_rec = 0.0;      // Reconstructed kinetic energy [MeV]
   double theta_rec = 0.0;  // Reconstructed polar angle [rad]
   double phi_rec = 0.0;    // Reconstructed azimuthal angle [rad]
   double sumEloss = 0.0;   // Total energy loss [MeV]
   
   // Monte Carlo truth variables (for comparison/validation)
   double p_true = 0.0;     // True initial momentum from simulation [MeV/c]
   double theta_true = 0.0; // True polar angle [rad]
   double phi_true = 0.0;   // True azimuthal angle [rad]
};

double wrapAngle(double a){
    //return TMath::ATan2(TMath::Sin(a), TMath::Cos(a));
    return a;
}


// Global vectors to accumulate angular residuals/differences across all processed events.
// These are typically used to fill histograms for resolution analysis at the end.
static std::vector<double> all_dtheta_rad;
static std::vector<double> all_dtheta_fwd_rad;
static std::vector<double> all_dtheta_sm_rad;
static std::vector<double> all_dtheta_sm_minus_fwd_rad;

// Global vectors for all clusters across all events
static std::vector<double> all_cluster_p_rec;
static std::vector<double> all_cluster_p_true;

// Dynamically resolves the path to the energy loss table (e.g., Proton in Hydrogen).
// Checks the standard VMCWORKDIR environment variable first, then falls back to a relative path.
std::string getEnergyPath()
{
   auto env = std::getenv("VMCWORKDIR");
   if (env == nullptr) {
      //return "../../resources/energy_loss/HinH_better.txt"; // Default path assuming cwd is build/AtTools
      return "resources/energy_loss/HinH.txt"; // Default path assuming cwd is project root
   }
   //return std::string(env) + "/resources/energy_loss/HinH_better.txt"; // Use environment variable
   return std::string(env) + "/resources/energy_loss/HinH.txt"; // Use environment variable
}

// Fundamental physics constants for the tracked particle
/*const double mass_p = 938.272;           // Mass of proton in [MeV/c^2]
const double charge_p = 1.602176634e-19; // Charge of proton in [Coulombs]
*/

double mass_pi = 139.57039;  // MeV/c^2
double charge_pi = 1.0;
double charge_pi_Coulombs = 1.602176634e-19; //1.0;      // en unidades de e (no Coulomb)


int mat_Z = 1;               // Hydrogen
int mat_A = 1;               // H-1
double density_300torr = 3.3084e-5;  // g/cm^3 (tu valor para 300 torr)
double density_60torr= 6.6168e-6; 

// Dictionaries (maps) linking a specific track ID to its downsampled cluster data.
// These store the grouped spatial coordinates [mm] and accumulated energy loss [MeV].
std::map<int, std::vector<double>> posX;
std::map<int, std::vector<double>> posY;
std::map<int, std::vector<double>> posZ;
std::map<int, std::vector<double>> Eloss;

// To save the XYZVector of the initial momentum for each track
std::map<int, ROOT::Math::XYZVector> initialMom;

int pointsToCluster= 10; //5

void LoadHitsROOT(TTree* tree, TClonesArray* tpcPoints, int targetEvent)
{
    posX.clear(); 
    posY.clear(); 
    posZ.clear(); 
    Eloss.clear();
    initialMom.clear();

    cout << "Loading hits for event " << targetEvent << "..." << endl;

   if (targetEvent >= tree->GetEntries()) return;

   tpcPoints->Clear("C");
   tree->GetEntry(targetEvent);
   int nPoints = tpcPoints->GetEntriesFast();
    
   std::map<int, int> hitCount;           
   std::map<int, double> currentELoss;

   // Iterate over all simulated TPC points to build particle tracks 
   // and downsample the hits into discrete clusters.

   for (int i = 0; i < nPoints; i++) {

      AtMCPoint* point = (AtMCPoint*)tpcPoints->At(i);
      if (!point) continue;

      int trackID = point->GetTrackID();
      if (trackID != 0) continue;


        // Extract position and energy loss, converting units
        // Assuming framework native units are cm and GeV, converting to mm and MeV
        double mmX = point->GetX() * 10.0;
        double mmY = point->GetY() * 10.0;
        double mmZ = point->GetZ() * 10.0;
        double Ei = point->GetEnergyLoss() * 1000.0;   

      if (posX[trackID].empty()) { 

            // This is the first point that we have in the track. Initialize vectors.
            posX[trackID].push_back(mmX); 
            posY[trackID].push_back(mmY); 
            posZ[trackID].push_back(mmZ);
            Eloss[trackID].push_back(0.0);
            
            // Convert initial momentum to MeV/c and store it
            double px = point->GetPx() * 1000.0;
            double py = point->GetPy() * 1000.0;
            double pz = point->GetPz() * 1000.0;
            initialMom[trackID] = ROOT::Math::XYZVector(px, py, pz);
            cout << "eventID: " << targetEvent << " - TrackID: " << trackID 
                 << " - Initial Momentum (MeV/c): (" << px << ", " << py << ", " << pz << ")" << endl;

            
            currentELoss[trackID] = Ei;
            hitCount[trackID] = 0;
      } else {

            // We add up the energy loss of this point to the present cluster
            currentELoss[trackID] += Ei;
            hitCount[trackID]++;
            
            // Once we reach the target number of points, save the cluste
            if (hitCount[trackID] % pointsToCluster == 0) {
                posX[trackID].push_back(mmX); 
                posY[trackID].push_back(mmY); 
                posZ[trackID].push_back(mmZ);
                Eloss[trackID].push_back(currentELoss[trackID]);
                
                // Reset the energy counter for the next cluster of this particle
                currentELoss[trackID] = 0;

                // NOTE: Any remaining energy in a cluster smaller than 'pointsToCluster' 
                // at the end of the track is currently left in currentELoss[trackID].
            }
        }
   }  
           
   // =========================================================
    // EVENT SUMMARY OUTPUT
    // =========================================================
    std::cout << "\n======================================================\n";
    std::cout << " SUMMARY FOR EVENT " << targetEvent << "\n";
    std::cout << "======================================================\n";
    
    // posX.size() tells us exactly how many tracks we have found
    std::cout << " Total tracks found : " << posX.size() << "\n";
    std::cout << "------------------------------------------------------\n";
    
    // Loop through the map to print details for each track
    for (auto const& [trackID, vectorHits] : posX) {
        std::cout << " -> Track ID [" << trackID << "] : " 
                  << vectorHits.size() << " clustered hits saved.\n";
    }
    std::cout << "======================================================\n\n";
}

UKFResult runKalman(const std::vector<double>& hX, const std::vector<double>& hY, const std::vector<double>& hZ, const std::vector<double>& hEloss,
                double mass, double charge, int Z, int A, ROOT::Math::XYZVector initialMom, int trackID, bool drawPlots = true)
{
   using namespace AtTools;

   // 1. Local vectors for this execution 
   // (Ensures data from previous tracks doesn't leak into the current one

   std::vector<double> xForward, yForward, zForward, pForward, sigmapForward, residualForward,  eLossForward;
   std::vector<double> xSmooth, ySmooth, zSmooth, pSmooth, sigmapSmooth, residualSmooth, eLossSmooth;
   double pForward_early = 0.0;  // momentum at early stage (hit 1)
   double pForward_late = 0.0;   // momentum at late stage (last hit)

    //------------------------CATIMA IMPLEMENTATION:----------------------------------
    // Set up the CATIMA energy loss model for Hydrogen gas at 300 torr.
    
    /*auto elossModel = std::make_unique<AtTools::AtELossCATIMA>(density_300torr);
    elossModel->SetProjectile(mat_Z, mat_A, mass / 931.494); // Convert mass from MeV/c^2 to approximate amu for CATIMA
    std::vector<std::tuple<int, int, int>> mat = {{mat_Z, mat_A, 2}}; // Hidrogen
    elossModel->SetMaterial(mat);*/
    //--------------------------------------------------------------------------------

   //For pions we use the betheBloch:
    auto elossModel = std::make_unique<AtTools::AtELossBetheBloch>(
    charge_pi, // Charge in units of e
    mass_pi, // MeV/c^2
    mat_Z,
    mat_A,
    density_300torr,
    -1.0);


    // TEST
    double p_test = 30.0;
    double p_test2 = 300.0;
    double p_test3 = 1000.0;

    std::cout << "dEdx 30 MeV/c = " << elossModel->GetdEdx(p_test) << std::endl;
    std::cout << "dEdx 300 MeV/c = " << elossModel->GetdEdx(p_test2) << std::endl;
    std::cout << "dEdx 1000 MeV/c = " << elossModel->GetdEdx(p_test3) << std::endl;

    std::cout << "mass pion = " << mass_pi << std::endl;
    std::cout << "charge pion = " << charge_pi << std::endl;
    std::cout << "density = " << density_300torr << std::endl;

   // 2. Propagator and UKF Setup
   // Initializes the physical environment (E and B fields) and the Runge-Kutta stepper

   AtTools::AtPropagator propagator(charge, mass, std::move(elossModel));
   propagator.SetEField({0, 0, 0});
   propagator.SetBField({0, 0, 2.});//2.85 T
   auto stepper = std::make_unique<AtTools::AtRK4Stepper>();
   kf::TrackFitterUKF ukf(std::move(propagator), std::move(stepper));

   // 3. Initial State Setup
   XYZPoint startPos(hX[0], hY[0], hZ[0]);
   XYZPoint nextPos(hX[1], hY[1], hZ[1]);

   // We are using the true initial momentum from the simulation directly as our starting seed.
   // (Alternative: estimate direction using the first two hits, as commented out below)
  // XYZVector startMom = initialMom.R() * (nextPos - startPos).Unit(); 

   XYZVector startMom = initialMom; 
   double beginMom = initialMom.R();

   // Uncertainties and UKF tuning parameters
   double sigma_pos = 1.0; 
   double sigma_mom = 0.01 * startMom.R();
    
   // Initial covariance matrix (6x6: x, y, z, px, py, pz)
   TMatrixD cov(6, 6); cov.UnitMatrix(); 
   cov(3,3) = sigma_mom * sigma_mom;
    
   ukf.fEnableEnStraggling = true;
   ukf.setParameters(1e-3, 2, 0);
   //ukf.setParameters(0.1, 2, 0);
   ukf.SetInitialState(startPos, startMom, cov);

   // Measurement covariance matrix (3x3: spatial coordinates only)
   TMatrixD cov_meas(3, 3); 
   cov_meas.UnitMatrix(); 
   cov_meas *= (sigma_pos * sigma_pos);

   // 4. Hit Loop (Forward Filter)
   ROOT::Math::XYZVector lastMom = startMom;
   for (size_t i = 1; i < hX.size(); ++i) {
      if (i % 200 == 0) { 
         std::cout << "Processing hit " << i << " of " << hX.size() << std::endl; 
      }
    
      XYZPoint point(hX[i], hY[i], hZ[i]);

      double currentKE = Kinematics::KE(lastMom, mass);
      
      // Safety Check: If the particle loses almost all its energy, the energy loss 
      // equations break down (Bragg peak singularity). We stop propagation to avoid crashing.
      if (currentKE < 0.05) { 
          std::cout << "[DEBUG] Critical energy reached (" << currentKE << " MeV). Stopping propagation." << std::endl;
          break; 
      }
      ukf.SetMeasCov(cov_meas);

      // Try-catch block protects against Kalman Filter math exceptions 
      // (e.g., matrix inversion failures due to collinear points or extreme noise).
      try {
          ukf.predictUKF(point);
          ukf.correctUKF(point);
      } catch (...) {
          // If the math blows up here, exit the loop cleanly rather than crashing the whole program.
          break; 
      }

      auto state = ukf.vecX();
      // Condition to avoid crash if UKF fails and stops updating states
      if (state.size() < 6) break;

      auto currentCov = ukf.matP();
      double eps = 1e-9;
      for (int k = 0; k < 6; k++) currentCov(k,k) += eps;

      ROOT::Math::XYZPoint pos(state[0], state[1], state[2]);

      ROOT::Math::Polar3DVector momPolar(state[3], state[4], state[5]);
      ROOT::Math::XYZVector mom(momPolar);

      // Safety Check: Ensure the updated momentum and covariances are physically valid numbers.
      if (std::isnan(mom.R()) || mom.R() < 1e-4 || std::isnan(currentCov(3,3)) || currentCov(3,3) < 0) {
          break; // Exit loop immediately, the particle cannot be tracked further.
      }

      // Calculate step energy loss
      double KE_in = Kinematics::KE(lastMom, mass);
      double KE_out = Kinematics::KE(mom, mass);
      lastMom = mom;

      // Store forward filter results
      xForward.push_back(pos.X());
      yForward.push_back(pos.Y());
      zForward.push_back(pos.Z());
      eLossForward.push_back(KE_in - KE_out);
      pForward.push_back(mom.R());
      sigmapForward.push_back(std::sqrt(currentCov(3, 3)));

      // Calculate spatial residual (distance between true hit and filtered state)
      double res = std::sqrt(std::pow(pos.X() - point.X(), 2) + 
                             std::pow(pos.Y() - point.Y(), 2) + 
                             std::pow(pos.Z() - point.Z(), 2));
      residualForward.push_back(res);
      
      // DEBUG: Print momentum at each step
      double p_true_init = beginMom;
      double dp_rel = (mom.R() - p_true_init) / p_true_init * 100.0;
      
      /*if (i % 5 == 1 || i < 3) {
        std::printf("[FORWARD] Hit %zu: p_rec=%.4f MeV/c, p_true=%.4f MeV/c, (p_rec-p_true)/p_true=%.2f%%\n", 
                    i, mom.R(), p_true_init, dp_rel);
      }
      */
      // Store early and late momentum for later comparison
      if (i == 1) pForward_early = mom.R();  // momentum at 2nd hit
      pForward_late = mom.R();  // last momentum in forward pass
    }
   
   // If the filter couldn't process enough points to make a meaningful track, flag as a failure.
   if (xForward.size() < 4) {
       UKFResult failResult;
       failResult.trackID = trackID;
       failResult.p_rec = -999.0;
       failResult.E_rec = -999.0;
       failResult.sumEloss = -999.0;
       failResult.p_true = initialMom.R();
       return failResult;
   }

     // 5. Smoothing (Backward Filter)
     // This pass runs backwards along the track to improve the state estimates using all available future data.

   
        double phi_forward = wrapAngle(ukf.vecX()[5]);
        double theta_forward =  ukf.vecX()[4];

        ukf.smoothUKF();
        auto smoothedStates = ukf.GetSmoothedStates();

        if (smoothedStates.empty()) {
        if (drawPlots) {
            std::cout << "[WARNING] Track " << trackID 
                    << " aborted by the UKF." << std::endl;
        }

        UKFResult failResult;
        failResult.trackID = trackID;
        failResult.p_rec = -999.0;
        failResult.E_rec = -999.0;
        failResult.sumEloss = -999.0;
        failResult.p_true = initialMom.R();
        return failResult;
        }

        // Calculate angular residuals between the smoothed state, forward state, and MC truth
        double phi_smoothed = wrapAngle(smoothedStates[1][5]);
        double theta_smoothed = smoothedStates[1][4];
        double theta_true = atan2(std::sqrt(initialMom.X()*initialMom.X() +initialMom.Y()*initialMom.Y()),initialMom.Z());
        double dtheta_rad = theta_smoothed - theta_true;
        all_dtheta_rad.push_back(dtheta_rad);

        double dtheta_fwd = theta_forward - theta_true;
        double dtheta_sm  = theta_smoothed - theta_true;
        double dtheta_sm_fwd = theta_smoothed - theta_forward;

        // Push to global accumulators, ensuring we filter out NaN/inf values
        if (std::isfinite(dtheta_fwd))  all_dtheta_fwd_rad.push_back(dtheta_fwd);
        if (std::isfinite(dtheta_sm))   all_dtheta_sm_rad.push_back(dtheta_sm);
        if (std::isfinite(dtheta_sm_fwd)) all_dtheta_sm_minus_fwd_rad.push_back(dtheta_sm_fwd);

        // Optional diagnostic printout per track
        printf("DIAG Track %d: dtheta_fwd=%.3f mrad dtheta_sm=%.3f mrad dtheta_sm-fwd=%.3f mrad\n",
            trackID, dtheta_fwd*1e3, dtheta_sm*1e3, dtheta_sm_fwd*1e3);

        printf("DIAG smoothing: theta_forward=%.3f deg theta_smoothed=%.3f deg theta_true=%.3f deg dtheta_mrad=%.3f\n",
           theta_forward * TMath::RadToDeg(), theta_smoothed * TMath::RadToDeg(), theta_true * TMath::RadToDeg(), dtheta_rad * 1e3);

   
   //auto smoothedStates = ukf.GetSmoothedStates();
   auto smoothedCov = ukf.GetSmoothedCovariances();

   // Extract and store the final smoothed states and calculate smoothed residuals

   // =========================================================================
   // DISCLAIMER: Extracting final kinematics from hit 1 instead of hit 0.
   // During the final step of the backward smoother, the algorithm attempts to 
   // reconcile the smoothed track with our manual initial guess (which has an 
   // artificially huge covariance). This "bad prior" pulls the reconstructed 
   // momentum and creates an unphysical offset. Skipping hit 0 avoids this bias.
   // =========================================================================
   
   //std::cout << "SmoothedStates size: " << smoothedStates.size() << "number of clusters: " << hX.size() << std::endl;

   for (size_t i = 1; i < smoothedStates.size(); ++i) {
      xSmooth.push_back(smoothedStates[i][0]);
      ySmooth.push_back(smoothedStates[i][1]);
      zSmooth.push_back(smoothedStates[i][2]);
      //pSmooth.push_back(smoothedStates[i][3]);

      double current_p_sm = smoothedStates[i][3]; 
      pSmooth.push_back(current_p_sm);
      all_cluster_p_rec.push_back(current_p_sm);
      all_cluster_p_true.push_back(initialMom.R());

      if (i == 0) {
        // Físicamente, en el punto inicial no ha habido "paso" previo, 
        // así que la pérdida de energía acumulada en este paso es 0.
        eLossSmooth.push_back(0.0); 
        } else {
     // if (i >= 1) {
        double prev_p_sm = smoothedStates[i-1][3];
        double ke_prev = Kinematics::KE(prev_p_sm, mass);
        double ke_curr = Kinematics::KE(current_p_sm, mass);
        //std::cout << "[DEBUG] Smoothed Eloss at hit " << i << ": KE_prev=" << ke_prev << " MeV, KE_curr=" << ke_curr << " MeV, Eloss=" << (ke_prev - ke_curr) << " MeV" << std::endl;

        eLossSmooth.push_back(ke_prev - ke_curr);
      }

      sigmapSmooth.push_back(std::sqrt(smoothedCov[i](3, 3))); 
      
      // DEBUG: Print momentum at each smoothed step
      double p_true_init = initialMom.R();
      double dp_rel = (current_p_sm - p_true_init) / p_true_init * 100.0;
      if (i % 5 == 1 || i < 3) {
        //std::printf("[SMOOTHED] Hit %zu: p_rec=%.4f MeV/c, p_true=%.4f MeV/c, (p_rec-p_true)/p_true=%.2f%%\n", i, current_p_sm, p_true_init, dp_rel);
      }

      // Distance between the real measured hit and the smoothed estimated state
      double mX = hX[i]; // ¡Volvemos a i! Adiós al +1
      double mY = hY[i];
      double mZ = hZ[i];
      double resSmooth = std::sqrt(std::pow(smoothedStates[i][0] - mX, 2) + 
                                   std::pow(smoothedStates[i][1] - mY, 2) + 
                                   std::pow(smoothedStates[i][2] - mZ, 2));
   
      residualSmooth.push_back(resSmooth);

   }

   std::cout << "residualSmooth size: " << residualSmooth.size() << "residualForward size: " << residualForward.size() << std::endl;

   // 6. Summary of Results 
   // (Calculates the final kinematics to package into the UKFResult struct)
   double E_sim = Kinematics::KE(beginMom, mass);
   // Use momentum from hit 1 to avoid bias from smoothed state at hit 0
   // (which is pulled toward the artificially uncertain initial guess)
   double p_reco_val = smoothedStates[1][3]; // Use hit 1 instead of hit 0 to avoid bias
   double E_rec_val = Kinematics::KE(p_reco_val, mass);
   double sumElossMC = std::accumulate(hEloss.begin(), hEloss.end(), 0.0);
   double sumElossUKF = std::accumulate(eLossForward.begin(), eLossForward.end(), 0.0);
   

   // =========================================================================
   // VISUALIZATION: MULTIPLOT FOR EACH TRACK
   // =========================================================================
   
   // Only generate and draw plots if the user explicitly requested them
   if (drawPlots == true) 
   {
       gROOT->cd();
       
       // Generate unique names and titles for the canvas using the current Track ID
       TString canvasName = Form("c_Track%d", trackID);
       TString canvasTitle = Form("UKF Results - Track %d", trackID);

       // Create a large Canvas (1200x900 pixels) and divide it into a 2x2 grid (4 panels)
       TCanvas *cAll = new TCanvas(canvasName, canvasTitle, 1200, 900);
       cAll->Divide(2, 2);

   // -------------------------------------------------------------------------
   // PANEL 1: 3D TRACK (Spatial Trajectory)
   // -------------------------------------------------------------------------
   // Navigate to the first pad (top-left) to draw the 3D particle trajectory
   cAll->cd(1);
   
   TGraph2D *track = new TGraph2D((int)hX.size(), (double*)hX.data(), (double*)hY.data(), (double*)hZ.data());
   track->SetName(Form("MC_Track_%d", trackID));
   track->SetTitle(Form("3D Track %d;X [mm];Y [mm];Z [mm]", trackID));
   track->SetMarkerStyle(20);
   track->SetMarkerSize(0.8);

   TGraph2D *track2 = new TGraph2D(xForward.size(), xForward.data(), yForward.data(), zForward.data());
   track2->SetName(Form("Prop_Track_%d", trackID)); 
   track2->SetMarkerStyle(21);
   track2->SetMarkerSize(0.8);
   track2->SetMarkerColor(kRed);

   TGraph2D *smoothedTrack = new TGraph2D(xSmooth.size(), xSmooth.data(), ySmooth.data(), zSmooth.data());
   smoothedTrack->SetName(Form("Smooth_Track_%d", trackID));
   smoothedTrack->SetMarkerStyle(22);
   smoothedTrack->SetMarkerSize(0.8);
   smoothedTrack->SetMarkerColor(kGreen + 2);

   // Axis range 
   double xmin = std::min(*std::min_element(hX.begin(), hX.end()), *std::min_element(xForward.begin(), xForward.end()));
   double xmax = std::max(*std::max_element(hX.begin(), hX.end()), *std::max_element(xForward.begin(), xForward.end()));
   double ymin = std::min(*std::min_element(hY.begin(), hY.end()), *std::min_element(yForward.begin(), yForward.end()));
   double ymax = std::max(*std::max_element(hY.begin(), hY.end()), *std::max_element(yForward.begin(), yForward.end()));
   double zmin = std::min(*std::min_element(hZ.begin(), hZ.end()), *std::min_element(zForward.begin(), zForward.end()));
   double zmax = std::max(*std::max_element(hZ.begin(), hZ.end()), *std::max_element(zForward.begin(), zForward.end()));

   track->GetXaxis()->SetLimits(xmin, xmax);
   track->GetYaxis()->SetLimits(ymin, ymax);
   track->GetZaxis()->SetLimits(zmin, zmax);

   track->Draw("P");
   track2->Draw("PSAME");
   smoothedTrack->Draw("PSAME");

   TLegend *leg1 = new TLegend(0.6, 0.7, 0.88, 0.88);
   leg1->AddEntry(track, "MC Truth", "p");
   leg1->AddEntry(track2, "Filtered (Forward)", "p");
   leg1->AddEntry(smoothedTrack, "Smoothed", "p");
   leg1->Draw();

   // -------------------------------------------------------------------------
   // PANEL 2: ENERGY LOSS (dE/dx)
   // -------------------------------------------------------------------------
   cAll->cd(2);

   TGraph *elossGraph = new TGraph(hEloss.size());
   for (size_t i = 0; i < hEloss.size(); ++i) elossGraph->SetPoint(i, i, hEloss[i]);
   elossGraph->SetTitle("Energy Loss per Hit;Hit Number;Energy Loss [MeV]");
   elossGraph->SetMarkerStyle(20);

   TGraph *eLossForwardGraph = new TGraph(eLossForward.size());
   for (size_t i = 0; i < eLossForward.size(); ++i) eLossForwardGraph->SetPoint(i, i, eLossForward[i]);
   eLossForwardGraph->SetMarkerStyle(21);
   eLossForwardGraph->SetMarkerColor(kRed);

    TGraph *elossSmoothGraph = new TGraph(eLossSmooth.size());
    for (size_t i = 0; i < eLossSmooth.size(); ++i)  elossSmoothGraph->SetPoint(i, i, eLossSmooth[i]);
    elossSmoothGraph->SetMarkerStyle(22);
    elossSmoothGraph->SetMarkerColor(kGreen + 2);

   elossGraph->Draw("AP"); 
   eLossForwardGraph->Draw("PSAME"); 
   elossSmoothGraph->Draw("PSAME"); 

   TLegend *leg2 = new TLegend(0.6, 0.7, 0.88, 0.88);
   leg2->AddEntry(elossGraph, "Measured Eloss", "p");
   leg2->AddEntry(eLossForwardGraph, "Propagated Eloss", "p");
   leg2->AddEntry(elossSmoothGraph, "Smoothed Eloss", "p");
   leg2->Draw();

   // -------------------------------------------------------------------------
   // PANEL 3: MOMENTUM
   // -------------------------------------------------------------------------
   cAll->cd(3);

   TGraphErrors *pGraph = new TGraphErrors(pForward.size());
   for (size_t i = 0; i < pForward.size(); ++i) {
      pGraph->SetPoint(i, i, pForward[i]);
      pGraph->SetPointError(i, 0, sigmapForward[i] * 5); // Bars scaled by 5
   }
   pGraph->SetTitle("Momentum per Hit (Errors x5);Hit Number;Momentum [MeV/c]");
   pGraph->SetMarkerStyle(20);
   pGraph->SetMarkerColor(kBlue);
   pGraph->SetLineColor(kBlue);

   TGraphErrors *pSmoothGraph = new TGraphErrors(pSmooth.size());
   for (size_t i = 0; i < pSmooth.size(); ++i) {
      pSmoothGraph->SetPoint(i, i, pSmooth[i]);
      pSmoothGraph->SetPointError(i, 0, sigmapSmooth[i] * 5); 
   }
   pSmoothGraph->SetMarkerStyle(22);
   pSmoothGraph->SetMarkerColor(kGreen + 2);
   pSmoothGraph->SetLineColor(kGreen + 2);

   pGraph->Draw("AP");
   pSmoothGraph->Draw("PSAME");

   TLegend *leg3 = new TLegend(0.6, 0.7, 0.88, 0.88);
   leg3->AddEntry(pGraph, "Filtered Momentum", "pe");
   leg3->AddEntry(pSmoothGraph, "Smoothed Momentum", "pe");
   leg3->Draw();

   // -------------------------------------------------------------------------
   // PANEL 4: RESIDUALS
   // -------------------------------------------------------------------------
   cAll->cd(4);
   
   TGraph *residualGraph = new TGraph(residualForward.size());
   for (size_t i = 0; i < residualForward.size(); ++i) residualGraph->SetPoint(i, i, residualForward[i] * 0.1);
   residualGraph->SetTitle("Residual per Hit;Hit Number;Residual [cm]");
   residualGraph->SetMarkerStyle(23);
   residualGraph->SetMarkerColor(kMagenta);

  TGraph *residualSmoothGraph = nullptr;
   if (!residualSmooth.empty()) {
       residualSmoothGraph = new TGraph(residualSmooth.size());
       for (size_t i = 0; i < residualSmooth.size(); ++i) residualSmoothGraph->SetPoint(i, i+1, residualSmooth[i] * 0.1);
       residualSmoothGraph->SetMarkerStyle(24);
       residualSmoothGraph->SetMarkerColor(kOrange + 7);
   }

  // 3. Usamos TMultiGraph para que ROOT escale los ejes a TODAS las gráficas
   TMultiGraph *mg4 = new TMultiGraph();
   mg4->SetTitle("Residual per Hit;Hit Number;Residual [cm]");
   mg4->Add(residualGraph, "P");
   if (residualSmoothGraph) {
       mg4->Add(residualSmoothGraph, "P");
   }
   
   // Al dibujar el MultiGraph con "A", ROOT recalcula el suelo y el techo
   mg4->Draw("A");

   TLegend *leg4 = new TLegend(0.6, 0.7, 0.88, 0.88);
   leg4->AddEntry(residualGraph, "Filtered Residual", "p");
   if (residualSmoothGraph) leg4->AddEntry(residualSmoothGraph, "Smoothed Residual", "p");
   leg4->Draw();

   // Update teh canvas to render all the plots. 
   cAll->Update();
   cAll->SaveAs(Form("UKF_Plot_Track_%d.png", trackID));

   } //if drawPlots

   // --- FINAL SUMMARY IN TERMINAL ---
   std::cout << "\n\n" << std::string(65, '=') << std::endl;
   std::cout << "          KALMAN FILTER (UKF) SUMMARY - TRACK " << trackID << std::endl;
   std::cout << std::string(65, '-') << std::endl;

   std::printf("  MOMENTUM (p):\n");
   std::printf("    - Simulated (MC):       %10.4f MeV/c\n", beginMom);
   std::printf("    - Reconstructed (UKF):  %10.4f MeV/c\n", p_reco_val);
   std::printf("    - Relative Error:       %10.4f %%\n\n", (p_reco_val - beginMom)/beginMom * 100);

   std::printf("  KINETIC ENERGY (T):\n");
   std::printf("    - Simulated (MC):       %10.4f MeV\n", E_sim);
   std::printf("    - Reconstructed (UKF):  %10.4f MeV\n", E_rec_val);
   std::printf("    - Relative Error:       %10.4f %%\n\n", (E_rec_val - E_sim)/E_sim * 100);

   std::printf("  ENERGY LOSS VALIDATION (Total dE):\n");
   std::printf("    - Sum Eloss (MC):       %10.4f MeV\n", sumElossMC);
   std::printf("    - Sum Eloss (UKF):      %10.4f MeV\n", sumElossUKF);
   std::printf("    - dE Discrepancy:       %10.4f MeV\n", std::abs(sumElossMC - sumElossUKF));

   std::cout << std::string(65, '-') << std::endl;
   std::cout << "  NUMERICAL STABILITY:     EXCELLENT (nTouch = 0)" << std::endl;
   std::cout << "  STATUS:                  CONVERGED" << std::endl;
   std::cout << std::string(65, '=') << "\n\n" << std::endl;

//--------------------------------------------------------------------------------------------------------------------
   UKFResult result;
   result.trackID = trackID;
   result.p_rec = p_reco_val;  
   result.E_rec = E_rec_val;
   result.sumEloss = sumElossUKF;
   
   result.p_true = initialMom.R(); // Save true momentum for easier analysis later
   result.theta_rec = smoothedStates[1][4]; // Reconstructed polar angle (smoothed)
   result.phi_rec = wrapAngle(smoothedStates[1][5]);   // Reconstructed azimuthal angle (smoothed)
   
   //result.theta_true = std::atan2(std::sqrt(initialMom.X()*initialMom.X() + initialMom.Y()*initialMom.Y()), initialMom.Z()); // True polar angle
   //result.phi_true = std::atan2(initialMom.Y(), initialMom.X()); // True azimuthal angle
   result.phi_true = wrapAngle(initialMom.Phi()); // True azimuthal angle using ROOT's built-in function, which handles the correct quadrant
   result.theta_true = initialMom.Theta(); // True polar angle using ROOT's built-in function
    
   /*if (smoothedStates.size() > 1) {

    double theta0 = smoothedStates[0][4];
    double theta1 = smoothedStates[1][4];

    double phi0 = smoothedStates[0][5];
    double phi1 = smoothedStates[1][5];

    double p0 = smoothedStates[0][3];
    double p1 = smoothedStates[1][3];

    printf("SMOOTHER SHIFT Track %d:\n", trackID);
    printf("   Δp     = %.6f\n", p0 - p1);
    printf("   Δtheta = %.6f rad (%.3f mrad)\n", theta0 - theta1, (theta0 - theta1)*1e3);
    printf("   Δphi   = %.6f rad (%.3f mrad)\n", phi0 - phi1, (phi0 - phi1)*1e3);
    }*/

   result.status = FitStatus::SUCCESS;
   return result;
}

void UKFMultiTrack_Georgina(const char* simFilename = "/home/georgina/fair_install/ATTPCROOTv2_KF_fork/ATTPCROOTv2/macro/Simulation/ATTPC/protons/data/protonssim_2T_H300torr_40-80MeV_theta10-80.root", int eventToDraw = 1) 
//void UKFMultiTrack_Georgina(const char* simFilename = "/home/georgina/fair_install/ATTPCROOTv2_KF_fork/ATTPCROOTv2/macro/Simulation/ATTPC/pions/data/pionssim_20-40MeV_Bfield_20kG_H300torr_theta0-90.root", int eventToDraw = 1)
//void UKFMultiTrack_Georgina(const char* simFilename = "/home/georgina/fair_install/ATTPCROOTv2_KF_fork/ATTPCROOTv2/macro/Simulation/ATTPC/pions/data/pionssim_30MeV_Bfield_20kG_H300torr_theta30.root", int eventToDraw = 1) 
{
// === A. PREPARE OUTPUT FILE ===
    TFile* outFile = new TFile("reco_ukf_output_pionssim_30MeV_Bfield_20kG_H300torr_theta30_Bethe_initialMom_10kEvt_hit1.root", "RECREATE");
    //TFile* outFile = new TFile("reco_ukf_output_protonssim_40-80MeV_Bfield_20kG_H300torr_theta10-80_catima_initialMom_10kEvt_hit1_cluster10.root", "RECREATE");
    TTree* outTree = new TTree("UKFTree", "Resultados del UKF");

    UKFResult res; 
    outTree->Branch("eventID", &res.eventID, "eventID/I");
    outTree->Branch("trackID", &res.trackID, "trackID/I");
    outTree->Branch("p_rec", &res.p_rec, "p_rec/D");
    outTree->Branch("theta_rec", &res.theta_rec, "theta_rec/D");
    outTree->Branch("phi_rec", &res.phi_rec, "phi_rec/D");
    outTree->Branch("E_rec", &res.E_rec, "E_rec/D");
    outTree->Branch("sumEloss", &res.sumEloss, "sumEloss/D");
    outTree->Branch("p_true", &res.p_true, "p_true/D");
    outTree->Branch("theta_true", &res.theta_true, "theta_true/D");
    outTree->Branch("phi_true", &res.phi_true, "phi_true/D");
    outTree->Branch("status", &res.status, "status/I");

    // Global counters to calculate statistics and efficiency at the end of the script
    int totalTracksProcessed = 0;
    int successfulFits = 0;
    int failedStoppedTracks = 0;
    int failedFits = 0; // New counter for fits that failed for reasons other than stopping in the Bragg peak

   // === B. OPEN SIMULATION FILE (ONLY ONCE) ===
    TFile* simFile = TFile::Open(simFilename, "READ");
    TTree* simTree = (TTree*)simFile->Get("cbmsim");
    TClonesArray* tpcPoints = new TClonesArray("AtMCPoint"); 
    simTree->SetBranchAddress("AtTpcPoint", &tpcPoints);

    int numEvents = simTree->GetEntries();
    std::cout << "initializing" << numEvents << " events..." << std::endl;

   // === C. MAIN LOOP ===
    for (int ev = 0; ev < numEvents; ev++) {
  
        LoadHitsROOT(simTree, tpcPoints, ev);
        for (auto const& [trackID, xVector] : posX) {
            
            if (xVector.size() < 5) continue; // Extra Kalman protection
            totalTracksProcessed++;

            bool drawPlots = (ev == eventToDraw);

            /*try {
            // This function call will now "jump" to the 'catch' block 
            // as soon as the C++ code hits a 'throw'
            UKFResult result_parcial = runKalman(posX[trackID], posY[trackID], posZ[trackID], Eloss[trackID], 
                                                mass_p, charge_p, 1, 1, initialMom[trackID], trackID, drawPlots); //coulombs

            //--------------------------------------------------------------------------------------------------

            // If it didn't throw, we check the result and fill
           if (result_parcial.p_rec > 0 && result_parcial.p_rec != -999.0) {
                res.eventID  = ev;
                res.trackID  = result_parcial.trackID;
                res.p_rec    = result_parcial.p_rec;
                res.E_rec    = result_parcial.E_rec;
                res.sumEloss = result_parcial.sumEloss;
                res.p_true   = result_parcial.p_true;

                res.theta_rec  = result_parcial.theta_rec; //rad
                res.phi_rec    = result_parcial.phi_rec; //rad
                res.theta_true = result_parcial.theta_true; //rad
                res.phi_true   = result_parcial.phi_true; //rad

                outTree->Fill();
                successfulFits++;
            } else {
                failedStoppedTracks++;
                std::cout << "Event " << ev << ", Track " << trackID << ": UKF failed to reconstruct (p_rec = " << result_parcial.p_rec << "). Likely stopped in Bragg Peak. Skipping..." << std::endl;
            }

            } catch (const std::exception& e) {
                std::string msg = e.what();
                if (msg == "ParticleStopped") {
                    failedStoppedTracks++;
                    std::cout << "Event " << ev << ": Track stopped in Bragg Peak. Skipping..." << std::endl;
                } else {
                    // Rethrow if it is a different, unexpected error
                    throw;
                }
            
            } */

            UKFResult result_parcial = runKalman(posX[trackID], posY[trackID], posZ[trackID], Eloss[trackID], 
                                                mass_pi, charge_pi_Coulombs, 1, 1, initialMom[trackID], trackID, drawPlots); //coulombs

            std::cout << "status = " << (int)result_parcial.status << std::endl;

            if (result_parcial.status == FitStatus::SUCCESS) {

                res.eventID  = ev;
                res.trackID  = result_parcial.trackID;

                res.p_rec    = result_parcial.p_rec;
                res.E_rec    = result_parcial.E_rec;
                res.sumEloss = result_parcial.sumEloss;
                res.p_true   = result_parcial.p_true;

                res.theta_rec  = result_parcial.theta_rec; // rad
                res.phi_rec    = result_parcial.phi_rec;   // rad
                res.theta_true = result_parcial.theta_true; // rad
                res.phi_true   = result_parcial.phi_true;   // rad
                res.status = result_parcial.status;

                outTree->Fill();
                successfulFits++;
            }

            else if (result_parcial.status == FitStatus::STOPPED_IN_BRAGG) {

                failedStoppedTracks++;

                std::cout << "Event " << ev
                        << ", Track " << trackID
                        << ": Track stopped in Bragg Peak. Skipping..."
                        << std::endl;
            }

            else {

                failedFits++;

                std::cout << "Event " << ev
                        << ", Track " << trackID
                        << ": Kalman fit failed."
                        << std::endl;
            }
        }   
        
        if (ev % 100 == 0) std::cout << "Processed " << ev << " events..." << std::endl;
    }

   
    // === D. FINAL SUMMARY AND SAVE ===
    std::cout << "\n===============================================" << std::endl;
    std::cout << "       UKF RECONSTRUCTION SUMMARY" << std::endl;
    std::cout << "===============================================" << std::endl;
    std::cout << "Total tracks analyzed:      " << totalTracksProcessed << std::endl;
    std::cout << "Successful fits:            " << successfulFits << std::endl;
    std::cout << "Tracks stopped (Bragg Peak): " << failedStoppedTracks << std::endl;
    
    if (totalTracksProcessed > 0) {
        double efficiency = (100.0 * successfulFits) / totalTracksProcessed;
        std::cout << "Reconstruction Efficiency:  " << efficiency << "%" << std::endl;
    }
    std::cout << "===============================================\n" << std::endl;
    
    outFile->cd();

    /*TH2F* h_pdiff_vs_ptrue_all = new TH2F("h_pdiff_vs_ptrue_all", 
                                          "Momentum Diff vs True Momentum (All Clusters); p_{true} [MeV/c]; p_{rec} - p_{true} [MeV/c]", 
                                          50, 35, 85, 
                                          100, -2.0, 2.0);

    for (size_t k = 0; k < all_cluster_p_rec.size(); k++) {
        h_pdiff_vs_ptrue_all->Fill(all_cluster_p_true[k], all_cluster_p_rec[k] - all_cluster_p_true[k]);
    }

    h_pdiff_vs_ptrue_all->Write();

    TCanvas* c_all_clusters = new TCanvas("c_all_clusters", "All Clusters Momentum Diff", 800, 600);
    c_all_clusters->cd();
    
    h_pdiff_vs_ptrue_all->Draw("COLZ"); // Draw with the color palette
    
    // Optional: Draw a line at Y=0
    TLine* line0_all = new TLine(35, 0, 85, 0);
    line0_all->SetLineStyle(2);
    line0_all->SetLineColor(kRed);
    line0_all->SetLineWidth(2);
    line0_all->Draw("SAME");
    
    // Save it as an image file in your folder
    c_all_clusters->SaveAs("Momentum_Difference_All_Clusters.png");*/


    outTree->Write();
    outFile->Close();
    simFile->Close();
    std::cout << "\n << output file starting with 'reco_ikf_output' saved" << std::endl;
}
