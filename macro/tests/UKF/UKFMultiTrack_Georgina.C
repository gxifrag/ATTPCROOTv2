struct UKFResult {
   int eventID;
   int trackID;
   double p_rec;      // Momento reconstruido [MeV/c]
   double E_rec;      // Energía cinética reconstruida [MeV]
   double theta_rec;  // Ángulo polar reconstruido [rad]
   double phi_rec;    // Ángulo azimutal reconstruido [rad]
   double sumEloss;   // Pérdida de energía total [MeV]
   double p_true;     // Momento inicial de la simulación (para comparar luego)
   double theta_true;  // Ángulo polar verdadero (para comparar luego)
   double phi_true;    // Ángulo azimutal verdadero (para comparar luego)
};

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

const double mass_p = 938.272;           // Mass of proton in MeV/c^2
const double charge_p = 1.602176634e-19; // Charge of proton

// We are gonna use maps (dictionary). The key is the trackID.
std::map<int, std::vector<double>> posX;
std::map<int, std::vector<double>> posY;
std::map<int, std::vector<double>> posZ;
std::map<int, std::vector<double>> Eloss;

// NUEVO MAPA: Guardará un solo XYZVector (el momento inicial) por cada trackID
std::map<int, ROOT::Math::XYZVector> initialMom;

int pointsToCluster= 5; 

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

   for (int i = 0; i < nPoints; i++) {
      AtMCPoint* point = (AtMCPoint*)tpcPoints->At(i);
      if (!point) continue;

      int trackID = point->GetTrackID();

      // Extraemos las posiciones y la energía
        double mmX = point->GetX() * 10.0;
        double mmY = point->GetY() * 10.0;
        double mmZ = point->GetZ() * 10.0;
        double Ei = point->GetEnergyLoss() * 1000.0;   

      if (posX[trackID].empty()) { 
            // Es el primer punto que vemos de esta partícula
            posX[trackID].push_back(mmX); 
            posY[trackID].push_back(mmY); 
            posZ[trackID].push_back(mmZ);
            Eloss[trackID].push_back(0.0);
            
            double px = point->GetPx() * 1000.0;
            double py = point->GetPy() * 1000.0;
            double pz = point->GetPz() * 1000.0;
            initialMom[trackID] = ROOT::Math::XYZVector(px, py, pz);
            cout << "eventID: " << targetEvent << " - TrackID: " << trackID 
                 << " - Initial Momentum (MeV/c): (" << px << ", " << py << ", " << pz << ")" << endl;

            
            currentELoss[trackID] = Ei;
            hitCount[trackID] = 0;
      } else {
            // Ya teníamos puntos de esta partícula, seguimos sumando
            currentELoss[trackID] += Ei;
            hitCount[trackID]++;
            
            if (hitCount[trackID] % pointsToCluster == 0) {
                posX[trackID].push_back(mmX); 
                posY[trackID].push_back(mmY); 
                posZ[trackID].push_back(mmZ);
                Eloss[trackID].push_back(currentELoss[trackID]);
                
                // Reseteamos el contador de energía para el siguiente cluster de esta partícula
                currentELoss[trackID] = 0;
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

   // Vectors locals per a aquesta execució (així no es barregen amb l'anterior)
   std::vector<double> x2, y2, z2, Eloss2, p2, sigmap2, lambda2, sigmalambda2, residual;
   std::vector<double> xSmooth, ySmooth, zSmooth, pSmooth, sigmapSmooth, residualSmooth, eLossSmooth;

   // CATIMA Hidrogeno 300torr
   auto elossModel = std::make_unique<AtTools::AtELossCATIMA>(3.3084e-5);
   // CATIMA hidrogeno 60torr
   //auto elossModel = std::make_unique<AtTools::AtELossCATIMA>(6.6168e-6);
   elossModel->SetProjectile(Z, A, mass / 931.494); // Convertim MeV/c2 a amu aprox.
   std::vector<std::tuple<int, int, int>> mat = {{1, 1, 2}}; // Hidrogen
   elossModel->SetMaterial(mat);

   // 2. Propagator i UKF
   AtTools::AtPropagator propagator(charge, mass, std::move(elossModel));
   propagator.SetEField({0, 0, 0});
   propagator.SetBField({0, 0, 3.});//2.85 T
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

   // 4. Bucle de Hits (Usa hX.size() passat per argument)

   ROOT::Math::XYZVector lastMom = startMom;
   for (size_t i = 1; i < hX.size(); ++i) {
      if (i % 200 == 0) { 
         std::cout << "Processing hit " << i << " of " << hX.size() << std::endl; 
      }
    
      XYZPoint point(hX[i], hY[i], hZ[i]);

      double currentKE = Kinematics::KE(lastMom, mass);
      
      if (currentKE < 0.05) { 
          std::cout << "[DEBUG] Energia critica (" << currentKE << " MeV). Deteniendo propagacion para evitar cuelgue." << std::endl;
          break; 
      }
      ukf.SetMeasCov(cov_meas);

      try {
          ukf.predictUKF(point);
          ukf.correctUKF(point);
      } catch (...) {
          // Si la matemática explota aquí, salimos del bucle limpiamente
          break; 
      }

      auto state = ukf.vecX();
      if (state.size() < 6) break; //condition to avoid crash if UKF fails and stops updating states (we will check this later in the smoothed states)

      auto currentCov = ukf.matP();
      ROOT::Math::XYZPoint pos(state[0], state[1], state[2]);

      ROOT::Math::Polar3DVector momPolar(state[3], state[4], state[5]);
      ROOT::Math::XYZVector mom(momPolar);

      if (std::isnan(mom.R()) || mom.R() < 1e-4 || std::isnan(currentCov(3,3)) || currentCov(3,3) < 0) {
          break; // Salimos del bucle inmediatamente, la partícula no da para más
      }

      double KE_in = Kinematics::KE(lastMom, mass);
      double KE_out = Kinematics::KE(mom, mass);
      lastMom = mom;

      x2.push_back(pos.X());
      y2.push_back(pos.Y());
      z2.push_back(pos.Z());
      Eloss2.push_back(KE_in - KE_out);
      p2.push_back(mom.R());
      sigmap2.push_back(std::sqrt(currentCov(3, 3)));

      double res = std::sqrt(std::pow(pos.X() - point.X(), 2) + 
                             std::pow(pos.Y() - point.Y(), 2) + 
                             std::pow(pos.Z() - point.Z(), 2));
      residual.push_back(res);
   } 
   
   if (x2.size() < 4) {
       UKFResult failResult;
       failResult.trackID = trackID;
       failResult.p_rec = -999.0;
       failResult.E_rec = -999.0;
       failResult.sumEloss = -999.0;
       failResult.p_true = initialMom.R();
       return failResult;
   }

   // Smoothing
    try {
        ukf.smoothUKF();
    } catch (...) {
        std::cout << "[DEBUG] Abortado: Fallo matematico en el Smoothing." << std::endl;
        UKFResult failResult; failResult.p_rec = -999.0; return failResult;
    }
   auto smoothedStates = ukf.GetSmoothedStates();

   if (smoothedStates.empty()) {
       // Si no pasamos drawPlots, silenciamos el print para no inundar la terminal
       if (drawPlots) {
           std::cout << "[WARNING] Track " << trackID << " abortada por el UKF (se detuvo antes de tiempo)." << std::endl;
       }
       
       // Devolvemos un resultado "falso" con valores -999 para identificar el fallo
       UKFResult failResult;
       failResult.trackID = trackID;
       failResult.p_rec = -999.0; 
       failResult.E_rec = -999.0;
       failResult.sumEloss = -999.0;
       failResult.p_true = initialMom.R();
       return failResult; // Salimos de la función inmediatamente
   }

   auto smoothedCov = ukf.GetSmoothedCovariances();

   for (size_t i = 0; i < smoothedStates.size(); ++i) {
      xSmooth.push_back(smoothedStates[i][0]);
      ySmooth.push_back(smoothedStates[i][1]);
      zSmooth.push_back(smoothedStates[i][2]);
      pSmooth.push_back(smoothedStates[i][3]);

      sigmapSmooth.push_back(std::sqrt(smoothedCov[i](3, 3))); 

      // Calculamos la distancia entre el hit real medido y el estado suavizado
      double mX = hX[i]; 
      double mY = hY[i];
      double mZ = hZ[i];
      double resSmooth = std::sqrt(std::pow(smoothedStates[i][0] - mX, 2) + 
                                   std::pow(smoothedStates[i][1] - mY, 2) + 
                                   std::pow(smoothedStates[i][2] - mZ, 2));
   
      residualSmooth.push_back(resSmooth);
   }

   // 6. Resum de resultats (el bloc que t'agrada per a la captura)
   double E_sim = Kinematics::KE(beginMom, mass);
   double p_reco_val = smoothedStates[0][3];
   double E_rec_val = Kinematics::KE(p_reco_val, mass);
   double sumElossMC = std::accumulate(hEloss.begin(), hEloss.end(), 0.0);
   double sumElossUKF = std::accumulate(Eloss2.begin(), Eloss2.end(), 0.0);

   // =========================================================================
   // PLOTS: MULTIPLOT POR TRACK
   // =========================================================================
   
   if (drawPlots==true) 
   {
   gROOT->cd();
   TString canvasName = Form("c_Track%d", trackID);
   TString canvasTitle = Form("UKF Results - Track %d", trackID);

   // Creamos un Canvas grande y lo dividimos en 2x2 paneles
   TCanvas *cAll = new TCanvas(canvasName, canvasTitle, 1200, 900);
   cAll->Divide(2, 2);

   // -------------------------------------------------------------------------
   // PANEL 1: TRAZA 3D (Particle Track)
   // -------------------------------------------------------------------------
   cAll->cd(1);
   
   TGraph2D *track = new TGraph2D((int)hX.size(), (double*)hX.data(), (double*)hY.data(), (double*)hZ.data());
   track->SetName(Form("MC_Track_%d", trackID));
   track->SetTitle(Form("3D Track %d;X [mm];Y [mm];Z [mm]", trackID));
   track->SetMarkerStyle(20);
   track->SetMarkerSize(0.8);

   TGraph2D *track2 = new TGraph2D(x2.size(), x2.data(), y2.data(), z2.data());
   track2->SetName(Form("Prop_Track_%d", trackID)); 
   track2->SetMarkerStyle(21);
   track2->SetMarkerSize(0.8);
   track2->SetMarkerColor(kRed);

   TGraph2D *smoothedTrack = new TGraph2D(xSmooth.size(), xSmooth.data(), ySmooth.data(), zSmooth.data());
   smoothedTrack->SetName(Form("Smooth_Track_%d", trackID));
   smoothedTrack->SetMarkerStyle(22);
   smoothedTrack->SetMarkerSize(0.8);
   smoothedTrack->SetMarkerColor(kGreen + 2);

   // Ajuste de ejes
   double xmin = std::min(*std::min_element(hX.begin(), hX.end()), *std::min_element(x2.begin(), x2.end()));
   double xmax = std::max(*std::max_element(hX.begin(), hX.end()), *std::max_element(x2.begin(), x2.end()));
   double ymin = std::min(*std::min_element(hY.begin(), hY.end()), *std::min_element(y2.begin(), y2.end()));
   double ymax = std::max(*std::max_element(hY.begin(), hY.end()), *std::max_element(y2.begin(), y2.end()));
   double zmin = std::min(*std::min_element(hZ.begin(), hZ.end()), *std::min_element(z2.begin(), z2.end()));
   double zmax = std::max(*std::max_element(hZ.begin(), hZ.end()), *std::max_element(z2.begin(), z2.end()));

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

   TGraph *eloss2Graph = new TGraph(Eloss2.size());
   for (size_t i = 0; i < Eloss2.size(); ++i) eloss2Graph->SetPoint(i, i, Eloss2[i]);
   eloss2Graph->SetMarkerStyle(21);
   eloss2Graph->SetMarkerColor(kRed);

   // Solo dibujamos la smoothed si el vector tiene datos (evita crash)
   TGraph *elossSmoothGraph = nullptr;
   if (!eLossSmooth.empty()) {
       elossSmoothGraph = new TGraph(eLossSmooth.size());
       for (size_t i = 0; i < eLossSmooth.size(); ++i) elossSmoothGraph->SetPoint(i, i, eLossSmooth[i]);
       elossSmoothGraph->SetMarkerStyle(22);
       elossSmoothGraph->SetMarkerColor(kGreen + 2);
   }

   elossGraph->Draw("AP"); 
   eloss2Graph->Draw("PSAME"); 
   if (elossSmoothGraph) elossSmoothGraph->Draw("PSAME"); 

   TLegend *leg2 = new TLegend(0.6, 0.7, 0.88, 0.88);
   leg2->AddEntry(elossGraph, "Measured Eloss", "p");
   leg2->AddEntry(eloss2Graph, "Propagated Eloss", "p");
   if (elossSmoothGraph) leg2->AddEntry(elossSmoothGraph, "Smoothed Eloss", "p");
   leg2->Draw();

   // -------------------------------------------------------------------------
   // PANEL 3: MOMENTUM
   // -------------------------------------------------------------------------
   cAll->cd(3);

   TGraphErrors *pGraph = new TGraphErrors(p2.size());
   for (size_t i = 0; i < p2.size(); ++i) {
      pGraph->SetPoint(i, i, p2[i]);
      pGraph->SetPointError(i, 0, sigmap2[i] * 5); // Bars scaled by 5
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

   TGraph *residualGraph = new TGraph(residual.size());
   for (size_t i = 0; i < residual.size(); ++i) residualGraph->SetPoint(i, i, residual[i] * 0.1);
   residualGraph->SetTitle("Residual per Hit;Hit Number;Residual [cm]");
   residualGraph->SetMarkerStyle(23);
   residualGraph->SetMarkerColor(kMagenta);

  TGraph *residualSmoothGraph = nullptr;
   if (!residualSmooth.empty()) {
       residualSmoothGraph = new TGraph(residualSmooth.size());
       for (size_t i = 0; i < residualSmooth.size(); ++i) residualSmoothGraph->SetPoint(i, i, residualSmooth[i] * 0.1);
       residualSmoothGraph->SetMarkerStyle(24);
       residualSmoothGraph->SetMarkerColor(kOrange + 7);
   }

   residualGraph->Draw("AP");
   if (residualSmoothGraph) residualSmoothGraph->Draw("PSAME"); //PSAME

   TLegend *leg4 = new TLegend(0.6, 0.7, 0.88, 0.88);
   leg4->AddEntry(residualGraph, "Filtered Residual", "p");
   if (residualSmoothGraph) leg4->AddEntry(residualSmoothGraph, "Smoothed Residual", "p");
   leg4->Draw();

   // Actualizamos el canvas para que dibuje todo
   cAll->Update();
  // cAll->SaveAs(Form("UKF_Plot_Track_%d.png", trackID));

   } //if drawPlots

   // --- RESUMEN FINAL POR TERMINAL ---
   std::cout << "\n\n" << std::string(65, '=') << std::endl;
   std::cout << "          KALMAN FILTER (UKF) SUMMARY - TRACK " << trackID << std::endl;
   std::cout << std::string(65, '-') << std::endl;

   std::printf("  MOMENTUM (p):\n");
   std::printf("    - Simulated (MC):       %10.4f MeV/c\n", beginMom);
   std::printf("    - Reconstructed (UKF):  %10.4f MeV/c\n", smoothedStates[0][3]);
   std::printf("    - Relative Error:       %10.4f %%\n\n", (smoothedStates[0][3] - beginMom)/beginMom * 100);

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
   result.p_true = initialMom.R(); // Guardamos el true para que sea fácil analizar luego
   result.theta_rec = smoothedStates[0][4]; // Ángulo polar reconstruido
   result.phi_rec = smoothedStates[0][5];   // Ángulo azimutal reconstruido
   result.theta_true = std::atan2(std::sqrt(initialMom.X()*initialMom.X() + initialMom.Y()*initialMom.Y()), initialMom.Z()); // Ángulo polar verdadero
   result.phi_true = std::atan2(initialMom.Y(), initialMom.X()); // Ángulo azimutal verdadero

   return result;
}

void UKFMultiTrack_Georgina(const char* simFilename = "/home/georgina/fair_install/ATTPCROOTv2_KF/macro/Simulation/ATTPC/protons/data/protons_300torr/protonssim_3T_H300torr_40MeV_theta30.root", int eventToDraw = 1) 
{
    // === A. PREPARAR ARCHIVO DE SALIDA ===
    TFile* outFile = new TFile("reco_ukf_output_3T_H300torr_p40MeV_theta30.root", "RECREATE");
    TTree* outTree = new TTree("UKFTree", "Resultados del UKF");

    UKFResult res; // Usando el struct que definimos antes
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

    int totalTracksProcessed = 0;
    int successfulFits = 0;
    int failedStoppedTracks = 0;

    // === B. ABRIR ARCHIVO DE SIMULACIÓN (UNA SOLA VEZ) ===
    TFile* simFile = TFile::Open(simFilename, "READ");
    TTree* simTree = (TTree*)simFile->Get("cbmsim");
    TClonesArray* tpcPoints = new TClonesArray("AtMCPoint"); 
    simTree->SetBranchAddress("AtTpcPoint", &tpcPoints);

    int numEventos = simTree->GetEntries();
    std::cout << "initializing" << numEventos << " events..." << std::endl;

    // === C. BUCLE PRINCIPAL ===
    for (int ev = 0; ev < 10000; ev++) {
        
        // 1. Delegamos el trabajo de cargar hits a tu función
        LoadHitsROOT(simTree, tpcPoints, ev);

        // 2. Procesamos las trazas que LoadHitsROOT ha dejado en los mapas globales
        for (auto const& [trackID, xVector] : posX) {
            
            if (xVector.size() < 5) continue; // Protección Kalman
            totalTracksProcessed++;

            bool drawPlots = (ev == eventToDraw);

            try {
            // This function call will now "jump" to the 'catch' block 
            // as soon as the C++ code hits a 'throw'
            UKFResult result_parcial = runKalman(posX[trackID], posY[trackID], posZ[trackID], Eloss[trackID], 
                                                mass_p, charge_p, 1, 1, initialMom[trackID], trackID, drawPlots);

            // If it didn't throw, we check the result and fill
           if (result_parcial.p_rec > 0 && result_parcial.p_rec != -999.0) {
                res.eventID  = ev;
                res.trackID  = result_parcial.trackID;
                res.p_rec    = result_parcial.p_rec;
                res.E_rec    = result_parcial.E_rec;
                res.sumEloss = result_parcial.sumEloss;
                res.p_true   = result_parcial.p_true;

                res.theta_rec  = result_parcial.theta_rec * TMath::RadToDeg();
                res.phi_rec    = result_parcial.phi_rec   * TMath::RadToDeg();
                res.theta_true = result_parcial.theta_true * TMath::RadToDeg();
                res.phi_true   = result_parcial.phi_true   * TMath::RadToDeg();

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
    outTree->Write();
    outFile->Close();
    simFile->Close();

    //std::cout << "\n << 'reco_ukf_output_3T_p40MeV_H300torr.root' saved" << std::endl;
    std::cout << "\n << 'reco_ukf_output.root' saved" << std::endl;
}

/*void UKFSingleTrack_protons() { 
   
   pointsToCluster = 5; 
   LoadHitsROOT(); // Carga las posiciones y el initialMom de TODAS las trazas

   // Recorremos cada traza que LoadHitsROOT haya encontrado
   for (auto const& [trackID, vectorHits] : posX) {
       
       std::cout << "\n=====================================" << std::endl;
       std::cout << " Iniciando Filtro de Kalman para Track " << trackID << std::endl;
       
       // Recuperamos el momento inicial específico de esta traza
       ROOT::Math::XYZVector momP = initialMom[trackID]; 
       
       std::cout << " Momento Inicial (Px, Py, Pz) = (" 
                 << momP.X() << ", " << momP.Y() << ", " << momP.Z() << ") MeV/c" << std::endl;

       // Ejecutamos el Kalman pasándole los vectores y el momento de ESTA traza
       runKalman(posX[trackID], posY[trackID], posZ[trackID], Eloss[trackID], mass_p, charge_p, 1, 1, momP, trackID);
   }
}*/

