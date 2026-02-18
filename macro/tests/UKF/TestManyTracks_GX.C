std::string getEnergyPath()
{
   auto env = std::getenv("VMCWORKDIR");
   if (env == nullptr) {
      return "../../resources/energy_loss/HinH_better.txt"; // Default path assuming cwd is build/AtTools
   }
   return std::string(env) + "/resources/energy_loss/HinH_better.txt"; // Use environment variable
}
using ROOT::Math::XYZPoint;
using ROOT::Math::XYZVector;

const double mass_p = 938.272;           // Mass of proton in MeV/c^2
const double charge_p = 1.602176634e-19; // Charge of proton

// Vectors to store the simulated points to compare to
std::vector<double> x_sim, y_sim, z_sim, Eloss_sim;
TH1F *hMom = nullptr;
TH1F *hMomSampled = nullptr;
TH1F *hMomError = nullptr;
TCanvas *c1 = new TCanvas();
TCanvas *c2 = new TCanvas();

// Parameters for model
const int pointsToCluster = 5;
const double sigma_pos = 1;                // Position uncertainty of 10 mm
const double sigma_mom = 0.1;              // Momentum uncertainty in percentage
const double sigma_mom_sample = 0.05;      // Sampled momentum uncertainty in percentage
const double sigma_theta = 1 * M_PI / 180; // Angular uncertainty of 1 degree
const double sigma_phi = 1 * M_PI / 180;   // Angular uncertainty of 1 degree
const double gasDensity = 3.3084e-5;        // g/cm^3
const double fAlpha = 1e-1;
const double fBeta = 2;
const double fKappa = 0;

// Global variables
kf::TrackFitterUKF *ukf = nullptr;

// Functions in file
/// Loads hits from an input file
// @brief Loads hits from an input file for a specific event and track
void LoadHits(int targetEvent, int targetTrack);
/// Run a single UKF for the passed initial state
void SingleUKF(XYZPoint initialPos, XYZVector initialMom, TMatrixD initialCov);
/// @brief Create the ukf pointer to reuse.
void CreateUKF();

TMatrixD CalculateInitialCov(double p);
TMatrixD CalculatePosCov();

/// @brief  Test many tracks, saving stat properties
/// @param n
/// @param bias
void TestManyTracks_GX(int n, double bias = 0)
{
   LoadHits(1,1); // Load hits for event 1, track 1 (proton=menos energy loss)
   CreateUKF();

   XYZPoint fTruePos(0.00140, 0.00232, 11.7002);
   XYZVector fTrueMom(-0.0889103, -0.101446, 0.132321); // In GeV/c
   fTrueMom *= 1e3; //In MeV/c like the code expects


   double fSigmaMom = fTrueMom.R() * sigma_mom_sample;

   hMom = new TH1F("hMom", "Reconstructed Momentum (MeV/c)", 100, fTrueMom.R() - 4 * fSigmaMom,
                   fTrueMom.R() + 4 * fSigmaMom);
   hMomSampled = new TH1F("hMom", "Reconstructed Momentum (MeV/c)", 100, fTrueMom.R() - 4 * fSigmaMom,
                          fTrueMom.R() + 4 * fSigmaMom);
   hMomError =
      new TH1F("hMomError", "Error (%)", 100, -4 * fSigmaMom / fTrueMom.R() * 100, 4 * fSigmaMom / fTrueMom.R() * 100);

   for (int i = 0; i < n; ++i) {

      if (i % 100 == 0)
         std::cout << "On iteration " << i << std::endl;

      double pSampled = gRandom->Gaus(fTrueMom.R(), sigma_mom_sample * fTrueMom.R());
      pSampled += bias;
      hMomSampled->Fill(pSampled);

      // pSampled = fTrueMom.R();

      ROOT::Math::Polar3DVector sampledMom(pSampled, fTrueMom.Theta(), fTrueMom.Phi());
      SingleUKF(fTruePos, XYZVector(sampledMom), CalculateInitialCov(pSampled));

      auto filtState = ukf->GetFilteredStates()[0];
      auto smoothState = ukf->GetSmoothedStates()[0];

      double pReco = smoothState[3];
      hMom->Fill(pReco);
      double error = (pReco - fTrueMom.R()) / fTrueMom.R() * 100;
      hMomError->Fill(error);

      std::cout << std::endl
                << std::endl
                << "With initial momentum " << pSampled << " reconstructed " << pReco << " Error: " << error << "%"
                << std::endl
                << std::endl;
   }

   // Draw results
   c1->cd();
   hMom->Draw("hist");
   hMomSampled->SetLineColor(kRed);
   hMomSampled->Draw("same hist");

   auto legend = new TLegend(0.65, 0.75, 0.88, 0.88);
   legend->AddEntry(hMom, "Reconstructed", "l");
   legend->AddEntry(hMomSampled, "Sampled", "l");
   legend->Draw();
   c2->cd();
   hMomError->Draw("hist");
}

void LoadHits(int targetEvent, int targetTrack)
{
   x_sim.clear(); y_sim.clear(); z_sim.clear(); Eloss_sim.clear();
   std::ifstream infile("hits_attpcsim_with_ids_Bfield.txt"); // Asegúrate de que el nombre sea correcto
   
   int ev, trk;
   double xi, yi, zi, Ei;
   double eLossAccum = 0;
   int hitCount = 0;

   while (infile >> ev >> trk >> xi >> yi >> zi >> Ei) {
      // Filtramos por el evento y la traza que nos interesa
      if (ev == targetEvent && trk == targetTrack) {
         
         // Agrupamos puntos según 'pointsToCluster' (ej. cada 5 hits)
         if (hitCount % pointsToCluster == 0) {
            x_sim.push_back(xi * 10.0); // Convertir cm a mm
            y_sim.push_back(yi * 10.0);
            z_sim.push_back(zi * 10.0);
            if (hitCount > 0) Eloss_sim.push_back(eLossAccum);
            eLossAccum = 0;
         }
         
         eLossAccum += Ei;
         hitCount++;
      }
   }

   if (x_sim.size() == 0) {
      std::cerr << "ERROR: No se encontraron hits para Ev=" << targetEvent << " Trk=" << targetTrack << std::endl;
   }
   std::cout << "Cargados " << x_sim.size() << " clusters para el proton." << std::endl;
}

void CreateUKF()
{
   auto elossModel2 = std::make_unique<AtTools::AtELossCATIMA>(gasDensity);
   elossModel2->SetProjectile(1, 1, 1);
   std::vector<std::tuple<int, int, int>> mat;
   mat.push_back({1, 1, 0});
   elossModel2->SetMaterial(mat);

   auto elossModel = std::make_unique<AtTools::AtELossTable>(0);
   elossModel->LoadSrimTable(getEnergyPath()); // Use the function to get the path
   elossModel->SetDensity(gasDensity);

   /*AtTools::AtPropagator propagator(charge_p, mass_p, std::move(elossModel));
   propagator.SetEField({0, 0, 0});    // No electric field
   //propagator.SetBField({0, 0, 2.85}); // Magnetic field
*/

   AtTools::AtPropagator propagator(charge_p, mass_p, std::move(elossModel2));
   propagator.SetEField({0, 0, 0});    // No electric field
   propagator.SetBField({0, 0, 2.85}); // Magnetic field
   //propagator.SetBField({0, 0, 0}); // No magnetic field
   // Setup stepper for UKF
   auto stepper = std::make_unique<AtTools::AtRK4Stepper>();

   // Setup UKF
   ukf = new kf::TrackFitterUKF(std::move(propagator), std::move(stepper));
   ukf->setParameters(fAlpha, fBeta, fKappa);
}

TMatrixD CalculateInitialCov(double p)
{
   TMatrixD cov(6, 6);
   cov.Zero();
   for (int i = 0; i < 3; ++i) {

      cov(i, i) = sigma_pos * sigma_pos; // Set diagonal covariance to some small number
   }
   cov(3, 3) = p * p * sigma_mom * sigma_mom; // Momentum uncertainty
   cov(4, 4) = sigma_theta * sigma_theta;     // Angular uncertainty
   cov(5, 5) = sigma_phi * sigma_phi;         // Angular uncertainty

   return cov;
}
TMatrixD CalculatePosCov()
{
   TMatrixD cov(3, 3);
   cov.Zero();
   for (int i = 0; i < 3; ++i) {

      cov(i, i) = sigma_pos * sigma_pos; // Set diagonal covariance to some small number
   }
   return cov;
}

void SingleUKF(XYZPoint intitialPos, XYZVector initialMom, TMatrixD intitialCov)
{
   ukf->SetInitialState(intitialPos, initialMom, intitialCov);

   for (int i = 1; i < x_sim.size(); ++i) {
      ukf->SetMeasCov(CalculatePosCov());
      XYZPoint point(x_sim[i], y_sim[i], z_sim[i]);

      ukf->predictUKF(point);
      ukf->correctUKF(point);
   }

   ukf->smoothUKF();
}
