#include "SpectralModel.h"

#include "TMath.h"

#include <cmath>

static const double kHbar = 197.3269804; // MeV*fm
static const double kNorm = 0.159154943; // 1/(2*pi), normalized Lorentzian prefactor

// ---- Penetrability ratio builders: P(x)/P(ER), dimensionless, = 1 at x = ER ----

static GammaFunc BuildAnalyticGammaFunc(int l, double s, double mu, double R)
{
   std::function<double(double, double)> ratio;

   if (l == 0) {
      ratio = [s](double x, double er) { return std::pow((x - s) / (er - s), 0.5); };
   } else if (l == 1) {
      ratio = [s, mu, R](double x, double er) {
         double k2R2_x = (2. * mu * (x - s) * R * R) / (kHbar * kHbar);
         double k2R2_er = (2. * mu * (er - s) * R * R) / (kHbar * kHbar);
         return std::pow((x - s) / (er - s), 1.5) * (2. * (x - s) / ((er - s) + (x - s))) * (1. + k2R2_er) /
                (1. + k2R2_x);
      };
   } else if (l == 2) {
      ratio = [s, mu, R](double x, double er) {
         double k2R2_x = (2. * mu * (x - s) * R * R) / (kHbar * kHbar);
         double k2R2_er = (2. * mu * (er - s) * R * R) / (kHbar * kHbar);
         return std::pow((x - s) / (er - s), 2.5) * (2. * (x - s) / ((er - s) + (x - s))) *
                (9. + 6. * k2R2_er + k2R2_er * k2R2_er) / (9. + 6. * k2R2_x + k2R2_x * k2R2_x);
      };
   } else {
      throw std::runtime_error("BuildAnalyticGammaFunc(): only l = 0,1,2 implemented");
   }

   return [ratio, s](double x, double er, double Gamma0) {
      if (x <= s)
         return 0.0;
      double Gamma = Gamma0 * ratio(x, er);
      return Gamma * kNorm / ((x - er) * (x - er) + Gamma * Gamma / 4.);
   };
}

/*static GammaFunc BuildTabulatedGammaFunc(TGraph *gPen, double s)
{
   return [gPen, s](double x, double er, double Gamma0) {
      if (x <= s)
         return 0.0;
      double P_x = gPen->Eval(x - s);
      double P_er = gPen->Eval(er - s);
      if (P_er <= 0)
         return 0.0;
      double Gamma = Gamma0 * (P_x / P_er);
      return Gamma * kNorm / ((x - er) * (x - er) + Gamma * Gamma / 4.);
   };
}*/

// ---- Registration ----

void SpectralModel::AddGaussian(const std::string &label, double ampInit, double meanInit, double sigmaInit)
{
   fGauss.push_back({label, ampInit, meanInit, sigmaInit});
}

void SpectralModel::AddBW(const std::string &label, int L, double s, double mu, double R, double ampInit,
                          double erInit, double sigmaInit, double gamma0Init)
{
   BWComponent c;
   c.label = label;
   c.L = L;
   c.gammaFunc = BuildAnalyticGammaFunc(L, s, mu, R);
   c.ampInit = ampInit;
   c.erInit = erInit;
   c.sigmaInit = sigmaInit;
   c.gamma0Init = gamma0Init;
   fBW.push_back(std::move(c));
}

/*void SpectralModel::AddBWTabulated(const std::string &label, int L, TGraph *penGraph, double s, double ampInit,
                                   double erInit, double sigmaInit, double gamma0Init)
{
   if (!penGraph)
      throw std::runtime_error("SpectralModel::AddBWTabulated(): null TGraph for " + label);
   BWComponent c;
   c.label = label;
   c.L = L;
   c.gammaFunc = BuildTabulatedGammaFunc(penGraph, s);
   c.ampInit = ampInit;
   c.erInit = erInit;
   c.sigmaInit = sigmaInit;
   c.gamma0Init = gamma0Init;
   fBW.push_back(std::move(c));
}*/

void SpectralModel::SetPhaseSpace(TGraph *psGraph, double ampInit)
{
   fPS = psGraph;
   fPSAmpInit = ampInit;
}

void SpectralModel::Finalize()
{
   fParNames.clear();
   fLabelToFirstIdx.clear();

   for (auto &g : fGauss) {
      fLabelToFirstIdx[g.label] = static_cast<int>(fParNames.size());
      fParNames.push_back(g.label + "_Amp");
      fParNames.push_back(g.label + "_Mean");
      fParNames.push_back(g.label + "_Sigma");
   }
   for (auto &b : fBW) {
      fLabelToFirstIdx[b.label] = static_cast<int>(fParNames.size());
      fParNames.push_back(b.label + "_Amp");
      fParNames.push_back(b.label + "_ER");
      fParNames.push_back(b.label + "_Sigma");
      fParNames.push_back(b.label + "_Gamma0");
   }
   if (fPS) {
      fLabelToFirstIdx["ps0"] = static_cast<int>(fParNames.size());
      fParNames.push_back("ps0_Amp");
   }
}

int SpectralModel::Idx(const std::string &label) const
{
   auto pos = label.rfind('_');
   if (pos == std::string::npos)
      throw std::runtime_error("SpectralModel::Idx(): malformed label " + label);
   std::string comp = label.substr(0, pos);
   std::string par = label.substr(pos + 1);

   auto it = fLabelToFirstIdx.find(comp);
   if (it == fLabelToFirstIdx.end())
      throw std::runtime_error("SpectralModel::Idx(): unknown component " + comp);
   int base = it->second;

   static const std::map<std::string, int> gaussOffset = {{"Amp", 0}, {"Mean", 1}, {"Sigma", 2}};
   static const std::map<std::string, int> bwOffset = {{"Amp", 0}, {"ER", 1}, {"Sigma", 2}, {"Gamma0", 3}};

   for (auto &g : fGauss)
      if (g.label == comp)
         return base + gaussOffset.at(par);
   for (auto &b : fBW)
      if (b.label == comp)
         return base + bwOffset.at(par);
   if (comp == "ps0")
      return base;

   throw std::runtime_error("SpectralModel::Idx(): could not resolve " + label);
}

const BWComponent *SpectralModel::FindBW(const std::string &label) const
{
   for (auto &b : fBW)
      if (b.label == label)
         return &b;
   return nullptr;
}

// ---- Cached convolution machinery ----

/*void SpectralModel::EnsureConvolution(const std::string &label, double ER, double sigma, double Gamma0)
{
   const BWComponent *comp = FindBW(label);
   if (!comp)
      throw std::runtime_error("SpectralModel::EnsureConvolution(): unknown label " + label);

   if (!fBWFuncs.count(label)) {
      // 2 parameters for this inner TF1: [0] = ER, [1] = Gamma0
      auto bwLambda = [comp](double *x, double *p) { return comp->gammaFunc(x[0], p[0], p[1]); };
      fBWFuncs[label] = std::make_shared<TF1>(("fBW_" + label).c_str(), bwLambda, fXMin, fXMax, 2);
   }
   if (!fGaussFuncs.count(label)) {
      // 1 parameter: [0] = sigma
      fGaussFuncs[label] =
         std::make_shared<TF1>(("fGauss_" + label).c_str(), "TMath::Gaus(x,0,[0],true)", fXMin, fXMax);
   }
   if (!fConvObjs.count(label)) {
      fConvObjs[label] = std::make_shared<TF1Convolution>(fBWFuncs[label].get(), fGaussFuncs[label].get());
      fConvObjs[label]->SetRange(fXMin, fXMax);
      if (fNConvPoints > 0)
         fConvObjs[label]->SetNofPointsFFT(fNConvPoints);
   }
   // Parameter order inside TF1Convolution(f1, f2): all of f1's pars, then all of f2's pars.
   // f1 = BW(ER, Gamma0), f2 = Gauss(sigma) -> [ER, Gamma0, sigma]
   fConvObjs[label]->SetParameters(ER, Gamma0, sigma);
}*/


void SpectralModel::EnsureConvolution(const std::string &label, double ER, double sigma, double Gamma0)
{
   const BWComponent *comp = FindBW(label);
   if (!comp)
      throw std::runtime_error("SpectralModel::EnsureConvolution(): unknown label " + label);

   const double pad = 5.0; // MeV — generous margin, tune later if needed

   if (!fBWFuncs.count(label)) {
      auto bwLambda = [comp](double *x, double *p) { return comp->gammaFunc(x[0], p[0], p[1]); };
      fBWFuncs[label] = std::make_shared<TF1>(("fBW_" + label).c_str(), bwLambda,
                                                fXMin - pad, fXMax + pad, 2);
   }
   if (!fGaussFuncs.count(label)) {
      fGaussFuncs[label] = std::make_shared<TF1>(("fGauss_" + label).c_str(),
                                                   "TMath::Gaus(x,0,[0],true)",
                                                   fXMin - pad, fXMax + pad);
   }
   if (!fConvObjs.count(label)) {
      fConvObjs[label] = std::make_shared<TF1Convolution>(fBWFuncs[label].get(), fGaussFuncs[label].get());
      fConvObjs[label]->SetRange(fXMin, fXMax);   // <-- THIS stays the visible/output range, unchanged
      if (fNConvPoints > 0)
         fConvObjs[label]->SetNofPointsFFT(fNConvPoints);
   }
   fConvObjs[label]->SetParameters(ER, Gamma0, sigma);
}

bool SpectralModel::ParametersChanged(const double *p, int n)
{
   bool changed = (fLastPars.size() != static_cast<size_t>(n));
   if (!changed) {
      for (int i = 0; i < n; i++) {
         if (fLastPars[i] != p[i]) {
            changed = true;
            break;
         }
      }
   }
   if (changed)
      fLastPars.assign(p, p + n);
   return changed;
}

// ---- Evaluation ----

double SpectralModel::operator()(double *x, double *p)
{
   int n = NPar();

   if (ParametersChanged(p, n)) {
      int idx = static_cast<int>(fGauss.size()) * 3; // skip gaussian block, no conv objects there
      for (auto &b : fBW) {
         double er = p[idx + 1];
         double sigma = p[idx + 2];
         double gamma0 = p[idx + 3];
         EnsureConvolution(b.label, er, sigma, gamma0);
         idx += 4;
      }
   }

   double val = 0;
   int idx = 0;

   for (auto &g : fGauss) {
      double amp = p[idx++];
      double mean = p[idx++];
      double sigma = p[idx++];
      val += amp * TMath::Gaus(x[0], mean, sigma, false);
   }

   for (auto &b : fBW) {
      double amp = p[idx];
      idx += 4; // ER, Sigma, Gamma0 already consumed above when (re)building the convolution
      val += amp * (*fConvObjs.at(b.label))(x, nullptr);
   }

   /*if (fPS) {
      double psAmp = p[idx++];
      val += psAmp * fPS->Eval(x[0]);
   }*/

   if (fPS) {
      double psAmp = p[idx++];
      double x_val = x[0];
      
      // Determine the range of the graph to prevent linear extrapolation
      double x_first = fPS->GetX()[0];
      double x_last = fPS->GetX()[fPS->GetN() - 1];

      // Only evaluate if within the data range, otherwise return 0
      if (x_val < x_first || x_val > x_last) {
         val += 0.0;
      } else {
         val += psAmp * fPS->Eval(x_val);
      }
   }

   return val;
}

TF1 *SpectralModel::Build(const char *name, double xMin, double xMax)
{
   fXMin = xMin;
   fXMax = xMax;

   auto *f = new TF1(name, this, xMin, xMax, NPar(), "SpectralModel");
   for (int i = 0; i < NPar(); i++)
      f->SetParName(i, fParNames[i].c_str());

   int idx = 0;
   for (auto &g : fGauss) {
      f->SetParameter(idx++, g.ampInit);
      f->SetParameter(idx++, g.meanInit);
      f->SetParameter(idx++, g.sigmaInit);
   }
   for (auto &b : fBW) {
      f->SetParameter(idx++, b.ampInit);
      f->SetParameter(idx++, b.erInit);
      f->SetParameter(idx++, b.sigmaInit);
      f->SetParameter(idx++, b.gamma0Init);
   }
   if (fPS)
      f->SetParameter(idx++, fPSAmpInit);

   return f;
}

TF1 *SpectralModel::GetComponentTF1(const std::string &label, const TF1 *fitted, double xMin, double xMax)
{
   for (auto &g : fGauss) {
      if (g.label != label)
         continue;
      int base = Idx(label + "_Amp");
      auto *f = new TF1((label + "_fit").c_str(), "gaus(0)", xMin, xMax);
      f->SetParameters(fitted->GetParameter(base), fitted->GetParameter(base + 1), fitted->GetParameter(base + 2));
      return f;
   }

   for (auto &b : fBW) {
      if (b.label != label)
         continue;
      int base = Idx(label + "_Amp");
      double amp = fitted->GetParameter(base);
      double er = fitted->GetParameter(base + 1);
      double sigma = fitted->GetParameter(base + 2);
      double gamma0 = fitted->GetParameter(base + 3);

      EnsureConvolution(label, er, sigma, gamma0); // make sure the cached conv matches these exact values
      auto convObj = fConvObjs.at(label);

      auto *f = new TF1(
         (label + "_fit").c_str(), [convObj, amp](double *x, double *p) { return amp * (*convObj)(x, nullptr); },
         xMin, xMax, 0);
      return f;
   }

   if (fPS && label == "ps0") {
      int base = Idx("ps0_Amp");
      TF1 *f = new TF1(
         "ps0_fit", [=](double *x, double *p) { return p[0] * fPS->Eval(x[0]); }, xMin, xMax, 1);
      f->SetParameter(0, fitted->GetParameter(base));
      return f;
   }

   throw std::runtime_error("SpectralModel::GetComponentTF1(): unknown component label " + label);
}

std::vector<std::string> SpectralModel::GaussLabels() const
{
   std::vector<std::string> out;
   for (auto &g : fGauss)
      out.push_back(g.label);
   return out;
}

std::vector<std::string> SpectralModel::BWLabels() const
{
   std::vector<std::string> out;
   for (auto &b : fBW)
      out.push_back(b.label);
   return out;
}