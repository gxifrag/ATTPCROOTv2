#ifndef SPECTRALMODEL_H
#define SPECTRALMODEL_H

#include "TF1.h"
#include "TF1Convolution.h"
#include "TGraph.h"

#include <functional>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

// Full normalized BW-with-penetrability lineshape, evaluated at x, given the
// resonance energy ER and the width AT resonance Gamma0 (i.e. Gamma(ER) = Gamma0).
// Penetrability enters as a multiplicative ratio P(x)/P(ER) inside Gamma(x),
// NOT as a convolution -- the convolution with the instrumental resolution
// happens afterwards, on the outside, via TF1Convolution.
using GammaFunc = std::function<double(double x, double ER, double Gamma0)>;

struct GaussComponent {
   std::string label; // e.g. "g0", "g1"
   double ampInit, meanInit, sigmaInit;
};

struct BWComponent {
   std::string label; // e.g. "bw0", "bw1", "bw2"
   int L;              // bookkeeping/printing only
   GammaFunc gammaFunc;
   double ampInit, erInit, sigmaInit, gamma0Init;
};

class SpectralModel {
public:
   explicit SpectralModel(int nConvFFTPoints = 1000) : fNConvPoints(nConvFFTPoints) {}

   void AddGaussian(const std::string &label, double ampInit, double meanInit, double sigmaInit);

   // Analytic penetrability (R-matrix, l = 0,1,2), same formula your colleague uses.
   // s = separation energy [MeV], mu = reduced mass [MeV/c^2], R = channel radius [fm].
   void AddBW(const std::string &label, int L, double s, double mu, double R, double ampInit, double erInit,
              double sigmaInit, double gamma0Init);

   // Fallback: your tabulated penetrability TGraph (P vs En = E - Sn), used as
   // the ratio P(E)/P(ER) so Gamma0 keeps the same "width at resonance" meaning
   // as in AddBW above (no gamma^2 conversion needed).
   void AddBWTabulated(const std::string &label, int L, TGraph *penGraph, double s, double ampInit, double erInit,
                       double sigmaInit, double gamma0Init);

   void SetPhaseSpace(TGraph *psGraph, double ampInit);

   // Call once, after all components are registered.
   void Finalize();

   int NPar() const { return static_cast<int>(fParNames.size()); }
   const std::string &ParName(int i) const { return fParNames.at(i); }
   int Idx(const std::string &label) const; // e.g. Idx("bw1_ER"), Idx("g0_Sigma")

   // TF1-compatible functor.
   double operator()(double *x, double *p);

   // Build a ready-to-fit TF1 with named parameters and initial values set.
   // Stores xMin/xMax internally -- needed for the TF1Convolution objects.
   TF1 *Build(const char *name, double xMin, double xMax);

   // After fitting: standalone TF1 for ONE component with best-fit params
   // plugged in (for plotting individual components).
   TF1 *GetComponentTF1(const std::string &label, const TF1 *fitted, double xMin, double xMax);

   std::vector<std::string> GaussLabels() const;
   std::vector<std::string> BWLabels() const;
   bool HasPhaseSpace() const { return fPS != nullptr; }

   TGraph* GetPS() const { return fPS; }

private:
   int fNConvPoints;
   double fXMin = 0, fXMax = 0;

   std::vector<GaussComponent> fGauss;
   std::vector<BWComponent> fBW;
   TGraph *fPS = nullptr;
   double fPSAmpInit = 0;

   std::vector<std::string> fParNames;
   std::map<std::string, int> fLabelToFirstIdx;

   // Cached TF1 / TF1Convolution per BW component label. Created once,
   // parameters updated in place -- this is what avoids the per-call
   // reconstruction cost you hit before with TSpline3.
   std::map<std::string, std::shared_ptr<TF1>> fBWFuncs;
   std::map<std::string, std::shared_ptr<TF1>> fGaussFuncs;
   std::map<std::string, std::shared_ptr<TF1Convolution>> fConvObjs;

   std::vector<double> fLastPars; // skip recomputing convolutions if p hasn't changed

   const BWComponent *FindBW(const std::string &label) const;
   void EnsureConvolution(const std::string &label, double ER, double sigma, double Gamma0);
   bool ParametersChanged(const double *p, int n);
};

#endif