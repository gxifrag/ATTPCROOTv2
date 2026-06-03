#include "AtEvent.h"
#include "AtHit.h"

#include "FairRun.h"
#include "FairRunAna.h"
#include "H5Cpp.h"
#include "H5File.h"
#include "TCanvas.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TH1F.h"
#include "TStopwatch.h"
#include "TString.h"
#include "TSystem.h"
#include "TTree.h"
#include "TTreePlayer.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include <H5Exception.h>
#include <hdf5.h>

#include <ios>
#include <iostream>
#include <istream>
#include <limits>
#include <map>
#include <vector>

// ROOT
#include "TApplication.h"
#include "TArrayD.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TGraph.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TRotation.h"
#include "TVectorD.h"

#include "Fit/Fitter.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/GenVector/AxisAngle.h"
#include "Math/GenVector/EulerAngles.h"
#include "Math/GenVector/Quaternion.h"
#include "Math/GenVector/Rotation3D.h"
#include "Math/GenVector/RotationX.h"
#include "Math/GenVector/RotationY.h"
#include "Math/GenVector/RotationZ.h"
#include "Math/GenVector/RotationZYX.h"
#include "Math/Minimizer.h"

hid_t _group;
hid_t _dataset;

enum class IO_MODE { READ, WRITE };

typedef struct ATHit_t {
   double x;
   double y;
   double z;
   int t;
   double A;
   int trackID;
   int pointIDMC;
   int trackIDMC;
   double energyMC;
   double elossMC;
   double angleMC;
   int AMC;
   int ZMC;

} ATHit_t;

const H5std_string MEMBER1("x");
const H5std_string MEMBER2("y");
const H5std_string MEMBER3("z");
const H5std_string MEMBER4("t");
const H5std_string MEMBER5("A");
const H5std_string MEMBER6("trackID");
const H5std_string MEMBER7("pointIDMC");
const H5std_string MEMBER8("energyMC");
const H5std_string MEMBER9("elossMC");
const H5std_string MEMBER10("angleMC");
const H5std_string MEMBER11("AMC");
const H5std_string MEMBER12("ZMC");

using namespace H5;

void usage();
