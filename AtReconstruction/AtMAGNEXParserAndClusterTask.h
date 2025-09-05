#ifndef ATMAGNEXPARSERANDCLUSTERTASK_H
#define ATMAGNEXPARSERANDCLUSTERTASK_H

#include <FairTask.h>

#include <Rtypes.h>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>

#include <memory>
#include <map>
#include <fstream>

#include "AtMAGNEXMap.h"

class TBuffer;
class TClass;
class TMemberInspector;

class AtMAGNEXParserAndClusterTask : public FairTask {
protected:
   Bool_t fIsPersistence{false};
   TString fInputFileName;
   TString fSiCFileName;
   TString fOutputBranchName;

   TString fChannelToColMapFileName;
   std::unique_ptr<std::ifstream> fChannelToColMapFile;
   std::unique_ptr<std::map<Int_t, Int_t>> fChannelToColMap;

   TString fDetectorParametersFileName;
   AtMAGNEXMap *fMap;
   Double_t f_ps_to_us = 1.0e-06;  // ps to us factor.
   Double_t fVDrift = 9.25;        // drift velocity in cm/us.
   Double_t fDetHeight = 185.;     // Detector height in mm.

   TClonesArray fEventArray;
   Int_t fEventNum{0};
   ULong64_t fEntryNum{0};
   ULong64_t fSiCEntryNum{0};
   ULong64_t fWindowSize{2000000};

   Int_t fStripCluster{1};
   ULong64_t fTimeCluster{100000};

   std::unique_ptr<TFile> fInputFile{nullptr};
   TTree *fInputTree{nullptr};
   UShort_t Board;
   UShort_t Channel;
   UShort_t FTS;
   ULong64_t CTS;
   ULong64_t Timestamp;
   UShort_t Charge;
   UInt_t Flags;
   UShort_t Pad;
   Double_t Charge_cal;
   UShort_t Row;
   UShort_t Section;

   std::unique_ptr<TFile> fSiCFile{nullptr};
   TTree *fSiCTree{nullptr};
   UShort_t Board_SiC;
   UShort_t Channel_SiC;
   UShort_t FTS_SiC;
   ULong64_t CTS_SiC;
   ULong64_t Timestamp_SiC;
   UShort_t Charge_SiC;
   UInt_t Flags_SiC;
   Double_t Charge_cal_SiC;

public:
   AtMAGNEXParserAndClusterTask(TString inputFileName, TString SiCFileName, TString channelToColMapFileName, TString detParFileName, TString outputBranchName = "AtHitClusterEventH");

   void SetPersistence(bool value) { fIsPersistence = value; }
   void SetWindowSize(ULong64_t value) { fWindowSize = value; }
   void SetDriftVelocity(Double_t value) { fVDrift = value; }
   void SetStripCluster(Int_t value) { fStripCluster = value; }
   void SetTimeCluster(ULong64_t value) { fTimeCluster = value; }

   virtual InitStatus Init() override;
   virtual void Exec(Option_t *opt) override;

   ClassDefOverride(AtMAGNEXParserAndClusterTask, 1);
};

#endif
