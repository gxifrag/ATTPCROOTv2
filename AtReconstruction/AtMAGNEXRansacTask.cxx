#include "AtMAGNEXRansacTask.h"

#include "AtEstimatorMethods.h"
#include "AtHitClusterEvent.h" // for AtHitClusterEvent
#include "AtPatternEvent.h"
#include "AtPatternTypes.h"
#include "AtSampleConsensus.h"
#include "AtSampleMethods.h"

#include <FairLogger.h>      // for LOG, Logger
#include <FairRootManager.h> // for FairRootManager

#include <TClonesArray.h> // for TClonesArray
#include <TObject.h>      // for TObject

#include <memory> // for allocator

ClassImp(AtMAGNEXRansacTask);

AtMAGNEXRansacTask::AtMAGNEXRansacTask()
   : fInputBranchName("AtHitClusterEventH"), fOutputBranchName("AtPatternEvent"), kIsPersistence(kFALSE), kIsReprocess(kFALSE),
     fPatternEventArray("AtPatternEvent", 1)
{
}

AtMAGNEXRansacTask::~AtMAGNEXRansacTask() = default;

void AtMAGNEXRansacTask::SetPersistence(Bool_t value)
{
   kIsPersistence = value;
}
void AtMAGNEXRansacTask::SetInputBranch(TString branchName)
{
   fInputBranchName = branchName;
}

void AtMAGNEXRansacTask::SetOutputBranch(TString branchName)
{
   fOutputBranchName = branchName;
}

void AtMAGNEXRansacTask::SetModelType(int model)
{
   fRANSACModel = model;
}
void AtMAGNEXRansacTask::SetDistanceThreshold(Float_t threshold)
{
   fRANSACThreshold = threshold;
}
void AtMAGNEXRansacTask::SetMinHitsLine(Int_t nhits)
{
   fMinHitsLine = nhits;
}
void AtMAGNEXRansacTask::SetNumItera(Int_t niterations)
{
   fNumItera = niterations;
}
void AtMAGNEXRansacTask::SetAlgorithm(Int_t val)
{
   fRANSACAlg = val;
}
void AtMAGNEXRansacTask::SetRanSamMode(Int_t mode)
{
   fRandSamplMode = mode;
};
void AtMAGNEXRansacTask::SetIsReprocess(Bool_t value)
{
   kIsReprocess = value;
}
void AtMAGNEXRansacTask::SetInputBranchName(TString inputName)
{
   fInputBranchName = inputName;
}
void AtMAGNEXRansacTask::SetOutputBranchName(TString outputName)
{
   fOutputBranchName = outputName;
}

InitStatus AtMAGNEXRansacTask::Init()
{

   FairRootManager *ioMan = FairRootManager::Instance();
   if (ioMan == nullptr) {
      LOG(error) << "Cannot find RootManager!";
      return kERROR;
   }

   fEventArray = dynamic_cast<TClonesArray *>(ioMan->GetObject(fInputBranchName));
   if (fEventArray == nullptr) {

      LOG(error) << "Cannot find AtHitClusterEvent array!";
      return kERROR;
   }

   ioMan->Register(fOutputBranchName, "MAGNEX", &fPatternEventArray, kIsPersistence);

   if (kIsReprocess) {

      ioMan->Register("AtHitClusterEventH", "MAGNEX", fEventArray, kIsPersistence);
   }

   return kSUCCESS;
}

void AtMAGNEXRansacTask::Exec(Option_t *opt)
{

   if (fEventArray->GetEntriesFast() == 0)
      return;

   fHitClusterEvent = dynamic_cast<AtHitClusterEvent *>(fEventArray->At(0));
   if (fHitClusterEvent->GetEventID() == fPreviousEventNum)
      return;
   fPreviousEventNum = fHitClusterEvent->GetEventID();

   LOG(debug) << "Running (MAGNEX) RANSAC with " << fHitClusterEvent->GetNumHits() << " hits in " << fHitClusterEvent->GetNumHitClusters() << " clusters.";
   LOG(debug) << "Running Unified RANSAC";

   auto sampleMethod = static_cast<RandomSample::SampleMethod>(fRandSamplMode);
   auto patternType = AtPatterns::PatternType::kLine;
   auto estimator = SampleConsensus::Estimators::kChi2;
   if (fRANSACAlg == 2)
      estimator = SampleConsensus::Estimators::kMLESAC;
   if (fRANSACAlg == 3)
      estimator = SampleConsensus::Estimators::kLMedS;
   if (fRANSACAlg == 4)
      estimator = SampleConsensus::Estimators::kWRANSAC;
   SampleConsensus::AtSampleConsensus ransac(estimator, patternType, sampleMethod);

   ransac.SetDistanceThreshold(fRANSACThreshold);
   ransac.SetMinHitsPattern(fMinHitsLine);
   ransac.SetNumIterations(fNumItera);
   ransac.SetChargeThreshold(fChargeThres);
   fPatternEventArray.Delete();
   auto patternEvent = ransac.Solve(fHitClusterEvent);
   new (fPatternEventArray[0]) AtPatternEvent(patternEvent);
}
