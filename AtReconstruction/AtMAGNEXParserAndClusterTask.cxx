#include "AtMAGNEXParserAndClusterTask.h"

#include <FairTask.h>
#include <FairLogger.h>
#include <FairRootManager.h>

#include <Rtypes.h>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <Math/Point2D.h>
#include <Math/Point3D.h>

#include <memory>
#include <map>
#include <fstream>
#include <utility>
#include <string>
#include <sstream>
#include <vector>
#include <iterator>

#include "AtHitClusterEvent.h"
#include "AtHitClusterFull.h"

using XYPoint = ROOT::Math::XYPoint;
using XYZPoint = ROOT::Math::XYZPoint;

AtMAGNEXParserAndClusterTask::AtMAGNEXParserAndClusterTask(TString inputFileName, TString SiCFileName, TString channelToColMapFileName, TString detParFileName, TString outputBranchName)
   : fInputFileName(std::move(inputFileName)), fSiCFileName(std::move(SiCFileName)), fOutputBranchName(std::move(outputBranchName)),
     fChannelToColMapFileName(std::move(channelToColMapFileName)), fDetectorParametersFileName(std::move(detParFileName)),
     fEventArray("AtHitClusterEvent", 1)
{
}

InitStatus AtMAGNEXParserAndClusterTask::Init()
{
   // Find RootManager and register the output branch.
   FairRootManager *ioMan = FairRootManager::Instance();
   if (ioMan == nullptr) {
      LOG(error) << "AtMAGNEXParserAndClusterTask::Init Error : Could not find RootManager!";
      return kERROR;
   }
   ioMan->Register(fOutputBranchName, "MAGNEX", &fEventArray, fIsPersistence);

   // Open the data file containing the TPC related data and set up its TTree.
   fInputFile = std::make_unique<TFile>(fInputFileName, "READ");
   if (fInputFile->IsZombie()) {
      LOG(error) << "AtMAGNEXParserAndClusterTask::Init Error : Could not open ROOT file " << fInputFileName << "!";
      return kERROR;
   } else {
      LOG(info) << "ROOT file " << fInputFileName << " opened successfully.";
   }

   fInputTree = (TTree *)fInputFile->Get("Data_R");
   LOG(info) << "Total number of entries in the tree: " << fInputTree->GetEntries() << ".";
   fInputTree->SetBranchAddress("Board", &Board);
   fInputTree->SetBranchAddress("Channel", &Channel);
   fInputTree->SetBranchAddress("FineTSInt", &FTS);
   fInputTree->SetBranchAddress("CoarseTSInt", &CTS);
   fInputTree->SetBranchAddress("Timestamp", &Timestamp);
   fInputTree->SetBranchAddress("Charge", &Charge);
   fInputTree->SetBranchAddress("Flags", &Flags);
   fInputTree->SetBranchAddress("Pads", &Pad);
   fInputTree->SetBranchAddress("Charge_cal", &Charge_cal);
   fInputTree->SetBranchAddress("Row", &Row);
   fInputTree->SetBranchAddress("Section", &Section);

   // Open the data file containing the SiC related data and set up its TTree.
   fSiCFile = std::make_unique<TFile>(fSiCFileName, "READ");
   if (fSiCFile->IsZombie()) {
      LOG(error) << "AtMAGNEXParserAndClusterTask::Init Error : Could not open ROOT file " << fSiCFileName << "!";
      return kERROR;
   } else {
      LOG(info) << "ROOT file " << fSiCFileName << " opened successfully.";
   }

   fSiCTree = (TTree *)fSiCFile->Get("Data_R");
   LOG(info) << "Total number of entries in the tree: " << fSiCTree->GetEntries() << ".";
   fSiCTree->SetBranchAddress("Board", &Board_SiC);
   fSiCTree->SetBranchAddress("Channel", &Channel_SiC);
   fSiCTree->SetBranchAddress("FineTSInt", &FTS_SiC);
   fSiCTree->SetBranchAddress("CoarseTSInt", &CTS_SiC);
   fSiCTree->SetBranchAddress("Timestamp", &Timestamp_SiC);
   fSiCTree->SetBranchAddress("Charge", &Charge_SiC);
   fSiCTree->SetBranchAddress("Flags", &Flags_SiC);
   fSiCTree->SetBranchAddress("Charge_cal", &Charge_cal_SiC);

   // Set the mapping from electronic channel to column (strip).
   fChannelToColMapFile = std::make_unique<std::ifstream>(fChannelToColMapFileName);
   fChannelToColMap = std::make_unique<std::map<Int_t, Int_t>>();
   if (!fChannelToColMapFile->is_open()) {
      LOG(error) << "AtMAGNEXParserAndClusterTask::Init Error : Could not open the channel to column (strip) file " << fChannelToColMapFileName << "!";
      return kERROR;
   } else {
      LOG(info) << "Channel to column (strip) file " << fChannelToColMapFileName << " opened successfully.";

      std::string line;
      for (Int_t i = 0; i < 8; ++i) {
         std::getline(*fChannelToColMapFile, line);
         LOG(debug) << line;
      }

      while (!fChannelToColMapFile->eof()) {
         Int_t channel, column;
         std::getline(*fChannelToColMapFile, line);
         std::istringstream buffer(line);
         buffer >> channel >> column;
         LOG(debug) << "Channel: " << channel << " Column: " << column << ".";
         fChannelToColMap->insert(std::make_pair(channel, column));
      }
   }

   // Set the AtMAGNEXMap.
   fMap = new AtMAGNEXMap(fDetectorParametersFileName);
   fMap->GeneratePadPlane();

   return kSUCCESS;
}

void AtMAGNEXParserAndClusterTask::Exec(Option_t *opt)
{
   while (true && fEntryNum != fInputTree->GetEntries()) {
      fInputTree->GetEntry(fEntryNum);
      Int_t Col = fChannelToColMap->find(Channel)->second;
      if (Col < fMap->GetColNum())
         break;
      fEntryNum++;
   }

   if (fEntryNum == fInputTree->GetEntries())
      return;

   Long64_t timeinit = Timestamp;
   Long64_t timeref = timeinit;

   auto *event = dynamic_cast<AtHitClusterEvent *>(fEventArray.ConstructedAt(0, "C"));
   event->SetEventID(fEventNum++);
   event->SetTimestamp(timeinit);

   std::vector<ULong64_t> SiC_entries;

   while (fSiCEntryNum < fSiCTree->GetEntries()) {
      fSiCTree->GetEntry(fSiCEntryNum);
      if (Timestamp_SiC + fWindowSize < timeinit) {
         ++fSiCEntryNum;
         continue;
      }
      if (Timestamp_SiC > timeinit)
         break;
      SiC_entries.push_back(fSiCEntryNum);
      timeref = Timestamp_SiC;
      ++fSiCEntryNum;
   }

   std::vector<std::vector<ULong64_t>> entriesByRow;
   for (Int_t i = 0; i < fMap->GetRowNum(); ++i) {
      std::vector<ULong64_t> entriesVector{};
      entriesByRow.push_back(entriesVector);
   }

   while (fEntryNum < fInputTree->GetEntries()) {
      fInputTree->GetEntry(fEntryNum);
      if (Timestamp - timeinit > fWindowSize)
         break;
      Int_t Col = fChannelToColMap->find(Channel)->second;
      if (Col < fMap->GetColNum())
         entriesByRow[Row].push_back(fEntryNum);
      ++fEntryNum;
   }

   Int_t hitNum{0}, clusterNum{0};
   for (Int_t i = 0; i < fMap->GetRowNum(); ++i) {
      while (entriesByRow[i].size()) {
         AtHitClusterFull *hitCluster = new AtHitClusterFull();
         hitCluster->SetClusterID(clusterNum++);

         ULong64_t clusterEntrySeed = entriesByRow[i].back();
         entriesByRow[i].pop_back();
         fInputTree->GetEntry(clusterEntrySeed);

         Int_t Col = fChannelToColMap->find(Channel)->second;

         XYPoint point = fMap->CalcPadCenter(fMap->PadID(Col, Row));
         Double_t x = point.X();
         Double_t y = fDetHeight - fVDrift * (Timestamp - timeref) * f_ps_to_us * 10.;
         Double_t z = point.Y();

         AtHit *firstHit = new AtHit(hitNum++, fMap->PadID(Col, Row), XYZPoint(x, y, z), Charge);
         firstHit->SetTimeStamp(Timestamp);

         hitCluster->AddHit(*firstHit);

         Int_t minCol = Col - fStripCluster;
         Int_t maxCol = Col + fStripCluster;

         ULong64_t minTS = Timestamp - fTimeCluster;
         ULong64_t maxTS = Timestamp + fTimeCluster;

         Int_t skipEntries = 0;
         while (skipEntries != entriesByRow[i].size()) {
            while (skipEntries != entriesByRow[i].size()) {
               Int_t index = entriesByRow[i].size() - 1 - skipEntries;

               ULong64_t entry = entriesByRow[i][index];
               fInputTree->GetEntry(entry);

               Col = fChannelToColMap->find(Channel)->second;
               if (minCol <= Col && Col <= maxCol && minTS < Timestamp && Timestamp < maxTS) {
                  point = fMap->CalcPadCenter(fMap->PadID(Col, Row));
                  x = point.X();
                  y = fDetHeight - fVDrift * (Timestamp - timeref) * f_ps_to_us * 10.;
                  z = point.Y();

                  AtHit *hit = new AtHit(hitNum++, fMap->PadID(Col, Row), XYZPoint(x, y, z), Charge);
                  hit->SetTimeStamp(Timestamp);

                  hitCluster->AddHit(*hit);

                  if (Col - fStripCluster < minCol)
                     minCol = Col - fStripCluster;
                  if (Col + fStripCluster > maxCol)
                     maxCol = Col + fStripCluster;

                  if (Timestamp - fTimeCluster < minTS)
                     minTS = Timestamp - fTimeCluster;
                  if (Timestamp + fTimeCluster > maxTS)
                     maxTS = Timestamp + fTimeCluster;

                  entriesByRow[i].erase(std::next(entriesByRow[i].begin(), index));
                  skipEntries = 0;
                  break;
               }

               skipEntries++;
            }
         }

         hitCluster->SetPosition(hitCluster->GetPositionCharge());
         event->AddHitCluster(hitCluster);
      }
   }

   LOG(info) << "Event " << fEventNum << " with timestamp " << timeinit << " has " << hitNum << " hits in " << clusterNum << " clusters.";
   for (auto SiCEntry : SiC_entries)
      LOG(info) << "Event " << fEventNum << " in coincidence with SiC entry " << SiCEntry << ".";
}

ClassImp(AtMAGNEXParserAndClusterTask)
