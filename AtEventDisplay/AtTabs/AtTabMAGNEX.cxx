#include "AtTabMAGNEX.h"

#include "AtTabMain.h"
#include "AtPatternEvent.h"
#include "AtHitClusterEvent.h"
#include "AtViewerManager.h"
#include "AtMap.h"

#include <FairLogger.h>

#include <TEveManager.h>
#include <TEveEventManager.h>
#include <TEvePointSet.h>
#include <Rtypes.h>
#include <TStyle.h>
#include <TH2Poly.h>
#include <TCanvas.h>

namespace DataHandling {
class AtSubject;
}

ClassImp(AtTabMAGNEX);

AtTabMAGNEX::AtTabMAGNEX() : AtTabMain()
{
   fHitClusterEventBranch = &AtViewerManager::Instance()->GetHitClusterEventBranch();
   fHitClusterEventBranch->Attach(this);
}


AtTabMAGNEX::~AtTabMAGNEX()
{
   fHitClusterEventBranch->Detach(this);
   AtTabMain::~AtTabMain();
}

void AtTabMAGNEX::InitTab()
{
   std::cout << " =====  AtTabMAGNEX::Init =====" << std::endl;

   gEve->AddEvent(fEveHitClusterEvent.get());
   gEve->AddEvent(fEveMeanHits.get());
   fEveMeanHits->AddElement(fMeanHitSet.get());
   gEve->AddEvent(fEvePatternEvent.get());

   fTabInfo->AddAugment(std::make_unique<AtTabInfoFairRoot<AtPatternEvent>>(*fPatternEventBranch));
   fTabInfo->AddAugment(std::make_unique<AtTabInfoFairRoot<AtHitClusterEvent>>(*fHitClusterEventBranch));

   gStyle->SetPalette(55);

   std::cout << " AtTabMAGNEX::Init : Initialization complete! "
             << "\n";
}

void AtTabMAGNEX::Update(DataHandling::AtSubject *sub)
{
   if (sub == fHitClusterEventBranch || sub == fEntry) {
      UpdatePadPlane();
      UpdateHitClusterEventElements();
   }
   if (sub == fPatternEventBranch || sub == fEntry) {
      AtTabMain::UpdatePatternEventElements();
   }

   if (sub == fHitClusterEventBranch || sub == fPatternEventBranch || sub == fEntry) {
      gEve->Redraw3D(false);
   }
}

void AtTabMAGNEX::UpdateRenderState()
{
   fEveHitClusterEvent->SetRnrState(true);
   fEveMeanHits->SetRnrState(false);
   fEvePatternEvent->SetRnrState(false);
}

void AtTabMAGNEX::UpdateHitClusterEventElements()
{
   if (fEveHitClusterEvent == nullptr)
      return;

   auto fHitClusterEvent = GetFairRootInfo<AtHitClusterEvent>();
   if (fHitClusterEvent == nullptr) {
      LOG(debug) << "Cannot update AtHitClusterEvent elements: no event available";
      return;
   }

   auto &hitClusters = fHitClusterEvent->GetHitClusters();
   ExpandNumHitClusters(hitClusters.size());

   fEveHitClusterEvent->RemoveElements();

   fMeanHitSet->Reset(hitClusters.size());
   fMeanHitSet->SetOwnIds(true);
   fMeanHitSet->SetMarkerColorAlpha(kPink, 0.2);
   fMeanHitSet->SetMarkerSize(2);
   fMeanHitSet->SetMarkerStyle(kFullCircle);

   for (Int_t i = 0; i < hitClusters.size(); ++i) {
      auto hitSet = fHitClusterSets.at(i).get();
      fHitAttr.Copy(*hitSet);
      hitSet->SetMarkerColor(GetTrackColor(i));
      SetPointsFromHitCluster(*hitSet, *hitClusters[i]);
      fEveHitClusterEvent->AddElement(hitSet);

      auto position = hitClusters[i]->GetPositionCharge();
      fMeanHitSet->SetNextPoint(position.X() / 10., position.Y() / 10., position.Z() / 10.);
      fMeanHitSet->SetPointId(new TNamed(Form("Mean Hit %d", i), ""));
   }

   gEve->ElementChanged(fMeanHitSet.get());
}

void AtTabMAGNEX::ExpandNumHitClusters(Int_t num)
{
   if (fHitClusterSets.size() < num)
      LOG(info) << "Expanding number of hit clusters to " << num << " from " << fHitClusterSets.size() << "." << std::endl;
   while (fHitClusterSets.size() < num) {
      Int_t clusterID = fHitClusterSets.size();

      auto hitClusterSet = std::make_unique<TEvePointSet>(TString::Format("HitCluster_%d", clusterID));
      hitClusterSet->SetDestroyOnZeroRefCnt(false);
      fHitAttr.Copy(*hitClusterSet);
      hitClusterSet->SetMarkerColor(GetTrackColor(clusterID));

      fHitClusterSets.push_back(std::move(hitClusterSet));
   }
}

void AtTabMAGNEX::SetPointsFromHitCluster(TEvePointSet &hitSet, const AtHitClusterFull &hitCluster)
{
   Int_t nHits = hitCluster.GetClusterSize();

   hitSet.Reset(nHits);
   hitSet.SetOwnIds(true);

   for (Int_t i = 0; i < nHits; ++i) {
      auto hit = hitCluster.GetHits()[i];

      if (hit.GetCharge() < fThreshold)
         continue;

      auto position = hit.GetPosition();
      hitSet.SetNextPoint(position.X() / 10., position.Y() / 10., position.Z() / 10.);
      hitSet.SetPointId(new TNamed(Form("Hit %d", i), ""));
   }

   gEve->ElementChanged(&hitSet);
}

void AtTabMAGNEX::UpdatePadPlane()
{
   if (fPadPlane)
      fPadPlane->Reset(nullptr);
   else
      return;

   LOG(debug) << "Updating pad plane.";

   auto fHitClusterEvent = GetFairRootInfo<AtHitClusterEvent>();
   if (fHitClusterEvent == nullptr) {
      LOG(debug) << "Could not fill pad plane histogram: no event available.";
      return;
   }
   auto &hitClusters = fHitClusterEvent->GetHitClusters();

   TString plane = AtViewerManager::Instance()->GetMap()->GetPadPlanePlane();

   for (auto &hitCluster : hitClusters) {
      Int_t nHits = (*hitCluster).GetClusterSize();
      for (Int_t i = 0; i < nHits; ++i) {
         auto hit = (*hitCluster).GetHits()[i];
         auto position = hit.GetPosition();

         if (plane == "XY") {
            fPadPlane->Fill(position.X(), position.Y(), hit.GetCharge());
         } else if (plane == "XZ") {
            fPadPlane->Fill(position.X(), position.Z(), hit.GetCharge());
         } else if (plane == "YZ") {
            fPadPlane->Fill(position.Y(), position.Z(), hit.GetCharge());
         }
      }
   }

   fCvsPadPlane->Modified();
   fCvsPadPlane->Update();
}
