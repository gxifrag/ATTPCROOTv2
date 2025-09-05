#include "AtHitClusterEvent.h"

#include "AtHitClusterFull.h"

#include <FairLogger.h>

#include <Rtypes.h>


using HitClusterPtr = std::unique_ptr<AtHitClusterFull>;
using HitClusterVector = std::vector<HitClusterPtr>;

ClassImp(AtHitClusterEvent);

AtHitClusterEvent::AtHitClusterEvent() : AtBaseEvent("AtHitClusterEvent") {}

AtHitClusterEvent::AtHitClusterEvent(const AtHitClusterEvent &copy)
   : AtBaseEvent(copy), fEventCharge(copy.GetEventCharge())
{
   for (const auto &hitCluster : copy.GetHitClusters())
      fHitClusterArray.push_back(hitCluster->CloneCluster());
}

AtHitClusterEvent &AtHitClusterEvent::operator=(AtHitClusterEvent object)
{
   swap(*this, object);
   return *this;
}

void AtHitClusterEvent::Clear(Option_t *opt)
{
   AtBaseEvent::Clear(opt);
   fEventCharge = 0;
   fHitClusterArray.clear();
}

void AtHitClusterEvent::AddHitCluster(AtHitClusterFull *cluster)
{
   fHitClusterArray.push_back(cluster->CloneCluster());
   if (fHitClusterArray.back()->GetClusterID() == -1)
      fHitClusterArray.back()->SetClusterID(fHitClusterArray.size() - 1);
   updateEventCharge(cluster->GetCharge());

   LOG(debug) << "Added cluster with ID " << fHitClusterArray.back()->GetClusterID() << " and " << fHitClusterArray.back()->GetClusterSize() << " hits to event " << fEventID;
}

Int_t AtHitClusterEvent::GetNumHits() const
{
   Int_t hitNum{};
   for (const auto &hitCluster : fHitClusterArray)
      hitNum += hitCluster->GetClusterSize();
   return hitNum;
}
