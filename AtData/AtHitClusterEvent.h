#ifndef ATHITCLUSTEREVENT_H
#define ATHITCLUSTEREVENT_H

#include "AtBaseEvent.h"
#include "AtHitClusterFull.h"

#include <Rtypes.h>

#include <memory>
#include <vector>
#include <utility>

class TBuffer;
class TClass;
class TMemberInspector;

class AtHitClusterEvent : public AtBaseEvent {
public:
   using HitClusterPtr = std::unique_ptr<AtHitClusterFull>;
   using HitClusterVector = std::vector<HitClusterPtr>;

private:
   Double_t fEventCharge = 0;

   HitClusterVector fHitClusterArray;

   void updateEventCharge(Double_t charge) { fEventCharge += charge; }

public:
   AtHitClusterEvent();
   AtHitClusterEvent(const AtHitClusterEvent &copy);
   AtHitClusterEvent(const AtBaseEvent &copy) : AtBaseEvent(copy) { SetName("AtHitClusterEvent"); }
   AtHitClusterEvent &operator=(const AtHitClusterEvent object);
   virtual ~AtHitClusterEvent() = default;

   friend void swap(AtHitClusterEvent &first, AtHitClusterEvent &second)
   {
      using std::swap;
      swap(dynamic_cast<AtBaseEvent &>(first), dynamic_cast<AtBaseEvent &>(second));

      swap(first.fEventCharge, second.fEventCharge);
      swap(first.fHitClusterArray, second.fHitClusterArray);
   }

   void Clear(Option_t *opt = nullptr) override;

   void AddHitCluster(AtHitClusterFull *cluster);

   Int_t GetNumHits() const;
   Int_t GetNumHitClusters() const { return fHitClusterArray.size(); }
   Double_t GetEventCharge() const { return fEventCharge; }

   const AtHitClusterFull &GetHitCluster(Int_t hitNo) const { return *fHitClusterArray.at(hitNo); }
   const HitClusterVector &GetHitClusters() const { return fHitClusterArray; }
   void ClearHitClusters() { fHitClusterArray.clear(); }

   ClassDefOverride(AtHitClusterEvent, 1);
};

#endif
