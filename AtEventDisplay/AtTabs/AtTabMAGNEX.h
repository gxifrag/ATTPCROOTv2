#ifndef ATTABMAGNEX_H
#define ATTABMAGNEX_H

#include "AtTabMain.h"
#include "AtViewerManagerSubject.h"
#include "AtHitClusterFull.h"

#include <Rtypes.h>
#include <TEveEventManager.h>
#include <TEvePointSet.h>

class TBuffer;
class TClass;
class TMemberInspector;
namespace DataHandling {
class AtSubject;
}

class AtTabMAGNEX : public AtTabMain {
protected:
   TEveEventManagerPtr fEveHitClusterEvent{std::make_unique<TEveEventManager>("AtHitClusterEvent")};
   TEveEventManagerPtr fEveMeanHits{std::make_unique<TEveEventManager>("Mean Points")};
   TEvePointSetPtr fMeanHitSet{std::make_unique<TEvePointSet>("MeanHits")};
   std::vector<TEvePointSetPtr> fHitClusterSets;

   DataHandling::AtBranch *fHitClusterEventBranch;

public:
   AtTabMAGNEX();
   ~AtTabMAGNEX();

   void InitTab() override;
   void Update(DataHandling::AtSubject *sub) override;

protected:
   void UpdateRenderState() override;

   void SetPointsFromHitCluster(TEvePointSet &hitSet, const AtHitClusterFull &hitCluster);

private:
   void UpdateHitClusterEventElements();
   void UpdatePadPlane() override;

   void ExpandNumHitClusters(Int_t num);

   ClassDefOverride(AtTabMAGNEX, 1);
};

#endif
