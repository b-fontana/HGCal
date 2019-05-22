#include "UserCode/HGCalRecHitsMask/plugins/HGCalRecHitsMask.h"

HGCalRecHitsMask::HGCalRecHitsMask(const edm::ParameterSet& iConfig):
  recHitsToken_(consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit","HGCEERecHits"))),
  mask_(iConfig.getParameter<unsigned int>("Mask"))
{
  produces<HGCRecHitCollection>(HGCRecHitMaskedCollection_name_);  
}

HGCalRecHitsMask::~HGCalRecHitsMask()
{
}

//
// member functions
//

// ------------ method called to produce the data  ------------
void
HGCalRecHitsMask::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  auto recHitsMasked = std::make_unique<HGCRecHitCollection>();

  edm::Handle<HGCRecHitCollection> recHitsHandle;
  iEvent.getByToken(recHitsToken_, recHitsHandle);
  const auto &recHits = *recHitsHandle;
  const int size = recHitsHandle->size();
  recHitsMasked->reserve(size);

  for(const auto &recHit : recHits){
    HGCSiliconDetId sid(recHit.detid());
    
    //check that the layer corresponds to the current geometry

    //mask
    if(cellMask(sid))
      continue;

    recHitsMasked->push_back(recHit);
  }
  iEvent.put(std::move(recHitsMasked), HGCRecHitMaskedCollection_name_);

  /* this is an EventSetup example
  //Read SetupData from the SetupRecord in the EventSetup
  ESHandle<SetupData> pSetup;
  iSetup.get<SetupRecord>().get(pSetup);
  */
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
HGCalRecHitsMask::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
HGCalRecHitsMask::endStream() {
}

// ------------ method called when starting to processes a run  ------------
void
HGCalRecHitsMask::beginRun(edm::Run const&, edm::EventSetup const& es)
{
  edm::ESHandle<CaloGeometry> geom;
  es.get<CaloGeometryRecord>().get(geom);
  
  myDet_=DetId::HGCalEE;
  mySubDet_=ForwardSubdetector::ForwardEmpty;
  gHGCal_ = dynamic_cast<const HGCalGeometry*>(geom->getSubdetectorGeometry(myDet_, mySubDet_));
  }

bool HGCalRecHitsMask::cellMask(const HGCSiliconDetId &s) const {
  if(!gHGCal_) 
    throw std::invalid_argument("No geometry was selected.");
  return gHGCal_->topology().dddConstants().maskCell(s, mask_);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HGCalRecHitsMask);
