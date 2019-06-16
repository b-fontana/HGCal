#include "UserCode/HGCalMaskProd/plugins/HGCalMaskProd.h"

HGCalMaskProd::HGCalMaskProd(const edm::ParameterSet& iConfig):
  recHitsToken_(consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit", "HGCEERecHits"))),
  layersAnalysed_(iConfig.getParameter<std::vector<int_layer>>("LayersAnalysed")),
  mask_(iConfig.getParameter<unsigned int>("Mask"))
{
  produces<HGCRecHitCollection>(collectionName);  
}

HGCalMaskProd::~HGCalMaskProd()
{
}

//
// member functions
//

// ------------ method called to produce the data  ------------
void
HGCalMaskProd::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  std::unique_ptr<HGCRecHitCollection> recHitsMaskColl = std::make_unique<HGCRecHitCollection>();

  edm::Handle<HGCRecHitCollection> recHitsHandle;
  iEvent.getByToken(recHitsToken_, recHitsHandle);
  const auto &recHits = *recHitsHandle;

  for(const auto &recHit : recHits){
    HGCSiliconDetId sid(recHit.detid());
    if(cellMask(sid))
      continue;
    recHitsMaskColl->push_back(recHit);
  }    
  iEvent.put(std::move(recHitsMaskColl), collectionName);

  /* this is an EventSetup example
  //Read SetupData from the SetupRecord in the EventSetup
  ESHandle<SetupData> pSetup;
  iSetup.get<SetupRecord>().get(pSetup);
  */
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
HGCalMaskProd::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
HGCalMaskProd::endStream() {
}

// ------------ method called when starting to processes a run  ------------
void
HGCalMaskProd::beginRun(edm::Run const&, edm::EventSetup const& es)
{
  edm::ESHandle<CaloGeometry> geom;
  es.get<CaloGeometryRecord>().get(geom);
  
  for(const auto &it : layersAnalysed_) {
    if(it > lastLayerEE_) 
      throw std::domain_error("The chosen layer does not correspond to the HGCalEE.");
    else {
      myDet_=DetId::HGCalEE;
      mySubDet_=ForwardSubdetector::ForwardEmpty;
    }
  }
    
  if(myDet_ == DetId::HGCalEE) { // || myDet_==DetId::HGCalHSi) 
    gHGCal_ = dynamic_cast<const HGCalGeometry*>(geom->getSubdetectorGeometry(myDet_, mySubDet_));
  }
  else {
    gHGCal_ = nullptr;
    throw std::domain_error("Currently only the HGCalEE is supported.");
  }
}

bool HGCalMaskProd::cellMask(const HGCSiliconDetId &s) const {
  if(!gHGCal_) 
    throw std::invalid_argument("No geometry was selected.");
  return gHGCal_->topology().dddConstants().maskCell(s, mask_);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HGCalMaskProd);
