#include "UserCode/HGCalMaskProd/plugins/HGCalMaskProd.h"

HGCalMaskProd::HGCalMaskProd(const edm::ParameterSet& ps):
  recHitsTokens_( {consumes<HGCRecHitCollection>(ps.getParameter<edm::InputTag>("recHitsCEEToken")),
	consumes<HGCRecHitCollection>(ps.getParameter<edm::InputTag>("recHitsHSiToken"))} ),
  mask_(ps.getParameter<unsigned int>("Mask"))
{
  for (unsigned int idet=0; idet<2; ++idet)
    produces<HGCRecHitCollection>(colName_[idet]);
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
  ///CEE, HFE (Silicon)///
  for (unsigned int idet=0; idet<2; ++idet)
    {
      std::unique_ptr<HGCRecHitCollection> recHitsMask = std::make_unique<HGCRecHitCollection>();
      edm::Handle<HGCRecHitCollection> recHitsHandle;
      iEvent.getByToken(recHitsTokens_[idet], recHitsHandle);
      const auto &recHits = *recHitsHandle;
      for(const auto &recHit : recHits){
	HGCSiliconDetId sid(recHit.detid()); //HGCScintillatorDetId
	if(cellMask(sid))
	  continue;
	recHitsMask->push_back(recHit);
      }
      iEvent.put(std::move(recHitsMask), colName_.at(idet));
    }
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

// ------------ method called when starting to process a run  ------------
void
HGCalMaskProd::beginRun(edm::Run const&, edm::EventSetup const& es)
{
  edm::ESHandle<CaloGeometry> geom;
  es.get<CaloGeometryRecord>().get(geom);
  Det_ = std::make_pair(DetId::HGCalEE, DetId::HGCalHSi);
  SubDet_ = ForwardSubdetector::ForwardEmpty;
  const CaloSubdetectorGeometry* g1 = geom->getSubdetectorGeometry(Det_.first, SubDet_);
  gEE_ = dynamic_cast<const HGCalGeometry*>(g1);
  const CaloSubdetectorGeometry* g2 = geom->getSubdetectorGeometry(Det_.second, SubDet_);
  gHSi_ = dynamic_cast<const HGCalGeometry*>(g2);
}

bool HGCalMaskProd::cellMask(const HGCSiliconDetId &s) const {
  if( !gEE_ || !gHSi_ ) 
    throw std::invalid_argument("No geometry was selected.");

  bool b = false;
  if( s.isEE() )
    b = gEE_->topology().dddConstants().maskCell(s, mask_);
  else if( s.isHE() )
    b = gHSi_->topology().dddConstants().maskCell(s, mask_);
  else
    throw std::invalid_argument("The RecHit has to belong to a silicon subdetector.");
  return b;
}

//define this as a plug-in
DEFINE_FWK_MODULE(HGCalMaskProd);
