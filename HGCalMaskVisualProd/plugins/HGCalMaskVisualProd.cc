#include "UserCode/HGCalMaskVisualProd/plugins/HGCalMaskVisualProd.h"

HGCalMaskVisualProd::HGCalMaskVisualProd(const edm::ParameterSet& iConfig):
  recHitsToken_(consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit", "HGCEERecHits"))),
  cellFilterCuts_ (std::make_pair(iConfig.getParameter<double>("lCellFilterCut"),
				  iConfig.getParameter<double>("hCellFilterCut"))),
  layersAnalysed_(iConfig.getParameter<std::vector<int_layer>>("LayersAnalysed")),
  mask_(iConfig.getParameter<unsigned int>("Mask"))
{
}

HGCalMaskVisualProd::~HGCalMaskVisualProd()
{
}


// ------------ method called to produce the data  ------------
void
HGCalMaskVisualProd::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<HGCRecHitCollection> recHitsHandle;
  iEvent.getByToken(recHitsToken_, recHitsHandle);
  const auto &recHits = *recHitsHandle;

  for(const auto &recHit : recHits){
    HGCSiliconDetId sid(recHit.detid());
    int_layer det_layer = static_cast<int_layer>(sid.layer());

    //mask
    if(cellMask(sid))
      continue;

    //store the data in case the RecHit was measured in one 
    //of the user's chosen layersAnalysed
    if(std::find(layersAnalysed_.begin(), layersAnalysed_.end(), det_layer) != layersAnalysed_.end()) {
      std::pair<int,int> cellUV(sid.cellUV());
      std::pair<int,int> waferUV(sid.waferUV());
      int waferId = linearUV(waferUV.first, waferUV.second);
      //fill only if the histogram/wafer satisfies the cellFilter() criteria
      if(histosRecHits_[det_layer].count(waferId)) 
	histosRecHits_[det_layer][waferId]->Fill(cellUV.first, cellUV.second, recHit.energy());
    }
  }
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
HGCalMaskVisualProd::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
HGCalMaskVisualProd::endStream() {
}

// ------------ method called when starting to processes a run  ------------
void
HGCalMaskVisualProd::beginRun(edm::Run const&, edm::EventSetup const& es)
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

  //create directories and files, one per layer being analysed
  for(const auto &itLayer : layersAnalysed_)
    layersAnalysedDirs_.insert( std::make_pair(itLayer, fs_->mkdir("layer"+std::to_string(itLayer))) );

  createHistograms();
}

 
// ------------ method called when ending the processing of a run  ------------
/*
  void
  HGCalMaskVisualProd::endRun(edm::Run const&, edm::EventSetup const&)
  {
  }
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
  void
  HGCalMaskVisualProd::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
  {
  }
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
  void
  HGCalMaskVisualProd::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
  {
  }
*/

bool HGCalMaskVisualProd::cellFilter(const HGCSiliconDetId &s) {
  if(!gHGCal_) 
    throw std::invalid_argument("No geometry was selected.");
  double z = gHGCal_->topology().dddConstants().waferZ(s.layer(), true);
  std::pair<double,double> r = gHGCal_->topology().dddConstants().rangeR(z, true);
  GlobalPoint p = gHGCal_->getPosition(s);  
  return p.perp() < (r.first + cellFilterCuts_.first) || 
    p.perp() > (r.second - cellFilterCuts_.second);
}

bool HGCalMaskVisualProd::cellMask(const HGCSiliconDetId &s) const {
  if(!gHGCal_) 
    throw std::invalid_argument("No geometry was selected.");
  return gHGCal_->topology().dddConstants().maskCell(s, mask_);
}

void HGCalMaskVisualProd::createHistograms() {
  if(!gHGCal_) 
    throw std::invalid_argument("No geometry was selected.");
  const std::vector<DetId>& ids = gHGCal_->getValidDetIds(); //size: 4126448

  for(const auto &it : ids) {
    HGCSiliconDetId sid(it);
    int_layer det_layer = static_cast<int_layer>(sid.layer());
    std::pair<int,int> uv = sid.waferUV();
    int N = gHGCal_->topology().dddConstants().getUVMax(sid.type());
    int N2 = 2*N;
    int waferId = linearUV(uv.first, uv.second);

    if(std::find(layersAnalysed_.begin(), layersAnalysed_.end(), det_layer) 
       != layersAnalysed_.end() ) {
    
      //create only histograms that correspond to a wafer that has at least one cell that
      //satisfies the selection
      if(cellFilter(sid)) {

	//if histosRecHits_[layer] does not exist, it will be initialized
	if (!histosRecHits_[det_layer].count(waferId)) { 
	  umap<int, TH2F> h_tmp;
	  std::string nn = std::to_string(uv.first)+","+std::to_string(uv.second)+","+
	    std::to_string(N)+",RecHits";
	  histosRecHits_[det_layer].insert(std::make_pair(waferId,
			  layersAnalysedDirs_[det_layer].make<TH2F>(nn.c_str(),nn.c_str(),
								    N2,0.,N2,
								    N2,0.,N2)));
	}

	//if histosGeom_[layer] does not exist, it will be initialized
	if (!histosGeom_[det_layer].count(waferId)) { 
	  umap<int, TH2F> h_tmp;
	  std::string nn = std::to_string(uv.first)+","+std::to_string(uv.second)+","+
	    std::to_string(N)+",Geom";
	  histosGeom_[det_layer].insert(std::make_pair(waferId,
			 layersAnalysedDirs_[det_layer].make<TH2F>(nn.c_str(),nn.c_str(),
								   N2,0.,N2,
								   N2,0.,N2)));
	}
      }

      //filling geometry histograms here avoids looping through the DetIds more than once
      std::pair<int,int> celluv = sid.cellUV();
      fillGeomHistograms(det_layer, waferId, celluv);
    }
  }
}

void HGCalMaskVisualProd::fillGeomHistograms(int_layer layer, int wId, std::pair<int,int> cUV) {
  //if the histogram already exists
  if (histosGeom_[layer].count(wId))
    histosGeom_[layer][wId]->Fill(cUV.first, cUV.second);	
}

//define this as a plug-in
DEFINE_FWK_MODULE(HGCalMaskVisualProd);
