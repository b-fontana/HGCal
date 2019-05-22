#include "UserCode/HGCalLateralStudies/plugins/HGCalLateralStudies.h"

HGCalLateralStudies::HGCalLateralStudies(const edm::ParameterSet& iConfig):
  recHitsToken_(consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit", "HGCEERecHits"))),
  cellFilterCuts_ (std::make_pair(iConfig.getParameter<double>("lCellFilterCut"),
				  iConfig.getParameter<double>("hCellFilterCut"))),
  layersAnalysed_(iConfig.getParameter<std::vector<int_layer>>("LayersAnalysed")),
  mask_(iConfig.getParameter<unsigned int>("Mask")),
  CellUVCollection_name_(iConfig.getParameter<std::string>("CellUVCoordinates")),
  WaferUVCollection_name_(iConfig.getParameter<std::string>("WaferUVCoordinates"))
{
  produces<CellUVCollection_>(CellUVCollection_name_);  
  produces<CellUVCollection_>(WaferUVCollection_name_);  
}

HGCalLateralStudies::~HGCalLateralStudies()
{
}

//
// member functions
//

// ------------ method called to produce the data  ------------
void
HGCalLateralStudies::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  std::unique_ptr<CellUVCollection_> celluv_coords = std::make_unique<CellUVCollection_>();
  std::unique_ptr<WaferUVCollection_> waferuv_coords = std::make_unique<WaferUVCollection_>();

  edm::Handle<HGCRecHitCollection> recHitsHandle;
  iEvent.getByToken(recHitsToken_, recHitsHandle);
  const auto &recHits = *recHitsHandle;
  const int size = recHitsHandle->size();
  celluv_coords->reserve(size);
  waferuv_coords->reserve(size);

  for(const auto &recHit : recHits){
    HGCSiliconDetId sid(recHit.detid());
    int_layer det_layer = static_cast<int_layer>(sid.layer());

    //mask and positive side of endcap checks
    if(cellMask(sid) || sid.zside()==1) 
      continue;

    //store the data in case the RecHit was measured in one of the user's chosen layersAnalysed
    if(std::find(layersAnalysed_.begin(), layersAnalysed_.end(), det_layer) != layersAnalysed_.end()) {
      std::pair<int,int> cellUV(sid.cellUV());
      std::pair<int,int> waferUV(sid.waferUV());
      int waferId = linearUV(waferUV.first, waferUV.second);
      //fill only if the histogram/wafer satisfies the cellFilter() criteria
      if(histosRecHits_[det_layer].count(waferId)) 
	histosRecHits_[det_layer][waferId]->Fill(cellUV.first, cellUV.second);
    }
  }
  FileUtils::close(outRecHits_);

  //iEvent.put(std::move(celluv_coords), CellUVCollection_name);
  //iEvent.put(std::move(waferuv_coords), WaferUVCollection_name);

  /* this is an EventSetup example
  //Read SetupData from the SetupRecord in the EventSetup
  ESHandle<SetupData> pSetup;
  iSetup.get<SetupRecord>().get(pSetup);
  */
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
HGCalLateralStudies::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
HGCalLateralStudies::endStream() {
}

// ------------ method called when starting to processes a run  ------------
void
HGCalLateralStudies::beginRun(edm::Run const&, edm::EventSetup const& es)
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
  for(const auto &itLayer : layersAnalysed_) {
    layersAnalysedDirs_.insert( std::make_pair(itLayer, fs_->mkdir("layer"+std::to_string(itLayer))) );
    outRecHits_.setName("data/HistoNamesRecHitsLayer"+
			std::to_string(itLayer)+'_'+std::to_string(mask_),
			static_cast<int>(itLayer));
    FileUtils::create(outRecHits_, static_cast<int>(itLayer));
    outGeom_.setName("data/HistoNamesGeomLayer"+
		     std::to_string(itLayer)+'_'+std::to_string(mask_),
		     static_cast<int>(itLayer));
    FileUtils::create(outGeom_, static_cast<int>(itLayer));
  }

  createHistograms();
}

 
// ------------ method called when ending the processing of a run  ------------
/*
  void
  HGCalLateralStudies::endRun(edm::Run const&, edm::EventSetup const&)
  {
  }
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
  void
  HGCalLateralStudies::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
  {
  }
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
  void
  HGCalLateralStudies::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
  {
  }
*/

bool HGCalLateralStudies::cellFilter(const HGCSiliconDetId &s) {
  if(!gHGCal_) 
    throw std::invalid_argument("No geometry was selected.");
  double z = gHGCal_->topology().dddConstants().waferZ(s.layer(), true);
  std::pair<double,double> r = gHGCal_->topology().dddConstants().rangeR(z, true);
  GlobalPoint p = gHGCal_->getPosition(s);  
  return p.perp() < (r.first + cellFilterCuts_.first) || 
    p.perp() > (r.second - cellFilterCuts_.second);
}

bool HGCalLateralStudies::cellMask(const HGCSiliconDetId &s) const {
  if(!gHGCal_) 
    throw std::invalid_argument("No geometry was selected.");
  /*
  const auto cor = ((cornerType_ > 0) ? geom->getCorners(id) : 
		    ((cornerType_ < 0) ? geom->get8Corners(id) :
		     geom->getNewCorners(id)));
  */
  return gHGCal_->topology().dddConstants().maskCell(s, 3);
}

void HGCalLateralStudies::createHistograms() {
  if(!gHGCal_) 
    throw std::invalid_argument("No geometry was selected.");
  const std::vector<DetId>& ids = gHGCal_->getValidDetIds(); //size: 4126448

  for(const auto &it : ids) {
    HGCSiliconDetId sid(it);
    int_layer det_layer = static_cast<int_layer>(sid.layer());
    std::pair<int,int> uv = sid.waferUV();
    int waferId = linearUV(uv.first, uv.second);

    if(std::find(layersAnalysed_.begin(), layersAnalysed_.end(), det_layer) 
       != layersAnalysed_.end() ) {
    
      //create only histograms that correspond to a wafer that has at least one cell that
      //satisfies the selection
      if(cellFilter(sid)) {

	//if histosRecHits_[layer] does not exist, it will be initialized
	if (!histosRecHits_[det_layer].count(waferId)) { 
	  umap<int, TH2F> h_tmp;
	  std::string nn = std::to_string(uv.first)+","+std::to_string(uv.second)+",RecHits";
	  histosRecHits_[det_layer].insert(std::make_pair(waferId,
			  layersAnalysedDirs_[det_layer].make<TH2F>(nn.c_str(),nn.c_str(),
								    25,0,24,25,0,24)));
	  FileUtils::reopen(outRecHits_, static_cast<int>(det_layer), std::ios_base::app);
	  std::string write_str = std::to_string(uv.first)+","+std::to_string(uv.second)+",RecHits";
	  FileUtils::write(outRecHits_, write_str);
	}

	//if histosGeom_[layer] does not exist, it will be initialized
	if (!histosGeom_[det_layer].count(waferId)) { 
	  umap<int, TH2F> h_tmp;
	  std::string nn = std::to_string(uv.first)+","+std::to_string(uv.second)+",Geom";
	  histosGeom_[det_layer].insert(std::make_pair(waferId,
			       layersAnalysedDirs_[det_layer].make<TH2F>(nn.c_str(),nn.c_str(),
									 25,0,24,25,0,24)));
	  FileUtils::reopen(outGeom_, static_cast<int>(det_layer), std::ios_base::app);
	  std::string write_str = std::to_string(uv.first)+","+std::to_string(uv.second)+",Geom";
	  FileUtils::write(outGeom_, write_str);
	}
      }

      //filling the geometry histograms here avoids looping through the DetIds more than once
      std::pair<int,int> celluv = sid.cellUV();
      fillGeomHistograms(det_layer, waferId, celluv);
    }
  }
  FileUtils::close(outGeom_); //the Geometry histograms are not needed in the produce() stage
}

void HGCalLateralStudies::fillGeomHistograms(int_layer layer, int wId, std::pair<int,int> cUV) {
  //if the histogram already exists
  if (histosGeom_[layer].count(wId))
    histosGeom_[layer][wId]->Fill(cUV.first, cUV.second);	
}

void HGCalLateralStudies::perpHistogram(GlobalPoint p) {
  static bool first_time = true;
  static TH1D *h_perp;
  if(first_time) {
    TFileDirectory subdir_perp = fs_->mkdir("perpHistogram");   
    h_perp = subdir_perp.make<TH1D>("h_perp", "h_perp", 50, 0, 200);
  }
  h_perp->Fill(p.perp());
  first_time = false;
}

//define this as a plug-in
DEFINE_FWK_MODULE(HGCalLateralStudies);
