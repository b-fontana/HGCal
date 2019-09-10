#ifndef HGCalMaskProd_h
#define HGCalMaskProd_h

// system include files
#include <iostream>
#include <fstream>
#include <vector>
#include <utility> //std::pair
#include <algorithm> //std::find
#include "TH2F.h"
#include <memory>

// cmssw include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHit.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "DataFormats/ForwardDetId/interface/HGCSiliconDetId.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/EDPutToken.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

//
// class declaration
//
class HGCalMaskProd : public edm::stream::EDProducer<> {
public:
  explicit HGCalMaskProd(const edm::ParameterSet&);
  ~HGCalMaskProd();
  
private:
  //aliases
  typedef unsigned int int_layer;
  template<class T,class U> using umap = std::unordered_map<T,U>;

  //variables
  const edm::Service<TFileService> fs_;
  std::vector< edm::EDGetTokenT<HGCRecHitCollection> > recHitsTokens_;
  const HGCalGeometry *gEE_, *gHSi_;
  std::pair<DetId::Detector, DetId::Detector> Det_; 
  ForwardSubdetector SubDet_;
  const int mask_ = 0;
  std::vector<std::string> colName_ = {"HGCEERecHits", "HGCHSiRecHits"};

  //functions
  virtual void beginStream(edm::StreamID) override;
  virtual void produce(edm::Event&, const edm::EventSetup&) override;
  virtual void endStream() override;
  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  bool cellMask(const HGCSiliconDetId&) const;
};

#endif //HGCalMaskProd_h
