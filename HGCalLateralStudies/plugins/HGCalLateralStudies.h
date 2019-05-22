#ifndef HGCalLateralStudies_h
#define HGCalLateralStudies_h

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

//user include files
#include "UserCode/HGCalLateralStudies/interface/FileUtils.h"

//
// class declaration
//
class HGCalLateralStudies : public edm::stream::EDProducer<> {
public:
  explicit HGCalLateralStudies(const edm::ParameterSet&);
  ~HGCalLateralStudies();
  
private:
  //aliases
  typedef unsigned int int_layer;
  template<class T,class U> using umap = std::unordered_map<T,U>;

  //variables
  const edm::EDGetTokenT<HGCRecHitCollection> recHitsToken_;  
  const edm::Service<TFileService> fs_;
  const HGCalGeometry* gHGCal_;
  DetId::Detector myDet_; 
  ForwardSubdetector mySubDet_;
  const std::pair<double, double> cellFilterCuts_;
  umap<int_layer, TFileDirectory> layersAnalysedDirs_; //one subdir per analysed layer
  const std::vector<int_layer> layersAnalysed_;
  umap<int_layer, umap<int, TH2F*> > histosRecHits_; //stores the position of RecHits
  umap<int_layer, umap<int, TH2F*> > histosGeom_; //stores the underlying geometry
  FileUtils::outData outRecHits_;
  FileUtils::outData outGeom_;
  const int_layer lastLayerEE_ = 28;
  const int mask_ = 0;
			  
  //variables (outputs)
  const std::string CellUVCollection_name_ = "CellUVCollection";
  typedef std::vector< std::pair<int,int> > CellUVCollection_;
  const std::string WaferUVCollection_name_ = "WaferUVCollection";
  typedef std::vector< std::pair<int,int> > WaferUVCollection_;

  //functions
  virtual void beginStream(edm::StreamID) override;
  virtual void produce(edm::Event&, const edm::EventSetup&) override;
  virtual void endStream() override;
  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  int linearUV(int u, int v) {return v*100+u;}
  bool cellFilter(const HGCSiliconDetId&);
  bool cellMask(const HGCSiliconDetId&) const;
  void createHistograms();
  void fillGeomHistograms(int_layer, int, std::pair<int,int>);

  //debug
  void perpHistogram(GlobalPoint);
};

#endif //HGCalLateralStudies_h
