#ifndef _HGCalMaskResolutionAna_h_
#define _HGCalMaskResolutionAna_h_

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESTransientHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHit.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

#include "UserCode/HGCalMaskResolutionAna/interface/SlimmedHit.h"
#include "UserCode/HGCalMaskResolutionAna/interface/SlimmedROI.h"

#include "TLorentzVector.h"
#include "TTree.h"

#include <string>

/**
   @class HGCalMaskResolutionAna
   @author P. Silva (CERN)
*/

class HGCalMaskResolutionAna : public edm::EDAnalyzer 
{
  
 public:
  
  explicit HGCalMaskResolutionAna( const edm::ParameterSet& );
  ~HGCalMaskResolutionAna();
  virtual void analyze( const edm::Event&, const edm::EventSetup& );
  void endJob();

 private:

  GlobalPoint projectHitPositionAt(float z,float eta,float phi);

  edm::EDGetTokenT<edm::HepMCProduct> mc_;
  edm::EDGetTokenT<HGCRecHitCollection> recHitsEE_;
  edm::EDGetTokenT<reco::GenParticleCollection> genParticles_;
  edm::EDGetTokenT<std::vector<SimTrack>> simTracks_;
  edm::EDGetTokenT<std::vector<SimVertex>> simVertices_;
  edm::EDGetTokenT<std::vector<int>> genBarCodes_;

  hgcal::RecHitTools recHitTools_;

  TTree *tree_;
  std::vector<SlimmedHit> *slimmedHits_;
  std::vector<SlimmedROI> *slimmedROIs_;
  TLorentzVector *genVertex_;

  std::vector<double> dEdXWeights_, thicknessCorrection_;
  bool byClosest_;
};
 

#endif
