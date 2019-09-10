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
#include "TH1F.h"

#include <string>

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
  std::vector< edm::EDGetTokenT<HGCRecHitCollection> > recHitsTokens_;
  edm::EDGetTokenT<reco::GenParticleCollection> genParticles_;

  hgcal::RecHitTools recHitTools_;

  std::vector<SlimmedHit>* slimmedHits_;
  std::vector<SlimmedROI>* slimmedROIs_;
  TLorentzVector* genVertex_;

  std::vector< TH1F* > zhist;
  std::vector< TH1F* > xhist;
  std::vector< TH1F* > yhist; 
  std::vector< TH1F* > zsidehist;
  std::vector< TH1F* > ehist;
  std::vector< TH1F* > layerhist;
  std::vector< TH1F* > offsetlayerhist;

  TTree* tree_;
  std::vector< std::string > treenames_;

  std::vector<double> thicknessCorrection_;
  std::vector<double> dEdXWeights_;
  std::vector<std::string> geometrySource_;
  std::vector<double> distancesSR1_, distancesSR2_, distancesSR3_;
  int nControlRegions_;
  std::string particle_;
  unsigned int nSubDets_;
  int particleId_;
  bool byClosest_;
};
 

#endif
