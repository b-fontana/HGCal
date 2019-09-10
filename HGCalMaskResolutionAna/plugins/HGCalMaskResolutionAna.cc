#include "UserCode/HGCalMaskResolutionAna/plugins/HGCalMaskResolutionAna.h"
#include "UserCode/HGCalMaskResolutionAna/interface/SimTools.h"

#include "DetectorDescription/OfflineDBLoader/interface/GeometryInfoDump.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "SimG4CMS/Calo/interface/CaloHitID.h"

#include "DetectorDescription/Core/interface/DDFilter.h"
#include "DetectorDescription/Core/interface/DDFilteredView.h"
#include "DetectorDescription/Core/interface/DDSolid.h"

#include "DataFormats/GeometryVector/interface/Basic3DVector.h"

#include "CLHEP/Units/GlobalSystemOfUnits.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Geometry/Plane3D.h"
#include "CLHEP/Geometry/Vector3D.h"

#include "TVector2.h"

#include <iostream>

HGCalMaskResolutionAna::HGCalMaskResolutionAna( const edm::ParameterSet &ps ) : 
  mc_( consumes<edm::HepMCProduct>(edm::InputTag("generatorSmeared")) ),
  recHitsTokens_( {consumes<HGCRecHitCollection>(ps.getParameter<edm::InputTag>("recHitsCEEToken")),
	consumes<HGCRecHitCollection>(ps.getParameter<edm::InputTag>("recHitsHSiToken")),
	consumes<HGCRecHitCollection>(ps.getParameter<edm::InputTag>("recHitsHScToken"))} ),
  genParticles_( consumes<std::vector<reco::GenParticle>>(edm::InputTag("genParticles"))),
  slimmedHits_(new std::vector<SlimmedHit>),
  slimmedROIs_(new std::vector<SlimmedROI>),
  genVertex_(new TLorentzVector),
  treenames_({"CEE", "HEF", "HEB"}),
  dEdXWeights_(ps.getParameter<std::vector<double>>("dEdXWeights")),
  geometrySource_(ps.getParameter<std::vector<std::string>>("geometrySource")),
  distancesSR1_(ps.getParameter<std::vector<double>>("distancesSR1")),
  distancesSR2_(ps.getParameter<std::vector<double>>("distancesSR2")),
  distancesSR3_(ps.getParameter<std::vector<double>>("distancesSR3")),
  nControlRegions_(ps.getParameter<int>("nControlRegions")),
  particle_(ps.getParameter<std::string>("particle")),
  byClosest_(ps.getParameter<bool>("byClosest"))
{  
  if(particle_ == "pion")
    {
      nSubDets_ = 3;
      particleId_ = 211;
    }
  else if(particle_ == "photon")
    {
      nSubDets_ = 1;
      particleId_ = 22;
    }
  else
    {
      std::cout << "The only particles supported are photons and pions!" << std::endl;
      std::exit(0);
    }
  assert(geometrySource_.size() == nSubDets_);

  edm::Service<TFileService> fs;
  for(size_t idet = 0; idet<nSubDets_; ++idet) { 
    zhist.push_back(fs->make<TH1F>(("z"+treenames_[idet]).c_str(), ("z"+treenames_[idet]).c_str(), 2000, -10, 650));
    xhist.push_back(fs->make<TH1F>(("x"+treenames_[idet]).c_str(), ("x"+treenames_[idet]).c_str(), 1000, -40, 40));
    yhist.push_back(fs->make<TH1F>(("y"+treenames_[idet]).c_str(), ("y"+treenames_[idet]).c_str(), 1000, -40, 40));
    zsidehist.push_back(fs->make<TH1F>(("zside"+treenames_[idet]).c_str(), ("zside"+treenames_[idet]).c_str(), 4, -1.1, 1.1));
    ehist.push_back(fs->make<TH1F>(("e"+treenames_[idet]).c_str(), ("e"+treenames_[idet]).c_str(), 2000, -10, 650));
    layerhist.push_back(fs->make<TH1F>(("layer"+treenames_[idet]).c_str(), ("layer"+treenames_[idet]).c_str(), 102, 0, 52));
    offsetlayerhist.push_back(fs->make<TH1F>(("offsetlayer"+treenames_[idet]).c_str(), ("offsetlayer"+treenames_[idet]).c_str(), 102, 0, 52));
  }

  tree_ = fs->make<TTree>((treenames_[0]+"_"+treenames_[1]+"_"+treenames_[2]).c_str(), 
			  (treenames_[0]+"_"+treenames_[1]+"_"+treenames_[2]).c_str());
  tree_->Branch("Hits", "std::vector<SlimmedHit>", &slimmedHits_);
  tree_->Branch("ROIs", "std::vector<SlimmedROI>", &slimmedROIs_);
  tree_->Branch("GenVertex", "TLorentzVector", &genVertex_);
}

//
HGCalMaskResolutionAna::~HGCalMaskResolutionAna()
{
  delete slimmedHits_;
  delete slimmedROIs_;
  delete genVertex_;
}

void HGCalMaskResolutionAna::endJob()
{
}

GlobalPoint HGCalMaskResolutionAna::projectHitPositionAt(float z,float eta,float phi){
  float theta = 2*TMath::ATan(exp(-eta));
  float rho = z*TMath::Tan(theta);
  GlobalPoint xyz(rho*TMath::Cos(phi), rho*TMath::Sin(phi), z);
  return xyz;
}
  
void HGCalMaskResolutionAna::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
  recHitTools_.getEventSetup(iSetup);
  size_t nparticles = 2;
  slimmedROIs_->clear();
  slimmedHits_->clear();
  std::vector< edm::Handle<HGCRecHitCollection> > recHitsHandles(nSubDets_);
  std::vector< edm::ESHandle<HGCalGeometry> > geomHandle(nSubDets_);

  //read the generated primary vertex
  edm::Handle<edm::HepMCProduct> mcHandle;
  iEvent.getByToken(mc_, mcHandle);
  HepMC::GenVertex *primaryVertex = *(mcHandle)->GetEvent()->vertices_begin();
  genVertex_->SetXYZT(primaryVertex->position().x() / 10.,
		      primaryVertex->position().y() / 10.,
		      primaryVertex->position().z() / 10.,
		      primaryVertex->position().t());
    
  //read gen level particles
  std::vector<size_t> selGenPartIdx;
  edm::Handle<std::vector<reco::GenParticle> > genParticlesHandle;
  iEvent.getByToken(genParticles_, genParticlesHandle);
  for(size_t i=0; i<genParticlesHandle->size(); ++i)  {    
    const reco::GenParticle& p = (*genParticlesHandle)[i];
    if(abs(p.pdgId())!=particleId_) continue;
    SlimmedROI particle_roi(p.pt(),p.eta(),p.phi(),p.mass(),p.pdgId());
    slimmedROIs_->push_back(particle_roi);
    selGenPartIdx.push_back(i);
  }
  if(slimmedROIs_->size()!=nparticles) 
    {
      std::cout << "Does not have two particles! It has " << slimmedROIs_->size() <<  "." << std::endl;
      std::exit(0);
    }

  //shift ROIs in phi to create control ROIs for noise/pileup control
  float dphi = 2.*TMath::Pi()/(static_cast<float>(nControlRegions_)+1.);
  for(size_t ir=0; ir<nparticles; ir++) 
    {
      TLorentzVector p4 = slimmedROIs_->at(ir).p4();
      for(int iphi=1; iphi<=nControlRegions_; iphi++) 
	{
	  float newPhi = TVector2::Phi_mpi_pi(p4.Phi()+iphi*dphi);
	  SlimmedROI noise(0,p4.Eta(),newPhi,0,-(ir+1));
	  slimmedROIs_->push_back(noise);
	}
    }

  ///CEE, HFE (Silicon), HBF (Scintillator)///
  for(size_t idet=0; idet<nSubDets_; ++idet) {
    iSetup.get<IdealGeometryRecord>().get(geometrySource_[idet], geomHandle[idet]);
  
    //collect rec hits in regions of interest
    iEvent.getByToken(recHitsTokens_[idet], recHitsHandles[idet]);
    int nMatched(0), nUnmatched(0);
    for(size_t i=0; i<recHitsHandles[idet]->size(); i++) {
      const HGCRecHit &h = recHitsHandles[idet]->operator[](i);
      const HGCalDetId did = h.detid();
      GlobalPoint xyz = recHitTools_.getPosition(did);
      TVector2 xy(xyz.x(), xyz.y());
      float z = xyz.z();
      zhist[idet]->Fill(z);
      //check if it is within ROI assuming linear propagation from genVertex
      int matchedROIidx(-1), signalRegion(-1), nInsideSR(0);
      for(size_t ir=0; ir<slimmedROIs_->size(); ir++)
	{
	  float eta = slimmedROIs_->operator[](ir).eta();
	  if(z==0) 
	    {
	      ehist[idet]->Fill(h.energy());
	      xhist[idet]->Fill(xyz.x());
	      yhist[idet]->Fill(xyz.y());
	      zsidehist[idet]->Fill(did.zside());
	      layerhist[idet]->Fill(did.layer());
	      offsetlayerhist[idet]->Fill(recHitTools_.getLayerWithOffset(did));
	      continue;
	    }
	  if(eta*z < 0)
	    continue;
	  float phi = slimmedROIs_->operator[](ir).phi();      
	  GlobalPoint xyzExp = projectHitPositionAt(z-genVertex_->Z(), eta, phi);      
	  TVector2 xyExp(xyzExp.x(), xyzExp.y());
	  float d = (xyExp-xy).Mod();

	  //get closest detid if required
	  if(byClosest_) 
	    {
	      DetId closest_did = geomHandle[idet].product()->getClosestCell(xyzExp);
	      GlobalPoint xyzClosest(recHitTools_.getPosition(closest_did));
	      TVector2 xyClosest(xyzClosest.x(),xyzClosest.y());        
	      d = (xyClosest-xy).Mod();
	    }

	  if(d > distancesSR3_[idet]) 
	    continue;
	  nMatched += 1;
	  matchedROIidx = ir; 
	  signalRegion = 3;
	  if(d <= distancesSR2_[idet]) 
	    signalRegion = 2;
	  if(d <= distancesSR1_[idet]) 
	      signalRegion = 1;
	 
	  nInsideSR += 1;
	  if(nInsideSR > 2)
	    {
	      std::cout << "A second signal region captured the RecHit!" << std::endl;
	      std::exit(0);
	    }
	}
      if(matchedROIidx < 0) 
	{
	  nUnmatched++;
	  continue;
	}

      //store information
      unsigned int layer = recHitTools_.getLayerWithOffset(did);
      float en = h.energy();
      double mip = dEdXWeights_.at(layer) * 0.001;  // convert to GeV, weights are in MeV
      
      float en_mips = en/static_cast<float>(mip);
      if(en_mips<1) continue;

      float t = h.time();
      SlimmedHit sh(layer, matchedROIidx, idet+1, signalRegion, en, en_mips, t, 0);
      slimmedHits_->push_back(sh);
    }
    std::cout << "\tTotal SR3: " << nMatched << "; Unmatched: " << nUnmatched << " (" <<
      static_cast<float>(nUnmatched)/recHitsHandles[idet]->size() << "); Detector: " << treenames_[idet] << std::endl;
  }
  tree_->Fill();
}

//define this as a plug-in
DEFINE_FWK_MODULE(HGCalMaskResolutionAna);
