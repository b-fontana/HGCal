#ifndef _slimmedroi_h_
#define _slimmedroi_h_

#include "TObject.h"
#include "TLorentzVector.h"

/**
   @short use to strip down the information of a Region Of Interest
 */
class SlimmedROI : public TObject
{
 public:
 SlimmedROI(): pt_(0), eta_(0), phi_(0), mass_(0), pdgid_(0) {}
 SlimmedROI(float pt, float eta, float phi,float mass,int pdgid): pt_(pt), eta_(eta), phi_(phi), mass_(mass), pdgid_(pdgid) {}
  SlimmedROI(const SlimmedROI &other)
    {
      pt_       = other.pt_;            
      eta_      = other.eta_; 
      phi_      = other.phi_; 
      mass_     = other.mass_;
      pdgid_       = other.pdgid_;    
    }
  virtual ~SlimmedROI() { }

  float pt()  { return pt_;   }
  float eta() { return eta_;  }
  float phi() { return phi_;  }
  float m()   { return mass_; }
  int pdgid()    { return pdgid_;   }
  TLorentzVector p4() { 
    TLorentzVector p(0,0,0,0);
    p.SetPtEtaPhiM(pt_,eta_,phi_,mass_);
    return p;
  }

  float pt_,eta_,phi_,mass_;
  int pdgid_;

  ClassDef(SlimmedROI,2)
};

typedef std::vector<SlimmedROI> SlimmedROICollection;


#endif
