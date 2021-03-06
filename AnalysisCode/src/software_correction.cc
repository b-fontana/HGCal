#include "UserCode/AnalysisCode/interface/software_correction.h"

//classify showers as complete or incomplete according to their energy per layer distributions
//0: complete; 1, 2, ...: each background category
int_ diff_ed(const vec1d<float_>& curr, const vec1d<float_>& standard, 
	     const vec1d<float_>& thresholds, const float_& minval) {
  assert(curr.size()==standard.size());
  float_ cumdiff = 0.;
  typename vec1d<float_>::const_iterator it1;
  for(it1 = curr.cbegin(); it1!=curr.cend(); ++it1) {
    int_ idx = std::distance(curr.cbegin(), it1);
    if(standard[idx] < minval)
      continue;
    cumdiff += fabs(standard[idx]-curr[idx]);
  }
  if(cumdiff<thresholds[0])
    return 0;
  typename vec1d<float_>::const_iterator it2;
  for(it2 = thresholds.cbegin(); it2!=(thresholds.cend()-1); ++it2) {
    int_ idx = std::distance(thresholds.cbegin(), it2);
    if(cumdiff>=*it2 && cumdiff<*(it2+1))
      return idx+1;
  }
  uint_ n = thresholds.size();
  if(cumdiff>=thresholds[n-1])
    return n;
  //it should never reach this point
  return -1;
}

SoftwareCorrector::SoftwareCorrector(std::string fname) {
  wn = {{"weight1_sr1", "weight2_sr1", "weight3_sr1"},
	{"weight1_sr2", "weight2_sr2", "weight3_sr2"},
	{"weight1_sr3", "weight2_sr3", "weight3_sr3"}};
  hn_sig = {"en1_layer_sign", "en2_layer_sign", "en3_layer_sign"};
  fw = new TFile(fname.c_str(), "READ");
  TH1F* h;
  for(int_ ireg=0; ireg<nreg; ++ireg) {
    for(uint_ iw=0; iw<discrvals.size(); ++iw) {
      h = static_cast<TH1F*>(fw->Get(wn[ireg][iw].c_str()));
      for(int_ il=1; il<=nlayers; ++il) {
	int_ b = h->FindBin(il);
	weights(ireg,iw,il-1) = h->GetBinContent(b);
      }
    }
  }
}

SoftwareCorrector::~SoftwareCorrector() {
  fw->Close();
  delete fw;
}

vec1d<TGraph*> SoftwareCorrector::build_weights_graphs(int_ region) {
  vec1d<TGraph*> g;
  for(int_ il=0; il<nlayers; ++il) {
    g[il] = new TGraph(discrvals.size());
    for(uint_ iw=0; iw<discrvals.size(); ++iw) {
      float_ x = discrvals[iw];
      float_ y = this->weights(region-1,iw,il);
      g[il]->SetPoint(iw, x, y);
    }
  }
  return g;
}

vec1d<float_> SoftwareCorrector::low_stats_factor(vec1d<uint_> limits, std::string mode) {
  vec1d<float_> f(nreg, 1.);
  TH1F* h;
  for(int_ ireg=0; ireg<nreg; ++ireg) {
    assert(limits.size()==wn.size());
    h = static_cast<TH1F*>(fw->Get(hn_sig[ireg].c_str()));
    if(h->Integral()==0) {
      std::cout << "The integral is zero." << std::endl;
      std::exit(0);
    }
    uint_ limita = 0;
    uint_ limitb = 0;
    if(mode=="left") {
      limita = h->FindBin(1);
      limitb = h->FindBin(limits[ireg]);
    }
    else if(mode == "right") {
      limita = h->FindBin(limits[ireg]);
      limitb = h->FindBin(h->GetNbinsX());
    }
    else {
      std::cout << "The mode introduced is not valid." << std::endl;
      std::exit(0);
    }
    f[ireg] = h->Integral(limita,limitb)/h->Integral();
  }
  return f;
}
