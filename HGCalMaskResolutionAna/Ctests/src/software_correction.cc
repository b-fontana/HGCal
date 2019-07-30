#include "../interface/software_correction.h"

SoftwareCorrection::SoftwareCorrection(std::string fname) {
  wn = {{"weight1_sr1", "weight2_sr1", "weight3_sr1"},
	{"weight1_sr2", "weight2_sr2", "weight3_sr2"},
	{"weight1_sr3", "weight2_sr3", "weight3_sr3"}};
  hn_sig = {"en1_layer_sign", "en2_layer_sign", "en3_layer_sign"};
  hn_bckg = {{"en1_layer_bckg1", "en1_layer_bckg2", "en1_layer_bckg3"},
	     {"en2_layer_bckg1", "en2_layer_bckg2", "en2_layer_bckg3"},
	     {"en3_layer_bckg1", "en3_layer_bckg2", "en3_layer_bckg3"}};
  fw = new TFile(fname.c_str(), "READ");
  TH1F* h;
  for(int_ ireg=0; ireg<nreg; ++ireg) {
    for(int_ iw=0; iw<discrvals.size(); ++iw) {
      h = static_cast<TH1F*>(fw->Get(wn[ireg][iw].c_str()));
      for(int_ j=0; j<nlayers; ++j) {
	int_ b = h->FindBin(j);
	weights[ireg][iw][j] = h->GetBinContent(b);
      }
    }
  }
}

vec1d<TGraph*> SoftwareCorrection::build_weights_graphs(int_ region) {
  vec1d<TGraph*> g;
  for(int_ il=0; il<nlayers; ++il) {
    g[il] = new TGraph(discrvals.size());
    for(int_ iw=0; iw<discrvals.size(); ++iw) {
      float_ x = discrvals[iw];
      float_ y = this->weights[region-1][iw][il];
      g[il]->SetPoint(iw, x, y);
    }
  }
  return g;
}

vec1d<float_> SoftwareCorrection::low_stats_factor(vec1d<float_> limits, std::string mode) {
  vec1d<float_> f;
  for(int_ ireg=0; ireg<nreg; ++ireg) {
    assert(limits.size()==wn.size());
    TH1F* h = static_cast<TH1F*>(this->fw->Get(this->hn_sig[ireg].c_str()));
    if(h->Integral()==0) {
      std::cout << "The integral is zero." << std::endl;
      std::exit(0);
    }
    uint_ limita;
    uint_ limitb;
    if(mode=="left") {
      limita = h->FindBin(1);
      limitb = h->FindBin(limits[ireg]);
    }
    else if(mode == "right") {
      limita = h->FindBin(limits[ireg]);
      limitb = h->FindBin(h->GetNbinsX());
    }
    f[ireg] = h->Integral(limita,limitb)/h->Integral();
    h->Delete();
  }
  return f;
}

//classify showers as complete or incomplete according to their energy per layer distributions
//0: complete; 1, 2, ...: each background category
int_ SoftwareCorrection::diff_ed(vec1d<float_> curr, vec1d<float_> standard, 
				 vec1d<float_> tresholds, float_ min_val) {

}
