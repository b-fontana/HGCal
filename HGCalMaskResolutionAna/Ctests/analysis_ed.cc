#include <iostream>
#include <vector>
#include <iterator>
#include "TFile.h"
#include "TTree.h"
#include "ROOT/RDataFrame.hxx"

#include "interface/utils.h"
#include "interface/calibration.h"

int_ main() {
  float_ mingenen = 20.;
  vec1d<float_> etareg = {2.7, 2.94};
  vec1d<float_> etareg_shift = VecOps(etareg).shift();
  uint_ n = etareg.size();
  const int_ nreg = 3;
  const int_ nlayers = 28;
  std::string label = "";
  std::string samples = "inner";
  uint_ mask = 3;
  std::string noPUFile = "../summaries/summary_mask3_inner.root";
  std::string outpath = "../pics_inner/mask3";
  std::array<float_, 3> bckgcuts = {0.5, 0.6, 0.7};
  
  Calibration calibration(mingenen, etareg, mask, label, samples, 
			  mask, noPUFile, outpath);
  calibration.nopu_calibration(6, true);
  vec1d<mapfunc> calib = calibration.calib;
  std::array<std::array<float_, nlayers>, nreg> frac_en = {0.};
  std::array<std::array<int_, nlayers>, nreg> countfrac_en = {0};

  //obtain the fractional energy distribution per layer of complete showers
  auto standard_ed = [&](float_ gen, float_ geta, 
			 vec1d<float_> en, vec1d<float_> noi,
			 vec2d<float_>  en_layer) {
    for(int_ ireg=1; ireg<=nreg; ++ireg) {
      float_ f1=1., f2=0.;
      typename vec1d<float_>::const_iterator it;
      for(it = etareg.cbegin(); it!=(etareg.cend()-1); ++it) {
	int_ idx = it-etareg.cbegin();
	std::string idstr;
	//in case it lies outside the limits of the calibration
	//the event is calibrated with the full calibration region
	if (geta<etareg[0] || geta>etareg[n-1]) {
	  idstr = "sr" + std::to_string(ireg) + "from" + 
	    etastr(std::to_string(etareg[0])) + "to" + etastr(std::to_string(etareg[n-1]));
	}
	else if (geta<*it || geta>*(it+1))
	  continue;
	else {
	  idstr = "sr" + std::to_string(ireg) + "from" + 
	    etastr(std::to_string(*it)) + "to" + etastr(std::to_string(*(it+1)));
	}

	//it should only get here once per event
	if (calib[0].size() > 0) {
	  check_key(calib[0], idstr);
	  f1 /= calib[0][idstr]->Eval(geta)+1.0;
	  if (calib[1].size() > 0) {
	    check_key(calib[1], idstr);
	    f1 /= calib[1][idstr]->Eval(f1*en[ireg])+1.0;
	    if (calib[2].size() > 0) {
	      check_key(calib[2], idstr);
	      f2 = calib[2][idstr]->Eval(noi[2]);
	    }
	  }
	}
	en[ireg-1] = f1*en[ireg-1] - f2;
      }
      assert(f1!=1.);
      for(int_ il=1; il<=nlayers; ++il) {
	float_ v = f1*en_layer[ireg-1][il-1] - f2;
	if( (samples=="inner" && geta<etareg[0]+0.05 || samples=="outer" and geta>1.6)
	    && en[ireg-1] != 0. ) {
	  frac_en[ireg-1][il-1] += (v/en[ireg-1]);
	  countfrac_en[ireg-1][il-1] += 1;
	}
      }
    }
  };

  std::string def1 = 
    "std::vector<float> en = {en_sr1_ROI, en_sr2_ROI, en_sr3_ROI};"
    "return en;";
  std::string def2 = 
    "std::vector<float> noise = {noise_sr1_ROI, noise_sr2_ROI, noise_sr3_ROI};"
    "return noise;";
  std::string def3 = 
    "std::vector<float> en_layer1 = {en_sr1_layer1, en_sr1_layer2, en_sr1_layer3,"
    "en_sr1_layer4, en_sr1_layer5, en_sr1_layer6, en_sr1_layer7, en_sr1_layer8,"
    "en_sr1_layer9, en_sr1_layer10, en_sr1_layer11, en_sr1_layer12, en_sr1_layer13,"
    "en_sr1_layer14, en_sr1_layer15, en_sr1_layer16, en_sr1_layer17, en_sr1_layer18,"
    "en_sr1_layer19, en_sr1_layer20, en_sr1_layer21, en_sr1_layer22, en_sr1_layer23,"
    "en_sr1_layer24, en_sr1_layer25, en_sr1_layer26, en_sr1_layer27, en_sr1_layer28};"
    "std::vector<float> en_layer2 = {en_sr2_layer1, en_sr2_layer2, en_sr2_layer3,"
    "en_sr2_layer4, en_sr2_layer5, en_sr2_layer6, en_sr2_layer7, en_sr2_layer8,"
    "en_sr2_layer9, en_sr2_layer10, en_sr2_layer11, en_sr2_layer12, en_sr2_layer13,"
    "en_sr2_layer14, en_sr2_layer15, en_sr2_layer16, en_sr2_layer17, en_sr2_layer18,"
    "en_sr2_layer19, en_sr2_layer20, en_sr2_layer21, en_sr2_layer22, en_sr2_layer23,"
    "en_sr2_layer24, en_sr2_layer25, en_sr2_layer26, en_sr2_layer27, en_sr2_layer28};"
    "std::vector<float> en_layer3 = {en_sr3_layer1, en_sr3_layer2, en_sr3_layer3,"
    "en_sr3_layer4, en_sr3_layer5, en_sr3_layer6, en_sr3_layer7, en_sr3_layer8,"
    "en_sr3_layer9, en_sr3_layer10, en_sr3_layer11, en_sr3_layer12, en_sr3_layer13,"
    "en_sr3_layer14, en_sr3_layer15, en_sr3_layer16, en_sr3_layer17, en_sr3_layer18,"
    "en_sr3_layer19, en_sr3_layer20, en_sr3_layer21, en_sr3_layer22, en_sr3_layer23,"
    "en_sr3_layer24, en_sr3_layer25, en_sr3_layer26, en_sr3_layer27, en_sr3_layer28};"
    "std::vector< std::vector<float> > en_layer = {en_layer1, en_layer2, en_layer3};"
    "return en_layer;";

  ROOT::EnableImplicitMT(10);
  ROOT::RDF::RSnapshotOptions opts;
  opts.fLazy = true;
  opts.fMode = "RECREATE";
  ROOT::RDataFrame d1("data", noPUFile.c_str());
  d1.Define("abs_geneta", "abs(geneta)")
    .Define("en", def1)
    .Define("noise", def2)
    .Define("en_layer", def3)
    .Foreach(standard_ed, {"genen", "abs_geneta", "en", "noise", "en_layer"});

  for(int_ ireg=0; ireg<nreg; ++ireg) {
    for(int_ il=0; il<nlayers; ++il) 
      frac_en[ireg][il] /= static_cast<float_>(countfrac_en[ireg][il]);
  }

  auto get_weights = [&](float_ gen, float_ geta, 
			 vec1d<float_> en, vec1d<float_> noi,
			 vec2d<float_> en_layer) {
    for(int_ ireg=1; ireg<=nreg; ++ireg) {
      float_ f1=1., f2=0.;
      typename vec1d<float_>::const_iterator it;
      for(it = etareg.cbegin(); it!=(etareg.cend()-1); ++it) {
	int_ idx = it-etareg.cbegin();
	std::string idstr;
	//in case it lies outside the limits of the calibration
	//the event is calibrated with the full calibration region
	if (geta<etareg[0] || geta>etareg[n-1]) {
	  idstr = "sr" + std::to_string(ireg) + "from" + 
	    etastr(std::to_string(etareg[0])) + "to" + etastr(std::to_string(etareg[n-1]));
	}
	else if (geta<*it || geta>*(it+1))
	  continue;
	else {
	  idstr = "sr" + std::to_string(ireg) + "from" + 
	    etastr(std::to_string(*it)) + "to" + etastr(std::to_string(*(it+1)));
	}

	//it should only get here once per event
	if (calib[0].size() > 0) {
	  check_key(calib[0], idstr);
	  f1 /= calib[0][idstr]->Eval(geta)+1.0;
	  if (calib[1].size() > 0) {
	    check_key(calib[1], idstr);
	    f1 /= calib[1][idstr]->Eval(f1*en[ireg])+1.0;
	    if (calib[2].size() > 0) {
	      check_key(calib[2], idstr);
	      f2 = calib[2][idstr]->Eval(noi[2]);
	    }
	  }
	}
	en[ireg-1] = f1*en[ireg-1] - f2;
      } 
      assert(f1!=1.);

      float_ deltaE = en[ireg-1]/gen-1.;
      std::array<float_, nlayers> ROI_en = {0.}; 
      for(int_ il=0; il<nlayers; ++il) {
	float_ v = f1*en_layer[ireg-1][il] - f2;
	if(en[ireg-1]==0)
	  ROI_en[il] = 0.;
	else
	  ROI_en[il] = v/en[ireg-1];
      }

      vec1d<float_> lshift;
      if(samples=="inner")
	lshift = {.65, .59, .48};
      else if(samples=="outer")
	lshift = {1., 1., 1.};
      
    }

  };
  

  //connect two trees
  /*TFile* f1 = new TFile(noPUFile.c_str());
  TTree* t1 = static_cast<TTree*>(f1->Get("data"));
  TFile* f2 = new TFile("calibrated_en.root");
  TTree* t2 = static_cast<TTree*>(f2->Get("calibrated_en"));
  t1->AddFriend(t2, "friend");
  */

  ROOT::EnableImplicitMT(10);
  ROOT::RDataFrame d2("data", noPUFile.c_str());
  d2.Define("abs_geneta", "abs(geneta)")
    .Define("en", def1)
    .Define("noise", def2)
    .Define("en_layer", def3)
    .Foreach(get_weights, {"genen", "abs_geneta", "en", "noise", "en_layer"});
  

  /*delete f1;
  delete f2;
  */
  return 0;
}
