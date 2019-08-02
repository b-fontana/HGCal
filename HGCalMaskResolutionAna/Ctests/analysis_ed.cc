#include <iostream>
#include <vector>
#include <iterator>
#include "TCanvas.h"
#include "TFile.h"
#include "TObject.h"
#include "TTree.h"
#include "ROOT/RDataFrame.hxx"
#include "TProfile.h"

#include "interface/utils.h"
#include "interface/calibration.h"
#include "interface/software_correction.h"

int_ main() {
  unsigned int ncores = std::thread::hardware_concurrency();

  //variables
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
  vec1d<float_> bckgcuts = {0.5, 0.6, 0.7};

  //Calibration
  Calibration calibration(mingenen, etareg, mask, label, samples, 
			  mask, noPUFile, outpath);
  calibration.nopu_calibration(6, true);
  vec1d<mapstr<TF1*>> calib = calibration.calib;


  //////////////////////////////////////////////////
  //////////////FIRST PART//////////////////////////
  //////////////////////////////////////////////////

  //create TTree that will used as a friend of the TTree being read by the RDataFrame
  //the goal is to avoid redoing some calculations in consecutive data loops
  vec1d<float> f1(nreg, 1.);
  vec1d<float> f2(nreg, 0.);
  TFile* filefriend = new TFile("filefriend.root", "RECREATE");
  filefriend->cd();
  TTree *tfriend1 = new TTree("tfriend1", "First Friend Tree");
  tfriend1->Branch("f1", &f1);
  tfriend1->Branch("f2", &f2);

  //vectors to create standard energy distribution for complete showers
  vec2d<float_> frac_en(nreg, vec1d<float_>(nlayers, 0.));
  vec2d<int_> countfrac_en(nreg, vec1d<int_>(nlayers, 0));
  
  //fills the 'frac_en' and 'countfrac_en' vectors
  auto standard_ed = [&](float_ gen, float_ geta, 
			 vec1d<float_> en, vec1d<float_> noi,
			 vec2d<float_>  en_layer) {
    for(int_ ireg=1; ireg<=nreg; ++ireg) {
      f1.at(ireg-1) = 1.;
      f2.at(ireg-1) = 0.;
      float_ encorr = 0.;
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
	  check_key(idstr, calib[0]);
	  f1.at(ireg-1) /= calib[0][idstr]->Eval(geta)+1.0;
	  if (calib[1].size() > 0) {
	    check_key(idstr, calib[1]);
	    f1.at(ireg-1) /= calib[1][idstr]->Eval(f1.at(ireg-1)*en[ireg-1])+1.0;
	    if (calib[2].size() > 0) {
	      check_key(idstr, calib[2]);
	      f2.at(ireg-1) = calib[2][idstr]->Eval(noi[2]);
	    }
	  }
	}
	encorr = f1.at(ireg-1)*en[ireg-1] - f2.at(ireg-1);
      }
      assert(f1.at(ireg-1)!=1.);
      for(int_ il=1; il<=nlayers; ++il) {
	float_ v = f1.at(ireg-1)*en_layer[ireg-1][il-1] - f2.at(ireg-1);
	if( (samples=="inner" && geta<etareg[0]+0.05 || samples=="outer" and geta>1.6)
	    && encorr != 0. ) {
	  frac_en[ireg-1][il-1] += (v/encorr);
	  countfrac_en[ireg-1][il-1] += 1;
	}
      }
    }
    tfriend1->Fill();
  };

  //create histograms to store the average weights (profiles)
  mapstr<TProfile*> hist;
  float_ bins[nlayers+1];
  for(int_ ib=0; ib<nlayers+1; ++ib)
    bins[ib] = 0.5 + ib;
  for(int_ ireg=0; ireg<nreg; ++ireg) {
    std::string s = "en"+std::to_string(ireg+1)+"_layers_signal";
    TProfile* p = new TProfile(s.c_str(), ";Layer;E_{reco}/E_{gen}", 
			       nlayers, bins);
    p->SetMarkerStyle(20);
    p->SetDirectory(0);
    hist.insert({{s, p}});
  }
  for(int_ ireg=0; ireg<nreg; ++ireg) {
    for(int_ iw=0; iw<bckgcuts.size(); ++iw) {
      std::string s = "en"+std::to_string(ireg+1)+"_layers_bckg"+std::to_string(iw+1);
      TProfile* p = new TProfile(s.c_str(), ";Layer;E_{reco}/E_{gen}", 
				 nlayers, bins);
      p->SetMarkerStyle(20);
      p->SetDirectory(0);
      hist.insert({{s, p}});
    }
  }

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

  //ROOT::EnableImplicitMT(ncores);
  ROOT::RDataFrame d1("data", noPUFile.c_str());
  d1.Define("abs_geneta", "abs(geneta)")
    .Define("en", def1)
    .Define("noise", def2)
    .Define("en_layer", def3)
    .Foreach(standard_ed, {"genen", "abs_geneta", "en", "noise", "en_layer"});

  filefriend->cd();
  tfriend1->Write();

  for(int_ ireg=0; ireg<nreg; ++ireg) {
    for(int_ il=0; il<nlayers; ++il) 
      frac_en[ireg][il] /= static_cast<float_>(countfrac_en[ireg][il]);
  }

  std::cout << "END FIRST PART" << std::endl;
  //////////////////////////////////////////////////
  //////////////SECOND PART/////////////////////////
  //////////////////////////////////////////////////
  vec1d<int_> showerid(nreg, -1);
  TFile* filefriend2 = new TFile("filefriend2.root", "RECREATE");
  filefriend2->cd();
  TTree *tfriend2 = new TTree("tfriend2", "Second Friend Tree");
  tfriend2->Branch("showerid", &showerid);

  //produces the weights that will be later used to correct the energy distributions
  auto get_weights = [&](float_ gen, float_ geta, 
			 vec1d<float_> en, vec1d<float_> noi,
			 vec2d<float_> en_layer,
			 vec1d<float_> f1, vec1d<float_> f2) {
    filefriend->cd();
    for(int_ ireg=1; ireg<=nreg; ++ireg) {
      float encorr = 0.;
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
	encorr = f1.at(ireg-1)*en[ireg-1] - f2.at(ireg-1);
      }

      float_ deltaE = encorr/gen-1.;
      vec1d<float_> ROI_en(nlayers, 0.); 
      for(int_ il=0; il<nlayers; ++il) {
	float_ v = f1.at(ireg-1)*en_layer[ireg-1][il] - f2.at(ireg-1);
	if(encorr==0)
	  ROI_en[il] = 0.;
	else
	  ROI_en[il] = v/encorr;
      }

      vec1d<float_> lshift;
      if(samples=="inner")
	lshift = {.65, .59, .48};
      else if(samples=="outer")
	lshift = {1., 1., 1.};
      assert(bckgcuts.size()==lshift.size());
      showerid.at(ireg-1) = diff_ed(ROI_en, frac_en[ireg-1], bckgcuts, 0.05);

      //calculate and calibrate the energy per layer
      float_ en_corr = 0;
      for(int_ il=1; il<=nlayers; ++il) {
	std::string sreg = std::to_string(ireg);
	float_ bin = hist["en"+sreg+"_layers_signal"]->FindBin(il);
	if(showerid.at(ireg-1)==0) {
	  hist["en"+sreg+"_layers_signal"]->Fill(bin, ROI_en[il-1]);
	}
	else {
	  std::string sbckg = std::to_string(showerid.at(ireg-1));
	  hist["en"+sreg+"_layers_bckg"+sbckg]->Fill(bin*lshift[showerid.at(ireg-1)-1], ROI_en[il-1]);
	}
      }
    }
    tfriend2->Fill();
  };

  //join the original tree with the one made in this macro
  TFile* fileorig = new TFile(noPUFile.c_str(), "READ");
  fileorig->cd();
  TTree* treeorig = static_cast<TTree*>(fileorig->Get("data"));
  treeorig->AddFriend(tfriend1, "friend1");

  ROOT::EnableImplicitMT(ncores);
  ROOT::RDataFrame d2(*treeorig);
  auto d2_def = d2.Define("abs_geneta", "fabs(geneta)")
    .Define("en", def1)
    .Define("noise", def2)
    .Define("en_layer", def3);
  d2_def.Foreach(get_weights, {"genen", "abs_geneta", "en", "noise", "en_layer", "friend1.f1", "friend1.f2"});

  filefriend2->cd();
  tfriend2->Write();
  filefriend2->Close();
  delete filefriend2;
  fileorig->Close();
  delete fileorig;
  filefriend->Close();
  delete filefriend;  
  
  TFile *fileweights = new TFile("fileweights.root","RECREATE");
  fileweights->cd();
  for(int ireg=1; ireg<=nreg; ++ireg) {
    std::string str_s = "en"+std::to_string(ireg)+"_layers_signal";
    hist[str_s]->Clone(("en"+std::to_string(ireg)+"_layer_sign").c_str())->Write();
    for(int iw=1; iw<=bckgcuts.size(); ++iw) {
      std::string str_b = "en"+std::to_string(ireg)+"_layers_bckg"+std::to_string(iw);
      std::string sclone = "weight"+std::to_string(iw)+"_sr"+std::to_string(ireg);
      TH1F* h_b = static_cast<TH1F*>(hist[str_b]->Clone( sclone.c_str() ));
      h_b->Divide(hist[str_s]);
      h_b->Write();
    }
  }
  fileweights->Close();
  delete fileweights;
  for(auto const& [k,v]: hist) delete v;
  
  std::cout << "END SECOND PART" << std::endl;

  return 0;
}
