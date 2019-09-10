#include <iostream>
#include <vector>
#include <iterator>
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "ROOT/RDataFrame.hxx"
#include "TH1F.h"

#include "UserCode/AnalysisCode/interface/utils.h"
#include "UserCode/AnalysisCode/interface/calibration.h"
#include "UserCode/AnalysisCode/interface/software_correction.h"
#include "UserCode/AnalysisCode/interface/parser.h"

int_ main(int_ argc, char_ **argv) {
  if(argc!=5) 
    {
      std::cout << "Please specify both the samples and the mask to be used." << std::endl;
      std::exit(0);
    }
  if(argv[3] != std::string("--samples"))
    {
      std::cout << "The first argument must specify the samples to be used." << std::endl;
      std::exit(0);
    }
  else if(argv[1] != std::string("--mask"))
    {
      std::cout << "The second argument must specify the mask to be used." << std::endl;
      std::exit(0);
    }

  //Calibration
  vec1d<std::string> varnames = {"mingenen", "etareg_"+std::string(argv[4]), "nreg", "nlayers", 
				 "input", "bckgcuts_"+std::string(argv[4])};
  CalibratorInputParameters p("params_pions.csv", varnames, std::string(argv[2]), std::string(argv[4]));
  Calibrator calibrator(p);
  calibrator.create_photon_calibration_values(6, false, false);
  vec1d<mapstr<TF1*>> calib = calibrator.calibration_values;

  //variables
  uint_ ncores = std::thread::hardware_concurrency();
  size_t n = p.etareg.size();

  //////////////////////////////////////////////////
  //////////////FIRST PART//////////////////////////
  //////////////////////////////////////////////////

  //vectors to create standard energy distribution for complete showers
  vec2d<float_> frac_en(p.nreg, vec1d<float_>(p.nlayers, 0.));
  vec2d<int_> countfrac_en(p.nreg, vec1d<int_>(p.nlayers, 0));

  //fills the 'frac_en' and 'countfrac_en' vectors
  auto standard_ed = [&](float_ gen, float_ geta, 
			 vec1d<float_> en, vec1d<float_> noi,
			 vec2d<float_>  en_layer) {
    float_ f1 = 1., f2 = 0.;
    for(uint_ ireg=1; ireg<=p.nreg; ++ireg) {
      f1 = 1.;
      f2 = 0.;
      float_ encalib = 0.;
      typename vec1d<float_>::const_iterator it;
      for(it = p.etareg.cbegin(); it!=(p.etareg.cend()-1); ++it) {
	std::string idstr;
	//in case it lies outside the limits of the calibration
	//the event is calibrated with the full calibration region
	if (geta<p.etareg[0] || geta>p.etareg[n-1]) {
	  idstr = "sr" + std::to_string(ireg) + "from" + 
	    etastr(std::to_string(p.etareg[0])) + "to" + etastr(std::to_string(p.etareg[n-1]));
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
	  f1 /= calib[0][idstr]->Eval(geta)+1.0;
	  if (calib[1].size() > 0) {
	    check_key(idstr, calib[1]);
	    f1 /= calib[1][idstr]->Eval(f1*en[ireg-1])+1.0;
	    if (calib[2].size() > 0) {
	      check_key(idstr, calib[2]);
	      f2 = calib[2][idstr]->Eval(noi[2]);
	    }
	  }
	}
	encalib = f1*en[ireg-1] - f2;
      }
      assert(f1!=1.);
      for(uint_ il=1; il<=p.nlayers; ++il) {
	float_ v = f1*en_layer[ireg-1][il-1] - f2;
	if( ((p.samples=="inner" && geta<p.etareg[0]+0.05) || (p.samples=="outer" && geta>p.etareg[n-1]-0.05))
	    && encalib != 0. ) {
	  frac_en[ireg-1][il-1] += (v/encalib);
	  countfrac_en[ireg-1][il-1] += 1;
	}
      }
    }
  };

  //create histograms to store the average weights
  mapstr<TH1F*> hist;
  float_ bins[p.nlayers+1];
  for(uint_ ib=0; ib<=p.nlayers; ++ib)
    bins[ib] = 0.5 + ib;
  for(uint_ ireg=0; ireg<p.nreg; ++ireg) {
    std::string s = "en"+std::to_string(ireg+1)+"_layers_signal";
    TH1F* htmp = new TH1F(s.c_str(), ";Layer;E_{reco}/E_{gen}", 
			       p.nlayers, bins);
    htmp->SetMarkerStyle(20);
    htmp->SetDirectory(0);
    hist.insert({{s, htmp}});
    for(uint_ iw=0; iw<p.bckgcuts.size(); ++iw) {
      std::string s = "en"+std::to_string(ireg+1)+"_layers_bckg"+std::to_string(iw+1);
      TH1F* htmp = new TH1F(s.c_str(), ";Layer;E_{reco}/E_{gen}", 
				 p.nlayers, bins);
      htmp->SetMarkerStyle(20);
      htmp->SetDirectory(0);
      hist.insert({{s, htmp}});
    }
  }
  vector3d<int_> histcount(p.nreg, p.bckgcuts.size()+1, p.nlayers);

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

  ROOT::EnableImplicitMT(ncores);
  ROOT::RDataFrame d1("data", p.noPUFile.c_str());
  d1.Define("abs_geneta", "abs(geneta)")
    .Define("en", def1)
    .Define("noise", def2)
    .Define("en_layer", def3)
    .Foreach(standard_ed, {"genen", "abs_geneta", "en", "noise", "en_layer"});

  for(uint_ ireg=0; ireg<p.nreg; ++ireg) {
    for(uint_ il=0; il<p.nlayers; ++il) 
      frac_en[ireg][il] /= static_cast<float_>(countfrac_en[ireg][il]);
  }

  //////////////////////////////////////////////////
  //////////////SECOND PART/////////////////////////
  //////////////////////////////////////////////////

  //produces the weights that will be later used to correct the energy distributions
  auto get_weights = [&](float_ gen, float_ geta, 
			 vec1d<float_> en, vec1d<float_> noi,
			 vec2d<float_> en_layer) {
    float f1=1., f2=0.;
    int sid=-1;
    for(uint_ ireg=1; ireg<=p.nreg; ++ireg) {
      f1 = 1.;
      f2 = 0.;
      float_ encalib = 0.;
      typename vec1d<float_>::const_iterator it;
      for(it = p.etareg.cbegin(); it!=(p.etareg.cend()-1); ++it) {
	std::string idstr;
	//in case it lies outside the limits of the calibration
	//the event is calibrated with the full calibration region
	if (geta<p.etareg[0] || geta>p.etareg[n-1]) {
	  idstr = "sr" + std::to_string(ireg) + "from" + 
	    etastr(std::to_string(p.etareg[0])) + "to" + etastr(std::to_string(p.etareg[n-1]));
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
	  f1 /= calib[0][idstr]->Eval(geta)+1.0;
	  if (calib[1].size() > 0) {
	    check_key(idstr, calib[1]);
	    f1 /= calib[1][idstr]->Eval(f1*en[ireg-1])+1.0;
	    if (calib[2].size() > 0) {
	      check_key(idstr, calib[2]);
	      f2 = calib[2][idstr]->Eval(noi[2]);
	    }
	  }
	}
	encalib = f1*en[ireg-1] - f2;
      }
      assert(f1!=1.);

      vec1d<float_> ROI_en(p.nlayers, 0.); 
      for(uint_ il=0; il<p.nlayers; ++il) {
	float_ v = f1*en_layer[ireg-1][il] - f2;
	if(encalib==0)
	  ROI_en[il] = 0.;
	else
	  ROI_en[il] = v/encalib;
      }

      vec1d<float_> lshift;
      if(p.samples=="inner")
	lshift = {.65,.59,.48};
      else if(p.samples=="outer")
	lshift = {1., 1., 1.};
      assert(p.bckgcuts.size()==lshift.size());
      sid = diff_ed(ROI_en, frac_en[ireg-1], p.bckgcuts, 0.05);

      //calculate and calibrate the energy per layer
      for(uint_ il=1; il<=p.nlayers; ++il) {
	std::string sreg = std::to_string(ireg);
	int_ bin = hist["en"+sreg+"_layers_signal"]->FindBin(il);
	if(sid==0) {
	  hist["en"+sreg+"_layers_signal"]->Fill(bin, ROI_en[il-1]);
	  histcount(ireg-1,0,il-1) += 1;
	}
	else if(sid>0 && sid<=static_cast<int>(p.bckgcuts.size())) {
	  std::string sbckg = std::to_string(sid);
	  int_ bin_shift = hist["en"+sreg+"_layers_signal"]->FindBin(il*lshift[sid-1]);
	  hist["en"+sreg+"_layers_bckg"+sbckg]->Fill(bin_shift, ROI_en[il-1]);
	  if(std::round(il*lshift[sid-1])==0) //avoid going for layer -1
	    histcount(ireg-1,sid,0) += 1;
	  else
	    histcount(ireg-1,sid,std::round(il*lshift[sid-1])-1) += 1;
	}
	else {
	  std::cout << "The showerid has some problem." << std::endl;
	  std::exit(0);
	}
      }
    }
  };

  ROOT::EnableImplicitMT(ncores);
  ROOT::RDataFrame d2("data", p.noPUFile.c_str());
  auto d2_def = d2.Define("abs_geneta", "fabs(geneta)")
    .Define("en", def1)
    .Define("noise", def2)
    .Define("en_layer", def3);
  d2_def.Foreach(get_weights, {"genen", "abs_geneta", "en", "noise", "en_layer"});
  
  std::string fwname = "root_files/fileweights_"+std::to_string(p.mask)+std::string(p.samples)+".root";
  TFile *fileweights = new TFile(fwname.c_str(),"RECREATE");
  fileweights->cd();
  for(uint_ ireg=1; ireg<=p.nreg; ++ireg) {
    //normalize (TProfile does not work with RDataFrame)
    std::string str_s = "en"+std::to_string(ireg)+"_layers_signal";
    for(uint_ il=1; il<=p.nlayers; ++il) {
      int_ bin = hist[str_s]->FindBin(il);
      hist[str_s]->SetBinContent(bin, hist[str_s]->GetBinContent(bin) / histcount(ireg-1,0,il-1));
    }
    for(uint_ iw=1; iw<=p.bckgcuts.size(); ++iw) {
      std::string str_b = "en"+std::to_string(ireg)+"_layers_bckg"+std::to_string(iw);
      for(uint_ il=1; il<=p.nlayers; ++il) {
	int_ bin = hist[str_b]->FindBin(il);
	hist[str_b]->SetBinContent(bin, hist[str_b]->GetBinContent(bin) / histcount(ireg-1,iw,il-1));
      }
    }
    //write to file
    TH1F* h_s = static_cast<TH1F*>(hist[str_s]->Clone(("en"+std::to_string(ireg)+"_layer_sign").c_str()));
    h_s->Write();
    for(uint_ iw=1; iw<=p.bckgcuts.size(); ++iw) {
      std::string str_b = "en"+std::to_string(ireg)+"_layers_bckg"+std::to_string(iw);
      std::string sclone = "weight"+std::to_string(iw)+"_sr"+std::to_string(ireg);
      hist[str_b]->Write();
      TH1F* h_b = static_cast<TH1F*>(hist[str_b]->Clone( sclone.c_str() ));
      h_b->Divide(hist[str_s]);
      h_b->Write();
    }
  }
  fileweights->Close();
  delete fileweights;
  for(auto const& v: hist) delete v.second;

  //////////////////////////////////////////////////
  //////////////THIRD PART//////////////////////////
  //////////////////////////////////////////////////

  std::string fawname = "root_files/final_" + std::to_string(p.mask) + std::string(p.samples)+"_ed.root";
  TFile *file = new TFile(fawname.c_str(), "RECREATE");
  file->cd();
  TTree *tree = new TTree("data", "Tree after weights");
  float_ geneta, genphi;
  vec1d<float_> deltaE(p.nreg), deltaE_corr(p.nreg);
  vec1d<float_> shower_number(p.nreg);
  tree->Branch("geneta", &geneta);
  tree->Branch("genphi", &genphi);
  tree->Branch("deltaE", &deltaE);
  tree->Branch("deltaE_corr", &deltaE_corr);
  tree->Branch("shower_number", &shower_number);

  SoftwareCorrection showercorr("root_files/fileweights_"+std::to_string(p.mask)+p.samples+".root");
  vector3d<float_> weights = showercorr.weights;
  vec1d<uint_> boundaries(p.nreg);
  vec1d<float_> lshift(p.nreg);
  std::string corr_mode;
  if(p.samples=="inner") {
    boundaries = {5, 5, 5};
    corr_mode = "left";
    lshift = {.65, .59, .48};
  }
  else if(p.samples=="outer") {
    boundaries = {23, 23, 23};
    corr_mode = "right";
    lshift = {1., 1., 1.};
  }
  vec1d<float_> lowfact = showercorr.low_stats_factor(boundaries, corr_mode);

  auto apply_weights = [&](float_ gen, float_ geta, float_ gphi,
			   vec1d<float_> en, vec1d<float_> noi,
			   vec2d<float_> en_layer) {
    float_ f1=1., f2=0.;
    int_ sid=-1;
    for(uint_ ireg=0; ireg<p.nreg; ++ireg) {
      f1 = 1.;
      f2 = 0.;
      float_ encalib = 0.;
      typename vec1d<float_>::const_iterator it;
      for(it = p.etareg.cbegin(); it!=(p.etareg.cend()-1); ++it) {
	std::string idstr;
	//in case it lies outside the limits of the calibration
	//the event is calibrated with the full calibration region
	if (geta<p.etareg[0] || geta>p.etareg[n-1]) {
	  idstr = "sr" + std::to_string(ireg+1) + "from" + 
	    etastr(std::to_string(p.etareg[0])) + "to" + etastr(std::to_string(p.etareg[n-1]));
	}
	else if (geta<*it || geta>*(it+1))
	  continue;
	else {
	  idstr = "sr" + std::to_string(ireg+1) + "from" + 
	    etastr(std::to_string(*it)) + "to" + etastr(std::to_string(*(it+1)));
	}

	//it should only get here once per event
	if (calib[0].size() > 0) {
	  check_key(idstr, calib[0]);
	  f1 /= calib[0][idstr]->Eval(geta)+1.0;
	  if (calib[1].size() > 0) {
	    check_key(idstr, calib[1]);
	    f1 /= calib[1][idstr]->Eval(f1*en[ireg])+1.0;
	    if (calib[2].size() > 0) {
	      check_key(idstr, calib[2]);
	      f2 = calib[2][idstr]->Eval(noi[2]);
	    }
	  }
	}
	encalib = f1*en[ireg] - f2;
      }
      assert(f1!=1.);

      vec1d<float_> ROI_en(p.nlayers, 0.); 
      for(uint_ il=0; il<p.nlayers; ++il) {
	float_ v = f1*en_layer[ireg][il] - f2;
	if(encalib==0)
	  ROI_en[il] = 0.;
	else
	  ROI_en[il] = v/encalib;
      }

      vec1d<float_> lshift;
      if(p.samples=="inner")
	lshift = {.65,.59,.48};
      else if(p.samples=="outer")
	lshift = {1., 1., 1.};
      assert(p.bckgcuts.size()==lshift.size());
      sid = diff_ed(ROI_en, frac_en[ireg], p.bckgcuts, 0.05);

      bool_ weight_limit;
      float_ encorr = 0.;
      for (uint_ il=1; il<=p.nlayers; ++il) {
	if(p.samples=="inner")
	  weight_limit = il > boundaries[ireg];
	else if(p.samples=="outer")
	  weight_limit = il < boundaries[ireg];
	float_ v = f1*en_layer.at(ireg).at(il-1) - f2;
	assert(f2==0.);
	assert(f1!=1.);  
	if(sid==0) {
	  encorr += v;
	}
	else if(sid>0 && sid<=static_cast<int_>(p.bckgcuts.size())) {
	  int_ w = sid-1;
	  if(weights(ireg,w,il-1) != 0 && weight_limit) {
	    int_ r = static_cast<int_>(std::round(il*lshift[w])-1);
	    if(r==-1) 
	      r = 0;
	    encorr += v/weights(ireg,w,r);
	  }
	}
	else {
	  std::cout << "There is a problem with the showerid." << std::endl;
	  std::exit(0);
	}
      }

      deltaE.at(ireg) = encalib/gen - 1.;
      if(sid>0) {
	encorr *= ( 1 / (1-lowfact[ireg]) );
	if(p.samples=="inner") encorr *= 1/0.09;
	else if(p.samples=="outer") encorr *= 1/0.08;
      }
      deltaE_corr.at(ireg) = encorr/gen - 1.;
      shower_number.at(ireg) = sid;
    }
    geneta = geta;
    genphi = gphi;
    tree->Fill();
  };

  ROOT::EnableImplicitMT(ncores);
  ROOT::RDataFrame dfinal("data", p.noPUFile.c_str());
  dfinal.Define("abs_geneta", "fabs(geneta)")
    .Define("en", def1)
    .Define("noise", def2)
    .Define("en_layer", def3)
    .Foreach(apply_weights, {"genen", "abs_geneta", "genphi", "en", "noise", "en_layer"});

  file->cd();
  tree->Write();
  file->Close();
  delete file;

  return 0;
}
