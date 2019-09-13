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
#include "UserCode/AnalysisCode/interface/photon_analysis.h"
#include "UserCode/AnalysisCode/interface/software_correction.h"
#include "UserCode/AnalysisCode/interface/parser.h"

int_ main(int_ argc, char_ **argv) {
  if(argc!=7) 
    {
      std::cout << "Please specify the mask, samples and method to be used." << std::endl;
      std::exit(0);
    }
  if(argv[1] != std::string("--mask"))
    {
      std::cout << "The first argument must specify the mask to be used." << std::endl;
      std::exit(0);
    }
  if(argv[3] != std::string("--samples"))
    {
      std::cout << "The second argument must specify the samples to be used." << std::endl;
      std::exit(0);
    }
  if(argv[5] != std::string("--method"))
    {
      std::cout << "The third argument must specify the samples to be used." << std::endl;
      std::exit(0);
    }

  PhotonsPartialWafersAnalysis ana(argv[2], argv[4], argv[6]);
  ana.create_photon_calibration_values(6, false, true);
  vec1d<mapstr<TF1*>> calib = ana.calibration_values;
  CalibratorInputParameters p = ana.get_input_parameters();
  size_t n = p.etareg.size();
  std::string slimmedtree = "slimmedT", slimmedfile = "slimmedF.root";

  //////////////////////////////////////////////////
  //////////////FIRST PART//////////////////////////
  //////////////////////////////////////////////////  
  ROOT::RDataFrame d1("data", p.noPUFile.c_str());
  auto df = d1.Define("abs_geneta", "abs(geneta)")
    .Define("en", ana.define_en)
    .Define("noise", ana.define_noise)
    .Define("en_layer", ana.define_en_layer)
    .Define("f1", ana.define_f1(), {"abs_geneta", "en"})
    .Define("f2", ana.define_f2(), {"abs_geneta", "noise"})
    .Define("en_calib", ana.define_calibrated_energy(), {"en", "f1", "f2"})
    .Define("en_layer_calib", ana.define_calibrated_energy_per_layer(), {"en_layer", "f1", "f2"});
  df.Snapshot(slimmedtree, slimmedfile, {"genen", "abs_geneta", "en_calib", "en_layer_calib"});
  auto compshowers_df = df.Filter(ana.complete_showers_filter(), {"abs_geneta", "en_calib"})
    .Define("frac_en", ana.complete_showers_energy_profile(), {"en_layer_calib", "en_calib", "f1", "f2"});
  
  vec2d<ROOT::RDF::RResultPtr<double>> frac_en_ptr(p.nreg, vec1d<ROOT::RDF::RResultPtr<double>>(p.nlayers));
  for(uint_ ireg=0; ireg<p.nreg; ++ireg) 
    {
      for(uint_ il=0; il<p.nlayers; ++il) 
	{
	  std::string varname = "frac_en_"+std::to_string(ireg)+"_"+std::to_string(il);
	  std::string varsintax = "frac_en["+std::to_string(ireg)+"]["+std::to_string(il)+"]";
	  frac_en_ptr[ireg][il] = compshowers_df.Define(varname, varsintax).Mean<float_>(varname); //lazy
	}
    }
  vec2d<float_> frac_en(p.nreg, vec1d<float_>(p.nlayers, 0.));
  for(uint_ ireg=0; ireg<p.nreg; ++ireg) 
    for(uint_ il=0; il<p.nlayers; ++il)   
      frac_en[ireg][il] = *frac_en_ptr[ireg][il];

  //////////////////////////////////////////////////
  //////////////SECOND PART/////////////////////////
  //////////////////////////////////////////////////
  //create histograms to store the average weights
  mapstr<TH1F*> hist;
  float_ bins[p.nlayers+1];
  for(uint_ ib=0; ib<=p.nlayers; ++ib)
    bins[ib] = 0.5 + ib;
  for(uint_ ireg=0; ireg<p.nreg; ++ireg) 
    {
      std::string s = "en"+std::to_string(ireg+1)+"_layers_signal";
      TH1F* htmp = new TH1F(s.c_str(), ";Layer;E_{reco}/E_{gen}", p.nlayers, bins);
      htmp->SetMarkerStyle(20);
      htmp->SetDirectory(0);
      hist.insert({{s, htmp}});
      for(uint_ iw=0; iw<p.bckgcuts.size(); ++iw) 
	{
	  std::string s = "en"+std::to_string(ireg+1)+"_layers_bckg"+std::to_string(iw+1);
	  TH1F* htmp = new TH1F(s.c_str(), ";Layer;E_{reco}/E_{gen}", 
				p.nlayers, bins);
	  htmp->SetMarkerStyle(20);
	  htmp->SetDirectory(0);
	  hist.insert({{s, htmp}});
	}
    }
  vector3d<int_> histcount(p.nreg, p.bckgcuts.size()+1, p.nlayers);

  std::exit(0);
  //produces the weights that will be later used to correct the energy distributions
  auto get_weights = [&](vec1d<float_> encalib, vec2d<float_> encalib_layer) {
    int sid=-1;
    for(uint_ ireg=1; ireg<=p.nreg; ++ireg) 
      {
	vec1d<float_> ROI_en(p.nlayers, 0.); 
	for(uint_ il=0; il<p.nlayers; ++il) 
	  {
	    if(encalib[ireg-1]==0)
	      ROI_en[il] = 0.;
	    else
	      ROI_en[il] = encalib_layer[ireg-1][il]/encalib[ireg-1];
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

  ROOT::RDataFrame d2("newtree", "newfile.root");
  d2.Foreach(get_weights, {"en_calib", "en_layer"});
  
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

  ROOT::RDataFrame dfinal("data", p.noPUFile.c_str());
  dfinal.Define("abs_geneta", "fabs(geneta)")
    .Define("en", ana.define_en)
    .Define("noise", ana.define_noise)
    .Define("en_layer", ana.define_en_layer)
    .Foreach(apply_weights, {"genen", "abs_geneta", "genphi", "en", "noise", "en_layer"});

  file->cd();
  tree->Write();
  file->Close();
  delete file;
  
  return 0;
}
