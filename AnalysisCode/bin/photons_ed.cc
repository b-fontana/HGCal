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

  /*initialize the master class for the analysis with photons.
    The calibration functions are calculated automatically.*/
  PhotonsPartialWafersAnalysis ana(argv[2], argv[4], argv[6]);
  CalibratorInputParameters p = ana.get_input_parameters();
  std::string slimmedtree = "slimmedT", slimmedfile = "slimmedF";

  /*create the standard energy distribution profile that will be used as the reference
    for all the events to decide whether they are signal or, if not, in which
    background category they belong.*/
  ROOT::EnableImplicitMT(std::thread::hardware_concurrency());
  ROOT::RDataFrame d1("data", p.noPUFile.c_str());
  auto df = d1.Define("abs_geneta", "abs(geneta)")
    .Define("en", ana.define_en)
    .Define("noise", ana.define_noise)
    .Define("en_layer", ana.define_en_layer)
    .Define("f1", ana.define_f1(), {"abs_geneta", "en"})
    .Define("f2", ana.define_f2(), {"abs_geneta", "noise"})
    .Define("en_calib", ana.define_calibrated_energy(), {"en", "f1", "f2"})
    .Define("en_layer_calib", ana.define_calibrated_energy_per_layer(), {"en_layer", "f1", "f2"});
  df.Snapshot(slimmedtree, slimmedfile+".root", {"genen", "abs_geneta", "genphi", "en_calib", "en_layer_calib"});
  auto compshowers_df = df.Filter(ana.complete_showers_filter(), {"abs_geneta", "en_calib"})
    .Define("frac_en", ana.complete_showers_energy_profile(), {"en_layer_calib", "en_calib", "f1", "f2"});
  
  using RDFptr = ROOT::RDF::RResultPtr<double>;
  vec2d<RDFptr> frac_en_ptr(p.nreg, vec1d<RDFptr>(p.nlayers));
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

  /*create the weights required for the "shower leakage" correction,
    which is based on the individual shower energy distribution.*/
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

  //produces the weights that will be later used to correct the energy distributions
  vector3d<int_> histcount(p.nreg, p.bckgcuts.size()+1, p.nlayers);
  std::mutex mut;
  auto get_weights = [&](vec2d<float_>& en_profile, const vec1d<int>& shower_idx) {
    for(uint_ ireg=1; ireg<=p.nreg; ++ireg) 
      {
	vec1d<float_> lshift;
	if(p.samples=="inner")
	  {
	    lshift = {.65,.59,.48};
	  }
	else if(p.samples=="outer")
	  lshift = {1., 1., 1.};
	assert(p.bckgcuts.size()==lshift.size());

	//calculate and calibrate the energy per layer
	for(uint_ il=1; il<=p.nlayers; ++il) {
	  std::string sreg = std::to_string(ireg);
	  int_ bin = hist["en"+sreg+"_layers_signal"]->FindBin(il);
	  if(shower_idx[ireg-1]==0) {
	    {
	      std::lock_guard lock(mut);
	      hist["en"+sreg+"_layers_signal"]->Fill(bin, en_profile[ireg-1][il-1]);
	      histcount(ireg-1,0,il-1) += 1;
	    }
	  }
	  else if(shower_idx[ireg-1]>0 && shower_idx[ireg-1]<=static_cast<int>(p.bckgcuts.size())) {
	    std::string sbckg = std::to_string(shower_idx[ireg-1]);
	    {
	      std::lock_guard lock(mut);
	      int_ bin_shift = hist["en"+sreg+"_layers_signal"]->FindBin(il*lshift[shower_idx[ireg-1]-1]);	     
	      hist["en"+sreg+"_layers_bckg"+sbckg]->Fill(bin_shift, en_profile[ireg-1][il-1]);
	    }
	    if(std::round(il*lshift[shower_idx[ireg-1]-1])==0) //avoid going for layer -1
	      {
		std::lock_guard lock(mut);
		histcount(ireg-1,shower_idx[ireg-1],0) += 1;
	      }
	    else
	      {
		std::lock_guard lock(mut);
		histcount.at(ireg-1,shower_idx[ireg-1],std::round(il*lshift[shower_idx[ireg-1]-1])-1) += 1;
	      }
	  }
	  else {
	    std::cout << "The showerid has some problem." << std::endl;
	    std::exit(0);
	  }
	}
      }
  };

  ROOT::RDataFrame d2(slimmedtree, slimmedfile+".root");
  auto define_df2 = d2.Define("en_profile", ana.define_shower_energy_profile(), {"en_calib", "en_layer_calib"})
    .Define("sid", ana.define_shower_index(frac_en), {"en_calib", "en_layer_calib", "en_profile"});
  define_df2.Snapshot(slimmedtree+"2", slimmedfile+"2.root", {"genen", "abs_geneta", "genphi", "en_calib", "en_layer_calib", "sid"});
  define_df2.Foreach(get_weights, {"en_profile", "sid"});
  
  std::string fwname = "root_files/fileweights_"+std::to_string(p.mask)+std::string(p.samples)+".root";
  TFile *fileweights = new TFile(fwname.c_str(),"RECREATE");
  fileweights->cd();
  for(uint_ ireg=1; ireg<=p.nreg; ++ireg) {
    //normalize
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

  /*perform the shower leakage correction and calculate uncorrected and corrected response*/
  SoftwareCorrector showercorr("root_files/fileweights_"+std::to_string(p.mask)+p.samples+".root");  
  std::string fname = "root_files/final_" + std::to_string(p.mask) + std::string(p.samples)+"_"+std::string(p.method)+".root";
  ROOT::RDataFrame dfinal(slimmedtree+"2", slimmedfile+"2.root");
  dfinal.Define("deltaE", ana.calculate_response(), {"genen", "en_calib"})
    .Define("en_calib_corr", ana.calculate_weight_corrected_energy(showercorr), {"en_calib", "en_layer_calib", "sid"})
    .Define("en_calib_corr2", ana.apply_low_stats_factor(showercorr), {"en_calib_corr", "sid"})
    .Define("deltaE_corr", ana.calculate_response(), {"genen", "en_calib_corr2"})
    .Snapshot("data", fname, {"genen", "abs_geneta", "genphi", "deltaE", "deltaE_corr", "sid"});
  
  return 0;
}
