#include <iostream>
#include <vector>
#include <iterator>
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "ROOT/RDataFrame.hxx"
#include "TH1F.h"

#include /*"UserCode/AnalysisCode/*/"../interface/utils.h"
#include /*"UserCode/AnalysisCode/*/"../interface/calibration.h"
#include /*"UserCode/AnalysisCode/*/"../interface/software_correction.h"
#include /*"UserCode/AnalysisCode/*/"../interface/parser.h"

int_ main(int_ argc, char_ **argv) {
  if(argc!=5) 
    {
      std::cout << "Please specify both the samples and the mask to be used." << std::endl;
      std::exit(0);
    }
  if(argv[1] != std::string("--samples"))
    {
      std::cout << "The first argument must specify the samples to be used." << std::endl;
      std::exit(0);
    }
  else if(argv[3] != std::string("--mask"))
    {
      std::cout << "The second argument must specify the mask to be used." << std::endl;
      std::exit(0);
    }

  //variables
  std::string samples = argv[2];
  uint_ mask = std::stoi(argv[4]);
  uint_ ncores = std::thread::hardware_concurrency();

  float_ mingenen;
  vec1d<float_> etareg;
  int_ nreg;
  int_ nlayers;
  std::string label;
  std::string noPUFile;
  std::string outpath;
  vec1d<float_> bckgcuts;

  std::ifstream infile("params.csv");
  for(CSVIterator it(infile); it != CSVIterator(); ++it)
    {
      CSVRow row = *it;
      if(row[0]=="mingenen") 
	{
	  if( row.size()!=2 ) {row.bad_row();}
	  mingenen = std::stof(row[1]);
	}
      else if(row[0]=="etareg")
	{
	  if( row.size()!=3 ) {row.bad_row();}
	  etareg = {std::stof(row[1]), std::stof(row[2])};
	}
      else if(row[0]=="nreg") 
	{
	  if( row.size()!=2 ) {row.bad_row();}
	  nreg = std::stoi(row[1]);
	}
      else if(row[0]=="nlayers") 
	{
	  if( row.size()!=2 ) {row.bad_row();}
	  nlayers = std::stoi(row[1]);
	}
      else if(row[0]=="input") 
	{
	  if( row.size()!=2 ) {row.bad_row();}
	  noPUFile = row[1]+"mask"+std::string(argv[4])+"_"+std::string(argv[2])+".root";
	}
      else if(row[0]=="output") 
	{
	  if( row.size()!=2 ) {row.bad_row();}
	  outpath = row[1]+std::string(argv[2])+"/mask"+std::string(argv[4]);
	}
      else if(row[0]=="bckgcuts") 
	{
	  if( row.size()!=4 ) {row.bad_row();}
	  bckgcuts = {std:stof(row[1]),std:stof(row[2]),std:stof(row[3])};
	}
    }

  uint_ n = etareg.size();

  //Calibration
  Calibration calibration(mingenen, etareg, mask, label, samples, 
			  mask, noPUFile, outpath);
  calibration.nopu_calibration(6, true);
  vec1d<mapstr<TF1*>> calib = calibration.calib;


  //////////////////////////////////////////////////
  //////////////FIRST PART//////////////////////////
  //////////////////////////////////////////////////

  //vectors to create standard energy distribution for complete showers
  vec2d<float_> frac_en(nreg, vec1d<float_>(nlayers, 0.));
  vec2d<int_> countfrac_en(nreg, vec1d<int_>(nlayers, 0));

  //fills the 'frac_en' and 'countfrac_en' vectors
  auto standard_ed = [&](float_ gen, float_ geta, 
			 vec1d<float_> en, vec1d<float_> noi,
			 vec2d<float_>  en_layer) {
    float_ f1 = 1., f2 = 0.;
    for(int_ ireg=1; ireg<=nreg; ++ireg) {
      f1 = 1.;
      f2 = 0.;
      float_ encalib = 0.;
      typename vec1d<float_>::const_iterator it;
      for(it = etareg.cbegin(); it!=(etareg.cend()-1); ++it) {
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
      for(int_ il=1; il<=nlayers; ++il) {
	float_ v = f1*en_layer[ireg-1][il-1] - f2;
	if( ((samples=="inner" && geta<etareg[0]+0.05) || (samples=="outer" && geta>etareg[n-1]-0.05))
	    && encalib != 0. ) {
	  frac_en[ireg-1][il-1] += (v/encalib);
	  countfrac_en[ireg-1][il-1] += 1;
	}
      }
    }
  };

  //create histograms to store the average weights (histiles)
  mapstr<TH1F*> hist;
  float_ bins[nlayers+1];
  for(int_ ib=0; ib<=nlayers; ++ib)
    bins[ib] = 0.5 + ib;
  for(int_ ireg=0; ireg<nreg; ++ireg) {
    std::string s = "en"+std::to_string(ireg+1)+"_layers_signal";
    TH1F* p = new TH1F(s.c_str(), ";Layer;E_{reco}/E_{gen}", 
			       nlayers, bins);
    p->SetMarkerStyle(20);
    p->SetDirectory(0);
    hist.insert({{s, p}});
    for(uint_ iw=0; iw<bckgcuts.size(); ++iw) {
      std::string s = "en"+std::to_string(ireg+1)+"_layers_bckg"+std::to_string(iw+1);
      TH1F* p = new TH1F(s.c_str(), ";Layer;E_{reco}/E_{gen}", 
				 nlayers, bins);
      p->SetMarkerStyle(20);
      p->SetDirectory(0);
      hist.insert({{s, p}});
    }
  }
  vector3d<int_> histcount(nreg, bckgcuts.size()+1, nlayers);

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

  for(int_ ireg=0; ireg<nreg; ++ireg) {
    for(int_ il=0; il<nlayers; ++il) 
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
    for(int_ ireg=1; ireg<=nreg; ++ireg) {
      f1 = 1.;
      f2 = 0.;
      float_ encalib = 0.;
      typename vec1d<float_>::const_iterator it;
      for(it = etareg.cbegin(); it!=(etareg.cend()-1); ++it) {
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

      vec1d<float_> ROI_en(nlayers, 0.); 
      for(int_ il=0; il<nlayers; ++il) {
	float_ v = f1*en_layer[ireg-1][il] - f2;
	if(encalib==0)
	  ROI_en[il] = 0.;
	else
	  ROI_en[il] = v/encalib;
      }

      vec1d<float_> lshift;
      if(samples=="inner")
	lshift = {.65,.59,.48};
      else if(samples=="outer")
	lshift = {1., 1., 1.};
      assert(bckgcuts.size()==lshift.size());
      sid = diff_ed(ROI_en, frac_en[ireg-1], bckgcuts, 0.05);

      //calculate and calibrate the energy per layer
      for(int_ il=1; il<=nlayers; ++il) {
	std::string sreg = std::to_string(ireg);
	int_ bin = hist["en"+sreg+"_layers_signal"]->FindBin(il);
	if(sid==0) {
	  hist["en"+sreg+"_layers_signal"]->Fill(bin, ROI_en[il-1]);
	  histcount(ireg-1,0,il-1) += 1;
	}
	else if(sid>0 && sid<=static_cast<int>(bckgcuts.size())) {
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

  //ROOT::EnableImplicitMT(ncores);
  ROOT::RDataFrame d2("data", noPUFile.c_str());
  auto d2_def = d2.Define("abs_geneta", "fabs(geneta)")
    .Define("en", def1)
    .Define("noise", def2)
    .Define("en_layer", def3);
  d2_def.Foreach(get_weights, {"genen", "abs_geneta", "en", "noise", "en_layer"});
  
  std::string fwname = "root_files/fileweights_"+std::to_string(mask)+std::string(samples)+".root";
  TFile *fileweights = new TFile(fwname.c_str(),"RECREATE");
  fileweights->cd();
  for(int_ ireg=1; ireg<=nreg; ++ireg) {
    //normalize (TProfile does not work with RDataFrame)
    std::string str_s = "en"+std::to_string(ireg)+"_layers_signal";
    for(int_ il=1; il<=nlayers; ++il) {
      int_ bin = hist[str_s]->FindBin(il);
      hist[str_s]->SetBinContent(bin, hist[str_s]->GetBinContent(bin) / histcount(ireg-1,0,il-1));
    }
    for(uint_ iw=1; iw<=bckgcuts.size(); ++iw) {
      std::string str_b = "en"+std::to_string(ireg)+"_layers_bckg"+std::to_string(iw);
      for(int_ il=1; il<=nlayers; ++il) {
	int_ bin = hist[str_b]->FindBin(il);
	hist[str_b]->SetBinContent(bin, hist[str_b]->GetBinContent(bin) / histcount(ireg-1,iw,il-1));
      }
    }
    //write to file
    TH1F* h_s = static_cast<TH1F*>(hist[str_s]->Clone(("en"+std::to_string(ireg)+"_layer_sign").c_str()));
    h_s->Write();
    for(uint_ iw=1; iw<=bckgcuts.size(); ++iw) {
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

  std::string fawname = "root_files/file_after_weights_" + std::to_string(mask) + std::string(samples)+".root";
  TFile *file = new TFile(fawname.c_str(), "RECREATE");
  file->cd();
  TTree *tree = new TTree("data", "Tree after weights");
  float_ geneta;
  vec1d<float_> deltaE(nreg), deltaE_corr(nreg);
  vec1d<float_> shower_number(nreg);
  tree->Branch("geneta", &geneta);
  tree->Branch("deltaE", &deltaE);
  tree->Branch("deltaE_corr", &deltaE_corr);
  tree->Branch("shower_number", &shower_number);

  SoftwareCorrection showercorr("root_files/fileweights_"+std::to_string(mask)+samples+".root");
  vector3d<float_> weights = showercorr.weights;
  vec1d<int_> boundaries(nreg);
  vec1d<float_> lshift(nreg);
  std::string corr_mode;
  if(samples=="inner") {
    boundaries = {5, 5, 5};
    corr_mode = "left";
    lshift = {.65, .59, .48};
  }
  else if(samples=="outer") {
    boundaries = {23, 23, 23};
    corr_mode = "right";
    lshift = {1., 1., 1.};
  }
  vec1d<float_> lowfact = showercorr.low_stats_factor(boundaries, corr_mode);

  auto apply_weights = [&](float_ gen, float_ geta, 
			   vec1d<float_> en, vec1d<float_> noi,
			   vec2d<float_> en_layer) {
    float_ f1=1., f2=0.;
    int_ sid=-1;
    for(int ireg=0; ireg<nreg; ++ireg) {
      f1 = 1.;
      f2 = 0.;
      float_ encalib = 0.;
      typename vec1d<float_>::const_iterator it;
      for(it = etareg.cbegin(); it!=(etareg.cend()-1); ++it) {
	std::string idstr;
	//in case it lies outside the limits of the calibration
	//the event is calibrated with the full calibration region
	if (geta<etareg[0] || geta>etareg[n-1]) {
	  idstr = "sr" + std::to_string(ireg+1) + "from" + 
	    etastr(std::to_string(etareg[0])) + "to" + etastr(std::to_string(etareg[n-1]));
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

      vec1d<float_> ROI_en(nlayers, 0.); 
      for(int_ il=0; il<nlayers; ++il) {
	float_ v = f1*en_layer[ireg][il] - f2;
	if(encalib==0)
	  ROI_en[il] = 0.;
	else
	  ROI_en[il] = v/encalib;
      }

      vec1d<float_> lshift;
      if(samples=="inner")
	lshift = {.65,.59,.48};
      else if(samples=="outer")
	lshift = {1., 1., 1.};
      assert(bckgcuts.size()==lshift.size());
      sid = diff_ed(ROI_en, frac_en[ireg], bckgcuts, 0.05);

      bool_ weight_limit;
      float_ encorr = 0.;
      for (int il=1; il<=nlayers; ++il) {
	if(samples=="inner")
	  weight_limit = il > boundaries[ireg];
	else if(samples=="outer")
	  weight_limit = il < boundaries[ireg];
	float_ v = f1*en_layer.at(ireg).at(il-1) - f2;
	assert(f2==0.);
	assert(f1!=1.);  
	if(sid==0) {
	  encorr += v;
	}
	else if(sid>0 && sid<=static_cast<int_>(bckgcuts.size())) {
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
	if(samples=="inner") encorr *= 1/0.09;
	else if(samples=="outer") encorr *= 1/0.08;
      }
      deltaE_corr.at(ireg) = encorr/gen - 1.;
      shower_number.at(ireg) = sid;
    }
    geneta = geta;
    tree->Fill();
  };

  //ROOT::EnableImplicitMT(ncores);
  ROOT::RDataFrame dfinal("data", noPUFile.c_str());
  dfinal.Define("abs_geneta", "fabs(geneta)")
    .Define("en", def1)
    .Define("noise", def2)
    .Define("en_layer", def3)
    .Foreach(apply_weights, {"genen", "abs_geneta", "en", "noise", "en_layer"});

  file->cd();
  tree->Write();
  file->Close();
  delete file;

  return 0;
}
