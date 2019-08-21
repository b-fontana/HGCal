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
      else if(row[0]=="etareg_"+std::string(samples)+"_fineeta")
	{
	  if( row.size()!=4 ) {row.bad_row();}
	  vec1d<double_> etareg_d = linspace(std::stof(row[1]), std::stof(row[2]), 
					     std::stoi(row[3]));
	  std::copy(etareg_d.begin(), etareg_d.end(), std::back_inserter(etareg));
	}
      else if(row[0]=="nreg") 
	{
	  if( row.size()!=2 ) {row.bad_row();}
	  nreg = std::stoi(row[1]);
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
      else if(row[0]=="bckgcuts_"+std::string(samples)) 
	{
	  if( row.size()!=4 ) {row.bad_row();}
	  bckgcuts = {std:stof(row[1]),std:stof(row[2]),std:stof(row[3])};
	}
    }

  uint_ n = etareg.size();

  //Calibration
  Calibration calibration(mingenen, etareg, nreg, label, samples, 
			  mask, noPUFile, outpath);
  calibration.nopu_calibration(6, false);
  vec1d<mapstr<TF1*>> calib = calibration.calib;

  std::string final = "root_files/final_" + std::to_string(mask) + std::string(samples)+"_fineeta.root";
  TFile *file = new TFile(final.c_str(), "RECREATE");
  file->cd();
  TTree *tree = new TTree("data", "data");
  float_ geneta, genphi;
  vec1d<float_> deltaE_corr(nreg);
  tree->Branch("geneta", &geneta);
  tree->Branch("genphi", &genphi);
  tree->Branch("deltaE_corr", &deltaE_corr);

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

  //produces the weights that will be later used to correct the energy distributions
  auto calibrate_all = [&](float_ gen, float_ geta, float_ gphi, vec1d<float_> en, vec1d<float_> noi) {
    float f1=1., f2=0.;
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
      deltaE_corr.at(ireg-1) = encalib/gen - 1;
      assert(f1!=1.);
    }
    geneta = geta;
    genphi = gphi;
    tree->Fill();
  };

  ROOT::EnableImplicitMT(ncores);
  ROOT::RDataFrame d2("data", noPUFile.c_str());
  auto d2_def = d2.Define("abs_geneta", "fabs(geneta)")
    .Define("en", def1)
    .Define("noi", def2);
  d2_def.Foreach(calibrate_all, {"genen", "abs_geneta", "genphi", "en", "noi"});

  file->cd();
  file->Write();
  delete file;

  return 0;
}
