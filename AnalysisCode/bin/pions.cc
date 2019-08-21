#include <iostream>
#include <vector>
#include <iterator>
#include "TLinearFitter.h"
#include "TVectorD.h"
#include "ROOT/RDataFrame.hxx"

#include "UserCode/AnalysisCode/interface/calibration.h"
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
  //uint_ mask = std::stoi(argv[4]);
  uint_ ncores = std::thread::hardware_concurrency();
  uint_ nsubdetectors = 3;
  float_ mingenen;
  vec1d<float_> etareg;
  //int_ nreg;
  std::string label;
  std::string noPUFile;
  std::string outpath;

  std::ifstream infile("params_pions.csv");
  for(CSVIterator it(infile); it != CSVIterator(); ++it)
    {
      CSVRow row = *it;
      if(row[0]=="mingenen") 
	{
	  if( row.size()!=2 ) {row.bad_row();}
	  mingenen = std::stof(row[1]);
	}
      else if(row[0]=="etareg_"+std::string(samples))
	{
	  if( row.size()!=4 ) {row.bad_row();}
	  if( row[3]<2 ) {row.bad_row();}
	  vec1d<double_> etareg_d = linspace(std::stof(row[1]), std::stof(row[2]), 
					     std::stoi(row[3]));
	  std::copy(etareg_d.begin(), etareg_d.end(), std::back_inserter(etareg));
	}
      /*
      else if(row[0]=="nreg") 
	{
	  if( row.size()!=2 ) {row.bad_row();}
	  nreg = std::stoi(row[1]);
	}
      */
      else if(row[0]=="input") 
	{
	  if( row.size()!=2 ) {row.bad_row();}
	  noPUFile = row[1]+"mask"+std::string(argv[4])+"_"+std::string(argv[2])+"_Pions.root";
	}
      else if(row[0]=="output") 
	{
	  if( row.size()!=2 ) {row.bad_row();}
	  outpath = row[1]+std::string(argv[2])+"/mask"+std::string(argv[4])+"_Pions";
	}
    }
  
  uint_ n = etareg.size();

  //Calibration
  /*
  Calibration calibration(mingenen, etareg, nreg, label, samples, 
			  mask, noPUFile, outpath);
  calibration.nopu_calibration(6, false);
  vec1d<mapstr<TF1*>> calib = calibration.calib;
  */
  int_ ireg=3;
  vec3d<double_> x;
  for(uint_ j=0; j<nsubdetectors; ++j) {
    vec2d<double_> v;
    x.push_back(v);
    for(uint_ i=0; i<n-1; ++i)
      {
	vec1d<double_> v;
	x[j].push_back(v);
      }
  }

  auto fill_calib_vector = [&](int_ idet, float_ geta, float_ en)
    {
      bool_ region_check = false;
      typename std::vector<float_>::const_iterator it;
      for(it = etareg.cbegin(); etareg.cend()-it>1; ++it) {
	if(geta <= *it || geta > *(it+1))
	  continue; 
	assert(!region_check);
	region_check = true;
	int_ idx = it - etareg.cbegin();
	x.at(idet).at(idx).push_back(en);
      }
    };
  auto cut = [](float_ var, float_ cut) {return var > cut;};
  
  std::string mingenen_str = "static_cast<float>(" + std::to_string(mingenen) + ")";
  std::vector<std::string> subdet = {"CEE", "HEF", "HEB"};
  ROOT::EnableImplicitMT(ncores);
  for(uint_ idet=0; idet<nsubdetectors; ++idet)
    {
      std::string idet_str = "static_cast<int>(" + std::to_string(idet) + ")";
      ROOT::RDataFrame d1(subdet[idet].c_str(), noPUFile.c_str());
      d1.Define("abs_geneta", "abs(geneta)")
	.Define("mingenen", mingenen_str)
	.Define("idet", idet_str)
	.Filter(cut, {"genen", "mingenen"})
	.Foreach(fill_calib_vector, {"idet", "abs_geneta", ("en_sr"+std::to_string(ireg)+"_ROI").c_str()});	  
    }

  //linear regression
  for(uint_ ieta=0; ieta<n-1; ++ieta)
    {
      TLinearFitter *lf = new TLinearFitter(nsubdetectors);
      lf->SetFormula("2*x[0] ++ 3*x[1] ++ 4*x[2]");

      assert(x[0][ieta].size() == x[1][ieta].size());
      assert(x[0][ieta].size() == x[2][ieta].size());
      int_ vsize = x[0][ieta].size()+x[1][ieta].size()+x[2][ieta].size();
      double_* xdata = new double_[vsize];
      double_* ydata = new double_[x[0][ieta].size()];
      //double_* edata = new double_[x[0][ieta].size()];
      for(uint_ idata=0; idata<x[0][ieta].size(); ++idata)
	{
	  xdata[0 + idata*nsubdetectors] = x[0][ieta].at(idata);
	  xdata[1 + idata*nsubdetectors] = x[1][ieta].at(idata);
	  xdata[2 + idata*nsubdetectors] = x[2][ieta].at(idata);
	  //edata[idata] = 0.01;
	  ydata[idata] = xdata[0+idata*nsubdetectors] + xdata[1+idata*nsubdetectors] + xdata[2+idata*nsubdetectors];
	}
      lf->AssignData(x[0][ieta].size(), nsubdetectors, xdata, ydata);//, edata);
      lf->Eval();
      TVectorD params;
      TVectorD errors;
      vec1d<double_> significances;
      lf->GetParameters(params);
      for(uint_ i=0; i<nsubdetectors; ++i)
	significances.push_back(lf->GetParSignificance(i));
      lf->GetErrors(errors);

      for(uint_ i=0; i<n-1; ++i)
	{
	  std::cout << "Eta region " << i+1 << std::endl;
	  for(uint_ j=0; j<nsubdetectors; ++j)
	    std::cout << "Subdetector " << j << ": " << params(j) << " +- " << errors(j) << "( significance : " << significances[j] << ")" << std::endl;
	}
      double_ chisquare = lf->GetChisquare();
      std::cout << "chisquare = " << chisquare << std::endl;
      std::cout << std::endl;

      //other formula//
      /*
      lf->SetFormula("pol3");
      lf->Eval();
      lf->GetParameters(params);
      for(uint_ i=0; i<significances.size(); ++i)
	significances[i] = lf->GetParSignificance(i);
      lf->GetErrors(errors);
      for(uint_ i=0; i<n-1; ++i)
	{
	  std::cout << "Eta region " << i+1 << std::endl;
	  for(uint_ j=0; j<nsubdetectors; ++j)
	    std::cout << "Subdetector " << j << ": " << params(j) << " +- " << errors(j) << "( significance : " << significances[j] << ")" << std::endl;
	}
      chisquare = lf->GetChisquare();
      std::cout << "chisquare = " << chisquare << std::endl;
      */

      delete lf;
      delete xdata;
      delete ydata;
      //delete edata;
    }

  return 0;
}
