#include <iostream>
#include <vector>
#include <iterator>
#include "TLinearFitter.h"
#include "TF1.h"
#include "TRandom.h"
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
  uint_ mask = std::stoi(argv[4]);
  float_ mingenen;
  vec1d<float_> etareg;
  uint_ nreg;
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
      else if(row[0]=="nreg") 
	{
	  if( row.size()!=2 ) {row.bad_row();}
	  nreg = std::stoi(row[1]);
	}
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

  //Calibration
  Calibration calibration(mingenen, etareg, nreg, label, samples, 
			  mask, noPUFile, outpath);
  calibration.pion_calibration(6, false, true);
  //vec1d<mapstr<TF1*>> calib = calibration.calib;

  return 0;
}
