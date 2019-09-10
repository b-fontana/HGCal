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
  InputParameters p;
  p.mask = std::stoi(argv[4]);
  p.samples = argv[2];

  std::ifstream infile("params_pions.csv");
  for(CSVIterator it(infile); it != CSVIterator(); ++it)
    {
      CSVRow row = *it;
      if(row[0]=="mingenen") 
	{
	  if( row.size()!=2 ) {row.bad_row();}
	  p.mingenen = std::stof(row[1]);
	}
      else if(row[0]=="etareg_"+std::string(p.samples))
	{
	  if( row.size()!=4 ) {row.bad_row();}
	  if( row[3]<2 ) {row.bad_row();}
	  vec1d<double_> etareg_d = linspace(std::stof(row[1]), std::stof(row[2]), 
					     std::stoi(row[3]));
	  std::copy(etareg_d.begin(), etareg_d.end(), std::back_inserter(p.etareg));
	}
      else if(row[0]=="etareg_central")
	{
	  if( row.size()!=4 ) {row.bad_row();}
	  if( row[3]<2 ) {row.bad_row();}
	  vec1d<double_> etareg_d = linspace(std::stof(row[1]), std::stof(row[2]), 
					     std::stoi(row[3]));
	  std::copy(etareg_d.begin(), etareg_d.end(), std::back_inserter(p.etareg_central));
	}
      else if(row[0]=="enreg_"+std::string(p.samples))
	{
	  if( row.size()!=4 ) {row.bad_row();}
	  if( row[3]<2 ) {row.bad_row();}
	  vec1d<double_> enreg_d = linspace(std::stof(row[1]), std::stof(row[2]), 
					    std::stoi(row[3]));
	  std::copy(enreg_d.begin(), enreg_d.end(), std::back_inserter(p.enreg));
	}
      else if(row[0]=="nreg") 
	{
	  if( row.size()!=2 ) {row.bad_row();}
	  p.nreg = std::stoi(row[1]);
	}
      else if(row[0]=="input") 
	{
	  if( row.size()!=2 ) {row.bad_row();}
	  p.noPUFile = row[1]+"mask"+std::string(argv[4])+"_"+std::string(argv[2])+"_Pions.root";
	}
      else if(row[0]=="input_raw") 
	{
	  if( row.size()!=2 ) {row.bad_row();}
	  p.noPUFile_raw = row[1]+"/Pions/mask"+std::string(argv[4])+"_central/hadd_mask"+std::string(argv[4])+"_central_nopu.root";
	}
      else if(row[0]=="output") 
	{
	  if( row.size()!=2 ) {row.bad_row();}
	  p.outpath = row[1]+std::string(argv[2])+"/mask"+std::string(argv[4])+"_Pions";
	}
    }

  //Calibration  
  Calibration calibration(p);
  calibration.pion_calibration(6, false, true);
  //vec1d<mapstr<TF1*>> calib = calibration.calib;

  return 0;
}
