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

int_ main(int_ argc, char_ **argv) 
{
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
  vec1d<std::string> varnames = {"mingenen", "etareg_"+std::string(argv[4]), "etareg_central", "enreg_"+std::string(argv[4]),
				 "nreg", "input", "input_raw", "output"};
  CalibratorInputParameters p("params_pions.csv", varnames, std::string(argv[2]), std::string(argv[4]));
  Calibrator calibrator(p);
  calibrator.create_pion_calibration_values(6, false, true);
  //vec1d<mapstr<TF1*>> calib = calibrator.calibration_values;

  return 0;
}
