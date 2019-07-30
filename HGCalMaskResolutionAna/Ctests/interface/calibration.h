#ifndef CALIBRATION_H
#define CALIBRATION_H

#include <iostream>
#include <map>
#include <vector>
#include <string>
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TH2D.h"
#include "TF1.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TFile.h"
#include "ROOT/RDataFrame.hxx"

#include "../interface/utils.h"
#include "../interface/types.h"

class Calibration 
{
 private:
  const float mingenen;
  vec1d<float> etareg;
  vec1d<float> etareg_shift;
  const std::string label;
  const std::string samples;
  const unsigned int mask;
  const std::string noPUFile;
  const std::string outpath;
  const int nreg;

  TF1* calibrate_spectrum(TH2D*, const std::string&, const std::string&, 
			  const std::string&, const bool&);
  vec3d<float> energies_for_calibration(const std::string&, const int&);  
  //void fill_layer_energies(const float&, const float&, const float&, const float&);
  void stop() {std::exit(0);}
 
 public: 
  vec1d<mapfunc> calib; 

  Calibration(const float, const std::vector<float>, const int,
	      const std::string, const std::string, const unsigned int,
	      const std::string, const std::string);

  void pu_calibration(const int&, const bool&);
  void nopu_calibration(const int&, const bool&);
};

#endif //CALIBRATION_H
