#ifndef CALIBRATION_H
#define CALIBRATION_H

#include <iostream>
#include <map>
#include <vector>
#include <string>
#include "ROOT/RDataFrame.hxx"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TH2D.h"
#include "TF1.h"
#include "TRandom.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TFile.h"
#include "TLinearFitter.h"
#include "TVectorD.h"

#include "UserCode/AnalysisCode/interface/utils.h"
#include "UserCode/AnalysisCode/interface/types.h"

class Calibration 
{
 private:
  const float_ mingenen;
  vec1d<float_> etareg;
  vec1d<float_> etareg_shift;
  const uint_ nreg;
  const std::string label;
  const std::string samples;
  const uint_ mask;
  const std::string noPUFile;
  const std::string outpath;
  const uint_ nsubdets;

  TF1* calibrate_spectrum(TH2D*, const std::string&, const std::string&, 
			  const std::string&, const bool_&);
  vec3d<float> energies_for_calibration(const std::string&, const int&);
  std::vector< std::tuple<float_, float_, float_> > pions_linear_regression(vec4d<float_>, const uint_&);
  //void fill_layer_energies(const float&, const float&, const float&, const float&);
  void stop() {std::exit(0);}
 
 public: 
  vec1d<mapstr<TF1*>> calib; 

  Calibration(const float_, const std::vector<float_>, const uint_,
	      const std::string, const std::string, const uint_,
	      const std::string, const std::string);

  void pion_calibration(const int_&, const bool_&, const bool_&);
  void photon_calibration(const int&, const bool_&, const bool_&);
};

#endif //CALIBRATION_H
