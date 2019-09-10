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
#include "TBufferJSON.h"

#include "UserCode/AnalysisCode/interface/utils.h"
#include "UserCode/AnalysisCode/interface/types.h"

class Calibration 
{
 private:
  const InputParameters p;
  vec1d<float_> etareg_shift;
  

  TF1* calibrate_spectrum(TH2D*, const std::string&, const std::string&, 
			  const std::string&, const bool_&);
  vec3d<float> values_for_calibration(const std::string&, const std::string&, 
				      const uint_&, const uint_&, std::optional<vec1d<float_>> etas = std::nullopt);
  vec3d<float> values_for_compensation(const vec1d<std::string>&, const uint_&);
  vec1d< tup3<float_> > pions_linear_regression(vec4d<float_>, const uint_&);
  vec1d< tup3<float_> > pions_linear_regression_python(vec4d<float_>, const uint_&);
  void compensation(const int&, const vec1d<tup3<float_>>&, const bool&, const bool&);
  //void fill_layer_energies(const float&, const float&, const float&, const float&);
  void stop() {std::exit(0);}
 
 public: 
  vec1d<mapstr<TF1*>> calib; 

  Calibration(const InputParameters&);

  void pion_calibration(const int_&, const bool_&, const bool_&);
  void photon_calibration(const int&, const bool_&, const bool_&);
};

#endif //CALIBRATION_H
