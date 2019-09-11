#ifndef CALIBRATION_H
#define CALIBRATION_H

#include <iostream>
#include <map>
#include <vector>
#include <string>
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
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
#include "UserCode/AnalysisCode/interface/parser.h"

class CalibratorInputParameters 
{
 public:
  CalibratorInputParameters() {};
  CalibratorInputParameters(const std::string&, const vec1d<std::string>&, const std::string&, const std::string&);
  //CalibratorInputParameters(CalibratorInputParameters&);
  void define_input_parameters(const std::string&, const vec1d<std::string>&, const std::string&, const std::string&);

  float_ mingenen = 20.;
  uint_ nreg = 3;
  uint_ nlayers = 0;
  std::string samples = "";
  uint_ mask = 3;
  uint_ nsubdets = 3;
  std::string noPUFile = "";
  std::string noPUFile_raw = "";
  vec1d<float_> etareg;
  vec1d<float_> etareg_central;
  vec1d<float_> enreg;
  vec1d<float_> bckgcuts;
  std::string label = "";
  std::string outpath = "";
};

class Calibrator
{
 private:
  CalibratorInputParameters p;
  vec1d<float_> etareg_shift;
  

  TF1* calibrate_spectrum(TH2D*, const std::string&, const std::string&, 
			  const std::string&, const bool_&);
  vec3d<float> get_values_for_calibration(const std::string&, const std::string&, const uint_&, const uint_&, 
				      std::optional<vec1d<float_>> etas = std::nullopt);
  vec3d<float> get_values_for_compensation(const std::string&, const uint_&, const float_&, const uint_& n_pions_per_event=2);
  vec1d< tup3<float_> > do_pions_linear_regression(vec4d<float_>, const uint_&);
  vec1d< tup3<float_> > do_pions_linear_regression_python(vec4d<float_>, const uint_&);
  void do_pion_compensation(const uint_&, const vec1d<tup3<float_>>&, const bool_&, const bool_&);
 
 public: 
  vec1d<mapstr<TF1*>> calibration_values; 

  Calibrator(const CalibratorInputParameters&);
  void create_pion_calibration_values(const int_&, const bool_&, const bool_&);
  void create_photon_calibration_values(const int&, const bool_&, const bool_&);
};

#endif //CALIBRATION_H
