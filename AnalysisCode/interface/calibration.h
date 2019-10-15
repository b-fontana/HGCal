#ifndef CALIBRATION_H
#define CALIBRATION_H

#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <algorithm>
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
#include "TSpline.h"

#include "UserCode/AnalysisCode/interface/utils.h"
#include "UserCode/AnalysisCode/interface/types.h"
#include "UserCode/AnalysisCode/interface/parser.h"

enum class CalibrationType {GenEta, RecoEn, GenPhi, PU, NTypes};

class CalibratorInputParameters 
{
 public:
  CalibratorInputParameters() {};
  CalibratorInputParameters(const std::string&, const vec1d<std::string>&, const std::string&, const std::string&, const std::string&, const std::string&);
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
  vec1d<float_> phireg;
  vec1d<float_> etareg_fineeta;
  vec1d<float_> enreg_fineeta;
  vec1d<float_> phireg_fineeta;
  vec1d<float_> bckgcuts;
  std::string label = "";
  std::string outpath = "";
  std::string method = "";
  std::string particle = "";
};

class Calibrator
{
 private:
  enum DataRow {Gen, Geta, En, Noi, Gphi, NElements};
  CalibratorInputParameters p;
  vec1d<float_> etareg_shift;

  vec1d<float_> reconstruct_raw_shower_energy(const uint_& n_particles_per_event, const uint_& ireg,
					      const ROOT::VecOps::RVec<float_>& hits_en,
					      const ROOT::VecOps::RVec<float_>& hits_sr,
					      const ROOT::VecOps::RVec<float_>& hits_roi_idx);
  void draw_spectrum(TH2D*, TGraphAsymmErrors*, TSpline3*, const std::string&, const std::string&, TF1*, const bool_&);
  vec3d<float> get_values_for_calibration(const std::string&, const std::string&, const uint_&, const uint_&, 
					  const CalibrationType& calib_type, opt<vec1d<float_>> calib_var_opt = std::nullopt);
  vec3d<float> get_values_for_compensation(const std::string&, const uint_&, const vec1d< tup3<float_> >, 
					   const float_&, const uint_& n_pions_per_event=2);
  vec1d< tup3<float_> > do_pions_linear_regression(vec4d<float_>, const uint_&);
  vec1d< tup3<float_> > do_pions_linear_regression_python(vec4d<float_>, const uint_&);
  void do_pion_compensation(const uint_&, const vec1d<tup3<float_>>&, const bool_&, const bool_&);
 
 public: 
  //vec1d<mapstr<TF1*>> calibration_values; 
  vec1d<mapstr<TSpline3*>> calibration_values; 

  Calibrator(const CalibratorInputParameters&);
  ~Calibrator();
  void create_pion_calibration_values(const vec1d<int_>&, const bool_&, const bool_&,
				      opt<vec1d<CalibrationType>> calib_types_opt=std::nullopt, 
				      opt<vec2d<float_>> calib_vars_opt=std::nullopt);
  void create_photon_calibration_values(const vec1d<int_>&, const bool_&, const bool_&, 
					opt<std::string> file_opt=std::nullopt, 
					opt<vec1d<CalibrationType>> calib_types_opt=std::nullopt, 
					opt<vec2d<float_>> calib_vars_opt=std::nullopt);
};

#endif //CALIBRATION_H
