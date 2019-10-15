#ifndef PHOTON_ANALYSIS_H
#define PHOTON_ANALYSIS_H

#include <iostream>
#include <map>
#include <vector>
#include <string>
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"

#include "UserCode/AnalysisCode/interface/calibration.h"
#include "UserCode/AnalysisCode/interface/utils.h"
#include "UserCode/AnalysisCode/interface/types.h"
#include "UserCode/AnalysisCode/interface/software_correction.h"
#include "UserCode/AnalysisCode/interface/parser.h"

class PhotonsPartialWafersAnalysis
{
 public:
  PhotonsPartialWafersAnalysis(std::string mask, std::string samples, std::string method);
  const CalibratorInputParameters& get_calibrator_parameters() const;  
  const vec1d<float_>& get_software_correction_shift() const;
  void do_photon_calibration(const vec1d<int_>&, const bool_&, const bool_&, opt<std::string> file_opt = std::nullopt, 
			     opt<vec1d<CalibrationType>> calib_types_opt = std::nullopt, 
			     opt<vec2d<float_>> calib_vars_opt = std::nullopt);

  //functions to pass to RDataFrame methods
  std::function<vec1d<float_>(const vec1d<float_>&, const vec1d<float_>&, const vec1d<float_>&)> define_calibrated_energy();
  std::function<vec2d<float_>(const vec2d<float_>&, const vec1d<float_>&, const vec1d<float_>&)> define_calibrated_energy_per_layer();
  std::function<vec1d<float_>(const float_&, const float_&, const vec1d<float_>&)> define_f1();
  std::function<vec1d<float_>(const float_&, const vec1d<float_>&)> define_f2();
  std::function<bool_(const float_&, const vec1d<float_>&)> complete_showers_filter();
  std::function<vec2d<float_>(const vec2d<float_>&, const vec1d<float_>&, 
			      const vec1d<float_>&, const vec1d<float_>&)> complete_showers_energy_profile();
  std::function<vec1d<int_>(const vec1d<float_>&, const vec2d<float_>&, const vec2d<float_>&)> define_shower_index(const vec2d<float_>&); 
  std::function<vec2d<float_>(const vec1d<float_>&, const vec2d<float_>&)> define_shower_energy_profile(); 
  std::function<vec1d<float_>(const vec1d<float_>&, const vec2d<float_>&, const vec1d<int_>&)> calculate_weight_corrected_energy(SoftwareCorrector&);
  std::function<vec1d<float_>(const float_&, const vec1d<float_>&)> calculate_response();
  std::function<vec1d<float_>(const vec1d<float_>&, const vec1d<int_>&)> apply_low_stats_factor(SoftwareCorrector&);

  std::string define_en = 
    "std::vector<float> en = {en_sr1_ROI, en_sr2_ROI, en_sr3_ROI};"
    "return en;";
  std::string define_noise = 
    "std::vector<float> noise = {noise_sr1_ROI, noise_sr2_ROI, noise_sr3_ROI};"
    "return noise;";
  std::string define_en_layer = 
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

 private:
  vec2d<float_> calib_vars; 
  vec1d<CalibrationType> calib_types;
  vec1d<mapstr<TSpline3*>> calibration_values;
  std::unique_ptr<CalibratorInputParameters> params;
  std::unique_ptr<Calibrator> calibrator;

  //software correction
  vec1d<uint_> boundaries;
  vec1d<float_> lshift;
  std::string corr_mode;

};

#endif //PHOTON_ANALYSIS_H
