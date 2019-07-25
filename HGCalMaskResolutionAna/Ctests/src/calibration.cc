#include "../interface/calibration.h"
#include "../interface/utils.h"

Calibration::Calibration() {
  float mingenen = 20.;
  std::vector<float> etareg = {2.7, 2.94};
  std::vector<float> etareg_shift = shift(etareg);
  //unsigned int n = etareg.size();
  std::string label = "";
  std::string samples = "inner";
  unsigned int mask = 3;
  std::string outpath = "."; //summaries/summary_mask"${imask}"_"${isample}".root
}

void Calibration::pu_calibration(const int& nq, const bool& plot) 
{
}

void Calibration::nopu_calibration(const int& nq, const bool& plot) 
{
}

TF1* Calibration::calibrate_spectrum(const TH1& h, const std::string& title, 
				     const std::string& proc, const std::string& func, 
				     const bool& plot)
{
}

Calibration::vec2d Calibration::energies_for_calibration(TTree* data, const int& ireg) 
{
}
