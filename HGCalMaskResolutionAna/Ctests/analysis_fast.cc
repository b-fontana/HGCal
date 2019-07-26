#include <iostream>
#include <vector>
#include <iterator>
#include "TFile.h"
#include "TTree.h"
#include "ROOT/RDataFrame.hxx"

#include "interface/utils.h"
#include "interface/calibration.h"

int main() {
  float mingenen = 20.;
  std::vector<float> etareg = {2.7, 2.94};
  int nreg = 3;
  std::string label = "";
  std::string samples = "inner";
  unsigned int mask = 3;
  std::string noPUFile = "../summaries/summary_mask3_inner.root";
  std::string outpath = "../pics_inner/mask3";
  
  Calibration calibration(mingenen, etareg, mask, label, samples, 
			  mask, noPUFile, outpath);
  calibration.nopu_calibration(6, true);
  std::vector<cmap> calib = calibration.calib;

  auto calibrate = [&calib](float gen, float geta, float en1, 
			   float en2, float en3, float noi) {
    std::vector<float> en = {en1, en2, en3};
    for(int ireg=1; ireg<4; ++ireg) {
      float f1=1., f2=0.;
      std::vector<float> etareg = {2.7, 2.94};
      std::vector<float> etareg_shift = shift(etareg);
      unsigned int n = etareg.size();

      typename std::vector<float>::iterator it;
      for(it = etareg.begin(); it != etareg.end(); ++it) {
	int idx = it-etareg.begin();
	std::string idstr;
	//in case it lies outside the limits of the calibration
	//the event is calibrated with the full calibration region
	if (geta<etareg[0] || geta>etareg[n-1]) {
	  idstr = "sr" + std::to_string(ireg) + "_from" + 
	    std::to_string(etareg[0]) + "to" + std::to_string(etareg[n-1]);
	}
	else if (geta<*it || geta>*(it+1))
	  continue;
	else {
	  idstr = "sr" + std::to_string(ireg) + "_from" + 
	    std::to_string(*it) + "to" + std::to_string(*(it+1));
	  if (calib[0].size() > 0) {
	    f1 /= calib[0][idstr].Eval(geta)+1.0;
	    if (calib[1].size() > 0) {
	      f1 /= calib[1][idstr].Eval(f1*en[ireg])+1.0;
	      if (calib[2].size() > 0)
		f2 = calib[2][idstr].Eval(noi);
	    }
	  }
	  en[ireg-1] = f1*en[ireg-1] - f2;
	}
      }
    }
  };

  //ROOT::EnableImplicitMT();
  ROOT::RDataFrame d("data", noPUFile.c_str());
  d.Define("abs_geneta", "abs(geneta)")
    .Foreach(calibrate, {"genen","abs_geneta","en_sr1_ROI","en_sr2_ROI","en_sr3_ROI","avgnoise_sr3"});

  return 0;
}
