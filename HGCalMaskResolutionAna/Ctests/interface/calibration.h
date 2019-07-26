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

typedef std::vector< std::vector<float> > vec2d;
typedef std::vector< std::vector< std::vector<float> > > vec3d;
typedef std::map<std::string, TGraph> cmap; //calibration storage map

class Calibration 
{
 private:
  float mingenen;
  std::vector<float> etareg;
  std::vector<float> etareg_shift;
  std::string label;
  std::string samples;
  unsigned int mask;
  std::string noPUFile;
  std::string outpath;
  int nreg;

  TF1* calibrate_spectrum(TH2D*, const std::string&, const std::string&, 
			  const std::string&, const bool&);
  vec3d energies_for_calibration(const std::string&, const int&);  
  //void fill_layer_energies(const float&, const float&, const float&, const float&);
  void stop() {std::exit(0);}
 
 public: 
  std::vector<cmap> calib; 

  Calibration(float, std::vector<float>, int, 
	      std::string, std::string, unsigned int, 
	      std::string, std::string);

  void pu_calibration(const int&, const bool&);
  void nopu_calibration(const int&, const bool&);
};

#endif //CALIBRATION_H
