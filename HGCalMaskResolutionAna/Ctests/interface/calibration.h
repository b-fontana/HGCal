#ifndef CALIBRATION_H
#define CALIBRATION_H

#include <iostream>
#include <map>
#include <vector>
#include "TGraph.h"
#include "TH1.h"
#include "TF1.h"
#include "TTree.h"

class Calibration 
{
 private:
  typedef std::vector< std::vector<double> > vec2d;
  typedef std::map<std::string, TGraph> cmap; //calibration storage map

  TF1* calibrate_spectrum(const TH1&, const std::string&, const std::string&, 
			  const std::string&, const bool&);
  vec2d energies_for_calibration(TTree*, const int&);  

 public:
  cmap calib;

  Calibration();
  void pu_calibration(const int&, const bool&);
  void nopu_calibration(const int&, const bool&);
};

#endif //CALIBRATION_H
