#ifndef SOFTWARE_CORRECTION_H
#define SOFTWARE_CORRECTION_H

#include <vector>
#include <string>

#include "TFile.h"
#include "TH1.h"
#include "TGraph.h"
#include "../interface/utils.h"
#include "../interface/types.h"

class SoftwareCorrection {
 private:
  int_ nlayers = 28;
  int_ nreg = 3;
  vec1d<float_> bckgcuts = {0.5, 0.6, 0.7};
  vec1d<float_> discrvals = VecOps(bckgcuts).average_array(0.9);
  vec3d<float_> weights;
  vec2d<std::string> wn;
  vec1d<std::string> hn_sig;
  vec2d<std::string> hn_bckg;
  TFile* fw;

 public:
  SoftwareCorrection(std::string);
  vec1d<TGraph*> build_weights_graphs(int);
  vec1d<float_> low_stats_factor(vec1d<float_>, std::string);

};

#endif //SOFTWARE_CORRECTION_H
