#ifndef SOFTWARE_CORRECTION_H
#define SOFTWARE_CORRECTION_H

#include <vector>
#include <string>

#include "TFile.h"
#include "TH1.h"
#include "TGraph.h"

#include "UserCode/AnalysisCode/interface/utils.h"
#include "UserCode/AnalysisCode/interface/types.h"

int_ diff_ed(const vec1d<float_>&, const vec1d<float_>&, const vec1d<float_>&, const float_&);

class SoftwareCorrector {
 private:
  int_ nlayers = 28;
  int_ nreg = 3;
  vec1d<float_> bckgcuts = {0.5, 0.6, 0.7};
  vec1d<float_> discrvals = VecOps(bckgcuts).average_array(0.9);
  vec2d<std::string> wn;
  vec1d<std::string> hn_sig;
  vec2d<std::string> hn_bckg;
  TFile* fw;

 public:
  vector3d<float_> weights{static_cast<size_t>(nreg),
      static_cast<size_t>(bckgcuts.size()),
      static_cast<size_t>(nlayers)};

  SoftwareCorrector()=default;
  SoftwareCorrector(std::string);
  ~SoftwareCorrector();
  vec1d<TGraph*> build_weights_graphs(int_);
  vec1d<float_> low_stats_factor(vec1d<uint_>, std::string);
};

#endif //SOFTWARE_CORRECTION_H
