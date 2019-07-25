#include <iostream>
#include <vector>
#include <iterator>
#include "TFile.h"
#include "TTree.h"
#include "ROOT/RDataFrame.hxx"

#include "interface/utils.h"

void calibrate(float gen, float geta, float en1, float en2, float en3) {
  std::vector<double> en = {en1, en2, en3};
  for(int ireg=1; ireg<4; ++ireg) {
    float f1=1., f2=0.;
    std::vector<float> etareg = {2.7, 2.94};
    std::vector<float> etareg_shift = shift(etareg);
    unsigned int n = etareg.size();

    typename std::vector<float>::iterator it;
    for(it = etareg.begin(); it != etareg.end(); ++it) {
      if(geta<etareg[0] || geta>etareg[n-1]) {
	//std::string idstr = 'sr{}_from{}to{}'.format(ireg, ES(etaregions[0]), 
	//				 ES(etaregions[-1]))

	int i = it-etareg.begin();
	std::cout << *it << "-" << etareg_shift.at(i) << std::endl;
      }
    }
  }
    /*
    etaregshift = np.roll(etareg, shift=-1)[:-1]
    for ieta1,ieta2 in zip(etareg[:-1], etaregshift):
      //in case it lies outside the limits of the calibration
      //the event is calibrated with the full calibration region
      if(geta < etareg[0] or geta > etareg[-1]):
	idstr = 'sr{}_from{}to{}'.format(ireg, ES(etaregions[0]), 
					 ES(etaregions[-1]))
	elif (geta < ieta1 or geta >= ieta2): 
	continue
	else:
	  idstr = 'sr{}_from{}to{}'.format(ireg, ES(ieta1), ES(ieta2))
	    if 'L0' in calib:
	    f1 /= calib['L0'][idstr].Eval(geta)+1.0
	    if 'L1' in calib:
	    f1 /= calib['L1'][idstr].Eval(f1*en[ireg])+1.0    
	    calib_en = f1*en[ireg] - f2 
    */
}

int main() {
  TFile *f = new TFile("summaries/summary_mask3_inner.root");
  TTree *tr = (TTree*)f->Get("data");
  ROOT::RDataFrame d("data", "summaries/summary_mask3_inner.root");
  d.Define("abs_geneta", "abs(geneta)").Foreach(calibrate, {"genen","abs_geneta","en_sr1_ROI","en_sr2_ROI","en_sr3_ROI"});
  return 0;
}
