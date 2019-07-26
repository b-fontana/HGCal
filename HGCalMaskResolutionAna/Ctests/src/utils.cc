#include "../interface/utils.h"

TGraphAsymmErrors* build_median_profile(TH2D* h) {
  TGraphAsymmErrors* medianGr = new TGraphAsymmErrors();
  std::string graph_name = std::string(h->GetName())+"_medianprof";
  medianGr->SetName( graph_name.c_str() );
  medianGr->SetLineWidth(2);
  medianGr->SetMarkerStyle(20);

  //median and 1 sigma quantiles
  double xq[3] = {0.16, 0.5, 0.84};
  //storage for the x values that correspond to the above quantiles
  double yq[3] = {0.0, 0.0, 0.0};
  
  for(int xbin=1; xbin<h->GetNbinsX()+1; ++xbin) {
    TH1D* tmp = h->ProjectionY("tmp", xbin, xbin);
    tmp->GetQuantiles(3, yq, xq);
    float xcen = h->GetXaxis()->GetBinCenter(xbin);

    int npts = medianGr->GetN();
    medianGr->SetPoint(npts, xcen, yq[1]);
    medianGr->SetPointError(npts, 0, 0, yq[1]-yq[0], yq[2]-yq[1]);
    tmp->Delete();
  }
  return medianGr;
}
