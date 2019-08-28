#include "UserCode/AnalysisCode/interface/calibration.h"

Calibration::Calibration(const float_ mingenen_, std::vector<float_> etareg_, 
			 const uint_ nreg_, const std::string label_, 
			 const std::string samples_, const uint_ mask_, 
			 const std::string noPUFile_, const std::string outpath_): 
  mingenen(mingenen_), etareg(etareg_), nreg(nreg_), label(label_), 
  samples(samples_), mask(mask_), noPUFile(noPUFile_), outpath(outpath_), nsubdets(3) {
  etareg_shift = VecOps(etareg).shift();
  for(uint_ i=0; i<nreg; ++i) {
    mapstr<TF1*> tmpmap;
    calib.push_back(tmpmap);
  }
}

/*Performs a 3d linear regression to obtain the calibration factor of
  each subdetector. For the result to be accurate the provided sample
  must contain data in the three subdetectors (central eta region).

  The function being fitted is of the form:
  y = a0*x + a1*y + a2*z
  where each variable (x, y and z) corresponds to the total shower energy deposited
  in their corresponding subdetector (CEE, HEF and HEB).

  //structure of x: 
  1D: subdetectors, size=3
  2D: etaregions being considered
  3D: values per events in the following order: {gen energy, gen eta, reco energy, noise}
*/
vec1d< tup3<float_> > Calibration::pions_linear_regression(vec4d<float_> x, const uint_& reg)
{
  vec1d< tup3<float_> > factors;
  TRandom randNum;
  TLinearFitter *lf;
  for(uint_ ieta=0; ieta<this->etareg.size()-1; ++ieta)
    {
      assert(x[0][ieta].size() == x[1][ieta].size());

      lf = new TLinearFitter();
      lf->SetFormula("x++y++z");
      /*
      for(uint_ idet=0; idet<this->nsubdets; ++idet)
	lf->SetParameter(idet, ("a"+std::to_string(idet)).c_str(), .33, .01, 0., 1);
      */
      int_ vsize = x[0][ieta].size()*this->nsubdets;
      double_* xdata = new double_[vsize];
      double_* ydata = new double_[x[0][ieta].size()];
      double_* edata = new double_[x[0][ieta].size()];
      for(uint_ idata=0; idata<x[0][ieta].size(); ++idata) 
	{
	  //TLinearFitter only works with doubles
	  xdata[0 + idata*this->nsubdets] = static_cast<double>(x[0][ieta][idata][2]);
	  xdata[1 + idata*this->nsubdets] = static_cast<double>(x[1][ieta][idata][2]);
	  xdata[2 + idata*this->nsubdets] = static_cast<double>(x[2][ieta][idata][2]);
	  edata[idata] = x[0][ieta][idata][0]/50.;
	  //gen information is repeated; I could have used x[1] or x[2]
	  assert(x[0][ieta][idata][0] == x[1][ieta][idata][0]);
	  assert(x[0][ieta][idata][0] == x[2][ieta][idata][0]);
	  ydata[idata] = x[0][ieta][idata][0] + randNum.Gaus()*edata[idata];
	}
      lf->AssignData(x[0][ieta].size(), this->nsubdets, xdata, ydata, edata);
      lf->Eval();
      TVectorD params;
      TVectorD errors;
      vec1d<double_> significances;
      lf->GetParameters(params);
      lf->GetErrors(errors);
      for(uint_ i=0; i<this->nsubdets; ++i) 
	significances.push_back(lf->GetParSignificance(i));
      factors.push_back( tup3<float_>(params(0), params(1), params(2)) );
      for(uint_ j=0; j<this->nsubdets; ++j) 
	{
	  std::cout << "Subdetector " << j << ": " << params(j) << " +- " << errors(j) << std::endl;
	  for(uint_ k=0; k<x[j][ieta].size(); ++k)
	    assert(x[j][ieta][k][2] == xdata[j + k*this->nsubdets]);
	}
      double_ chisquare = lf->GetChisquare();
      std::cout << "chisquare = " << chisquare << std::endl;
      
      delete[] xdata;
      delete[] ydata;
      delete[] edata;
      delete lf;
    }
  return factors;
}

vec1d< tup3<float_> > Calibration::pions_linear_regression_python(vec4d<float_> x, const uint_& reg)
{
  for(uint_ idet=0; idet<this->nsubdets; ++idet) 
    {
      for(uint_ i=0; i<this->etareg.size()-1; ++i) 
	{
	  std::string fname = "linear_regression_det"+std::to_string(idet)+"_eta"+std::to_string(i)+".json";
	  std::ofstream o(fname);
	  auto json = TBufferJSON::ToJSON(&x[idet][i]);
	  o << json << std::endl;
	}
    }
  
  std::exit(0);
  vec1d< tup3<float_> > factors;
  return factors;
}

void Calibration::pion_calibration(const int& nq, const bool& pu, const bool& plot) 
{
  if(pu==true)
    {
      std::cout << "PU calibration is not yet implemented." << std::endl;
      std::exit(0);
    }

  std::vector<std::string> subdet = {"CEE", "HEF", "HEB"};
  vec4d<float_> x;
  for(uint_ ireg=1; ireg<=nreg; ++ireg) 
    {
      x.clear();
      std::string hstr = "SR" + std::to_string(ireg);
      for(uint_ idet=0; idet<this->nsubdets; ++idet)
	{
	  vec3d<float_> vtmp; 
	  x.push_back(vtmp);
	  x[idet] = Calibration::energies_for_calibration(subdet[idet].c_str(), ireg);
	}

      //obtain the factors from a 3d linear regression
      std::cout << "----------------------------" << std::endl;
      vec1d< tup3<float_> > flr = Calibration::pions_linear_regression(x, ireg);
      std::cout << "----------------------------" << std::endl;

      typename std::vector<float_>::const_iterator it;
      for(it=etareg.cbegin(); it!=(etareg.cend()-1); ++it) {
	int_ idx = std::distance(etareg.cbegin(), it);
	float_ etatmp1 = *it;      
	float_ etatmp2 = *(it+1);
	std::string str1 = std::to_string(etatmp1).replace(1,1,"p").erase(5,10);
	std::string str2 = std::to_string(etatmp2).replace(1,1,"p").erase(5,10);
	std::string idstr = "sr" + std::to_string(ireg) + "from" + str1 + "to" + str2;

	uint_ ss = x[0][idx].size();
	double_ xq0[ss], xq1[ss];
	double_ probs[nq+1];
	for(int_ j=0; j<nq+1; ++j)
	  probs[j] = static_cast<double_>(j)/nq;
	float_ deltaE;
	vec1d<float_> genen, geneta, recen;
	for(uint_ j=0; j<ss; ++j) {
	  //std::cout << j << ", " << x[0][idx].at(j)[2] << ", " << x[1][idx][j][2] << ", " << x[2][idx][j][2] << std::endl;
	  /*The calibration should only be done for detectors that have particles that produced completely contained showers
	    Calibrating individual subdetectors would then be a bad idea (see Sec.6.2.11 Calorimetry, Wigmans, 2nd ed.).
	  */
	  //sum the reconstructed energy in the three subdetectors
	  genen.push_back(x[0][idx][j][0]);
	  geneta.push_back(x[0][idx][j][1]);
	  recen.push_back( std::get<0>(flr[idx])*x[0][idx][j][2] + std::get<1>(flr[idx])*x[1][idx][j][2] 
			   + std::get<2>(flr[idx])*x[2][idx][j][2] );
	  xq0[j] = static_cast<double_>(genen.at(j));
	  xq1[j] = static_cast<double_>(geneta.at(j));
	}
	double quantiles0[nq+1], quantiles1[nq+1]; 
	TMath::Quantiles(ss, nq+1, xq0, quantiles0, probs, kFALSE);
	TMath::Quantiles(ss, nq+1, xq1, quantiles1, probs, kFALSE);

	//relative calibration versus eta
	std::string hn = idstr + "_" + this->samples + "_mask" + std::to_string(this->mask);
	TH2D* htmp1 = new TH2D(("resVSeta_"+hn).c_str(), ";|#eta|;#DeltaE/E", 
			       nq, quantiles1, 250, -1, 3);
	for(uint_ i=0; i<ss; ++i) {
	  deltaE = (recen[i]/genen[i]) - 1.;
	  htmp1->Fill(geneta[i], deltaE);
	}
	this->calib[0][idstr] = Calibration::calibrate_spectrum(htmp1, hstr, label+"(PU=0)",
								"pol2", plot);
	delete htmp1;

	//relative calibration versus energy
	TH2D *htmp2 = new TH2D(("resVSen_"+hn).c_str(), ";Reco energy [GeV];#DeltaE/E", 
			       nq, quantiles0, 250, -1, 3);
	for(uint_ i=0; i<ss; ++i) {
	  recen[i] /= (this->calib[0][idstr]->Eval(geneta[i])+1.0);
	  deltaE = (recen[i]/genen[i]) - 1.;
	  htmp2->Fill(recen[i], deltaE);
	}
	this->calib[1][idstr] = Calibration::calibrate_spectrum(htmp2, hstr, label+" (PU=0)",
								"pol1", plot);
	delete htmp2;
      }
    }
}

void Calibration::photon_calibration(const int& nq, const bool& pu, const bool& plot) 
{
  if(pu==true)
    {
      std::cout << "PU calibration is not yet implemented." << std::endl;
      std::exit(0);
    }

  for(uint_ ireg=1; ireg<=nreg; ++ireg) {      
    std::string hstr = "SR" + std::to_string(ireg);
    vec3d<float_> x = Calibration::energies_for_calibration("data", ireg);
    typename std::vector<float_>::const_iterator it;

    //loop over eta calibration regions
    for(it=etareg.cbegin(); it!=(etareg.cend()-1); ++it) {
      int_ idx = std::distance(etareg.cbegin(), it);

      float_ etatmp1 = *it;      
      float_ etatmp2 = *(it+1);
      std::string str1 = std::to_string(etatmp1).replace(1,1,"p").erase(5,10);
      std::string str2 = std::to_string(etatmp2).replace(1,1,"p").erase(5,10);
      std::string idstr = "sr" + std::to_string(ireg) + "from" + str1 + "to" + str2;

      uint_ ss = x[idx].size();
      double_ xq0[ss], xq1[ss];
      double_ probs[nq+1];
      for(int_ j=0; j<nq+1; ++j)
	probs[j] = static_cast<double_>(j)/nq;
 
      for(uint_ j=0; j<ss; ++j) {
	xq0[j] = static_cast<double_>(x[idx][j][2]); //reco energy quantiles
	xq1[j] = static_cast<double_>(x[idx][j][1]); //gen eta quantiles
      }
      double quantiles0[nq+1], quantiles1[nq+1]; 
      TMath::Quantiles(ss, nq+1, xq0, quantiles0, probs, kFALSE);
      TMath::Quantiles(ss, nq+1, xq1, quantiles1, probs, kFALSE);

      //relative calibration versus eta
      std::string hn = idstr + "_" + this->samples + "_mask" + std::to_string(this->mask);
      TH2D* htmp1 = new TH2D(("resVSeta_"+hn).c_str(), ";|#eta|;#DeltaE/E", 
			     nq, quantiles1, 100, -1, 1);
      float_ genen, geneta, recen, deltaE;
      for(uint_ i=0; i<ss; ++i) {
	genen = x[idx][i][0];
	geneta = x[idx][i][1];
	recen = x[idx][i][2];
	deltaE = (recen/genen) - 1.;
	htmp1->Fill(geneta, deltaE);
      }
      this->calib[0][idstr] = Calibration::calibrate_spectrum(htmp1, hstr, label+"(PU=0)",
							      "pol2", plot);
      delete htmp1;

      //relative calibration versus energy
      TH2D *htmp2 = new TH2D(("resVSen_"+hn).c_str(), ";Reco energy [GeV];#DeltaE/E", 
			     nq, quantiles0, 100, -1, 1);
      for(uint_ i=0; i<ss; ++i) {
	genen = x[idx][i][0];
	geneta = x[idx][i][1];
	recen = x[idx][i][2];
	recen /= (this->calib[0][idstr]->Eval(geneta)+1.0);
	deltaE = (recen/genen) - 1.;
	htmp2->Fill(recen, deltaE);
      }
      this->calib[1][idstr] = Calibration::calibrate_spectrum(htmp2, hstr, label+" (PU=0)",
							      "pol1", plot);
      delete htmp2;
    }
  }
}

TF1* Calibration::calibrate_spectrum(TH2D* h, const std::string& title, 
				     const std::string& proc, const std::string& func, 
				     const bool& plot)
{
  TGraphAsymmErrors* prof = build_median_profile(h);
  prof->Fit(func.c_str());
  std::string clone_str = h->GetName() + title + "_calib";
  TF1* calibGr = static_cast<TF1*>(prof->GetListOfFunctions()->At(0)->Clone(clone_str.c_str()) );
 
  if(plot) {
    TCanvas c("c", "c", 500, 500);
    c.SetTopMargin(0.05);
    c.SetBottomMargin(0.1);
    c.SetLeftMargin(0.12);
    c.SetRightMargin(0.1);
    h->Draw("colz");
    prof->Draw("e2p");

    TLatex tex;
    tex.SetTextFont(42);
    tex.SetTextSize(0.04);
    tex.SetNDC();
    tex.DrawLatex(0.12, 0.96, "#bf{CMS} #it{simulation preliminary}");
    tex.DrawLatex(0.15, 0.88, title.c_str());
    tex.SetTextAlign(31);
    tex.DrawLatex(0.97, 0.96, proc.c_str());
    c.SaveAs( (this->outpath + h->GetName() + "_" + title + ".png").c_str() );
  }

  prof->Delete();
  return calibGr;
}

vec3d<float_> Calibration::energies_for_calibration(const std::string& tname, 
						    const int& ireg) 
{
  vec3d<float_> x;
  //fill energies vector with 2d vectors, one per etaregion considered in the calibration
  for(uint_ i=0; i<this->etareg.size()-1; ++i) {
    vec2d<float_> etareg_row;
    x.push_back(etareg_row);
  }

  std::vector<float_> newetareg(this->etareg.cbegin(), this->etareg.cend());
  auto fill_layer_energies = [&x, &newetareg](const float& gen, const float& geta, 
					      const float& en, const float& noi) {
    bool_ region_check = false;
    typename std::vector<float_>::const_iterator it;
    for(it = newetareg.cbegin(); newetareg.cend()-it>1; ++it) {
      if(geta <= *it || geta > *(it+1))
	continue; 
      assert(!region_check);
      region_check = true;
      //float_ avgnoise = noi * self.sr_area[ireg-1]/self.sr_area[2]

      std::vector<float_> row = {gen, geta, en, noi}; 
      int_ idx = it - newetareg.cbegin();
      x.at(idx).push_back(row);
    }
  };
  auto cutmin = [](float_ var, float_ cut) {return var > cut;};

  std::string mingenen_str = "static_cast<float>(" + std::to_string(mingenen) + ")";
  uint_ ncores = std::thread::hardware_concurrency();
  ROOT::EnableImplicitMT(ncores);
  ROOT::RDataFrame d(tname, noPUFile);
  d.Define("abs_geneta", "fabs(geneta)")
    .Define("mingenen", mingenen_str)
    .Filter(cutmin, {"genen", "mingenen"})
    .Foreach(fill_layer_energies, {"genen", "abs_geneta",
	  ("en_sr"+std::to_string(ireg)+"_ROI").c_str(), "noise_sr3_ROI"});
  return x;
}
