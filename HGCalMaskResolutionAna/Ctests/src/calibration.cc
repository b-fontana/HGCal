#include "../interface/calibration.h"

Calibration::Calibration(float mingenen_, std::vector<float> etareg_, int nreg_, 
			 std::string label_, std::string samples_, unsigned int mask_, 
			 std::string noPUFile_, std::string outpath_): 
  mingenen(mingenen_), etareg(etareg_), nreg(nreg_), label(label_), 
  samples(samples_), mask(mask_), noPUFile(noPUFile_), outpath(outpath_) {
  etareg_shift = shift(etareg);
  for(int i=0; i<nreg; ++i) {
    cmap tmpmap;
    calib.push_back(tmpmap);
  }
}

void Calibration::pu_calibration(const int& nq, const bool& plot) 
{
}

void Calibration::nopu_calibration(const int& nq, const bool& plot) 
{
  for(int ireg=1; ireg<=nreg; ++nreg) {
    vec3d x = Calibration::energies_for_calibration("data", ireg);
    typename std::vector<float>::iterator it;
    for(it=etareg.begin(); etareg.end()-it>1; ++it) {
      std::cout << "ETAREGIONS ITERATOR: " << *it << " " << *(it+1) << std::endl;
      int idx = std::distance(etareg.begin(), it);
      std::cout << "DISTANCE: " << idx << std::endl;

      std::string str1 = std::to_string(*it).replace(1,1,"p").erase(5,10);
      std::string str2 = std::to_string(*(it+1)).replace(1,1,"p").erase(5,10);
      std::string idstr = "sr" + std::to_string(ireg) + "from" + str1 + "to" + str2;
      std::cout << "STRING: " << idstr << std::endl;

      unsigned int ss = x[idx].size();
      std::cout << "SIZE: " << ss << std::endl;
      double xq0[ss], xq1[ss];
      double probs[nq+1];
      for(int j=0; j<nq+1; ++j) {
	probs[j] = static_cast<double>(j)/nq;
      }
      for(int j=0; j<ss; ++j) {
	xq0[j] = static_cast<double>(x[idx][j][0]);
	xq1[j] = static_cast<double>(x[idx][j][1]);
      }
      double quantiles0[nq+1], quantiles1[nq+1]; 
      TMath::Quantiles(ss, nq+1, xq0, quantiles0, probs, kFALSE);
      TMath::Quantiles(ss, nq+1, xq1, quantiles1, probs, kFALSE);

      //relative calibration versus eta
      std::string hn = idstr + "_" + this->samples + "_mask" + std::to_string(this->mask);
      TH2D* htmp = new TH2D(("resVSeta_"+hn).c_str(),  
			    ";|#eta|;#DeltaE/E", nq, quantiles1, 100, -1, 1);
      float genen, geneta, recen, deltaE;
      for(int i=0; i<ss; ++i) {
	genen = x[idx][i][0];
	geneta = x[idx][i][1];
	recen = x[idx][i][2];
	deltaE = (recen/genen) - 1.;
	htmp->Fill(geneta, deltaE);
      }
      std::string htmp_str = "SR" + std::to_string(ireg);
      this->calib[0][idstr] = Calibration::calibrate_spectrum(htmp, htmp_str, 
							      label+"(PU=0)",
							      "pol2", plot);
      htmp->Delete();
      delete htmp;

      //relative calibration versus energy
      htmp = new TH2D(("resVSen_"+hn).c_str(), ";Reconstructed energy [GeV];#DeltaE/E", 
		      nq, quantiles0, 100, -1, 1);
      for(int i=0; i<ss; ++i) {
	genen = x[idx][i][0];
	geneta = x[idx][i][1];
	recen = x[idx][i][2];
	recen /= (this->calib[0][idstr].Eval(geneta)+1.0);
	deltaE = (recen/genen) - 1.;
	htmp->Fill(recen, deltaE);
      }
      this->calib[1][idstr] = Calibration::calibrate_spectrum(htmp, htmp_str,
							      label+" (PU=0)",
							      "pol1", plot);
      htmp->Delete();
      delete htmp;
    }
  }
}

TF1* Calibration::calibrate_spectrum(TH2D* h, const std::string& title, 
				     const std::string& proc, const std::string& func, 
				     const bool& plot)
{
  TGraphAsymmErrors* prof = build_median_profile(h);
  std::cout << func << std::endl;
  prof->Fit(func.c_str());
  std::string clone_str = h->GetName() + title + "_calib";
  TF1* calibGr = static_cast<TF1*>( prof->GetListOfFunctions()->At(0)->Clone(clone_str.c_str()) );
  std::exit(0);
  if(plot) {
    TCanvas* c = new TCanvas("c", "c", 500, 500);
    c->SetTopMargin(0.05);
    c->SetBottomMargin(0.1);
    c->SetLeftMargin(0.12);
    c->SetRightMargin(0.1);
    h->Draw("colz");
    prof->Draw("e2p");

    TLatex* tex = new TLatex();
    tex->SetTextFont(42);
    tex->SetTextSize(0.04);
    tex->SetNDC();
    tex->DrawLatex(0.12, 0.96, "#bf{CMS} #it{simulation preliminary}");
    tex->DrawLatex(0.15, 0.88, title.c_str());
    tex->SetTextAlign(31);
    tex->DrawLatex(0.97, 0.96, proc.c_str());
    c->SaveAs( (this->outpath + h->GetName() + "_" + title + ".png").c_str() );
    prof->Delete();
    c->Delete();
    delete c;
    delete tex;
  }
  delete prof;
  return calibGr;
}

//correct SR detail
vec3d Calibration::energies_for_calibration(const std::string& tname, 
					    const int& ireg) 
{
  vec3d x;
  //fill energies vector with 2d vectors, one per etaregion considered in the calibration
  for(int i=0; i<this->etareg.size()-1; ++i) {
    vec2d etareg_row;
    x.push_back(etareg_row);
  }

  std::vector<float> newetareg(this->etareg.begin(), this->etareg.end());
  auto fill_layer_energies = [&x, &newetareg](const float& gen, const float& geta, 
					     const float& en, const float& noi) {
    bool region_check = false;
    typename std::vector<float>::iterator it;
    for(it = newetareg.begin(); newetareg.end()-it>1; ++it) {
      if(geta < *it || geta > *(it+1))
	continue; 
      if(region_check) {
	throw "This event was already assigned to an eta region.";
	std::exit(0);
	//Calibration::stop();
      }
      region_check = true;
      float avgnoise = noi * 1;//self.sr_area[ireg-1]/self.sr_area[2]

      std::vector<float> row = {gen, geta, en, noi}; 
      int idx = std::distance(newetareg.begin(), it);
      x.at(idx).push_back(row);
    }
  };
  auto cut = [](float var, float cut) {return var > cut;};

  std::string mingenen_str = "static_cast<float>(" + std::to_string(mingenen) + ")";
  //ROOT::EnableImplicitMT();
  ROOT::RDataFrame d(tname, noPUFile);
  d.Define("abs_geneta", "abs(geneta)")
    .Define("mingenen", mingenen_str)
    .Filter(cut, {"genen", "mingenen"})
    .Foreach(fill_layer_energies, {"genen","abs_geneta","en_sr1_ROI","noise_sr3_ROI"});
  return x;
}
