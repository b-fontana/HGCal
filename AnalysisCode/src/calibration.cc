#include "UserCode/AnalysisCode/interface/calibration.h"

Calibrator::Calibrator(const CalibratorInputParameters& p_): p(p_)
{
  for(uint_ i=0; i<to_underlying(CalibrationType::NTypes); ++i)
    {
      mapstr<TSpline3*> tmpmap;
      calibration_values.push_back(tmpmap);
   }
}

Calibrator::~Calibrator()
{
  for(uint_ i=0; i<this->calibration_values.size(); ++i)
    {
      for(auto const& x : this->calibration_values[i])
	delete x.second;
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
vec1d< tup3<float_> > Calibrator::do_pions_linear_regression(vec4d<float_> x, const uint_& reg)
{
  vec1d< tup3<float_> > factors;
  TRandom randNum;
  TLinearFitter *lf;
  for(uint_ ieta=0; ieta<this->p.etareg_central.size()-1; ++ieta)
    {
      assert(x[0][ieta].size() == x[1][ieta].size());
      lf = new TLinearFitter();
      lf->SetFormula("x++y++z");
      int_ vsize = x[0][ieta].size()*this->p.nsubdets;
      double_* xdata = new double_[vsize];
      double_* ydata = new double_[x[0][ieta].size()];
      double_* edata = new double_[x[0][ieta].size()];
      for(uint_ idata=0; idata<x[0][ieta].size(); ++idata)
	{
	  //TLinearFitter.AssignData only works with doubles
	  xdata[0 + idata*this->p.nsubdets] = static_cast<double>(x[0][ieta][idata][2]);
	  xdata[1 + idata*this->p.nsubdets] = static_cast<double>(x[1][ieta][idata][2]);
	  xdata[2 + idata*this->p.nsubdets] = static_cast<double>(x[2][ieta][idata][2]);
	  edata[idata] = x[0][ieta][idata][0]/50.;
	  //gen information is repeated; I could have used x[1] or x[2]
	  assert(x[0][ieta][idata][0] == x[1][ieta][idata][0]);
	  assert(x[0][ieta][idata][0] == x[2][ieta][idata][0]);
	  ydata[idata] = x[0][ieta][idata][0] + randNum.Gaus()*edata[idata];
	}
      lf->AssignData(x[0][ieta].size(), this->p.nsubdets, xdata, ydata, edata);
      lf->Eval();
      TVectorD params;
      TVectorD errors;
      vec1d<double_> significances;
      lf->GetParameters(params);
      lf->GetErrors(errors);
      for(uint_ i=0; i<this->p.nsubdets; ++i)
	significances.push_back(lf->GetParSignificance(i));
      factors.push_back( tup3<float_>(params(0), params(1), params(2)) );
      for(uint_ j=0; j<this->p.nsubdets; ++j) 
	{
	  std::cout << "Subdetector " << j << ": " << params(j) << " +- " << errors(j) << std::endl;
	  for(uint_ k=0; k<x[j][ieta].size(); ++k)
	    assert(x[j][ieta][k][2] == xdata[j + k*this->p.nsubdets]);
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

/*
Requires more work on importing the result of the regression
*/
vec1d< tup3<float_> > Calibrator::do_pions_linear_regression_python(vec4d<float_> x, const uint_& reg)
{
  for(uint_ idet=0; idet<this->p.nsubdets; ++idet)
    {
      for(uint_ i=0; i<this->p.etareg_central.size()-1; ++i)
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

/*return perhaps a std::function that describes the Fermi function,
  plus store new TF1 in the calib dictionary for later evaluation*/
void Calibrator::do_pion_compensation(const uint_& ireg, const vec1d<tup3<float_>>& flr, const bool_& pu, const bool_& draw_plot)
{
  if(pu==true)
    {
      std::cout << "PU calibration is not yet implemented." << std::endl;
      std::exit(0);
    }
  for(uint_ i=0; i<to_underlying(CalibrationType::NTypes); ++i)
    {
      if(this->calibration_values[i].size() == 0 && i!=to_underlying(CalibrationType::PU))
	{
	  std::cout << "The calibration functions have to be defined before calling do_pion_compensation()!" << std::endl;
	  std::exit(0);
	}
    }

  vec1d<float_> energy_mip_limits = {5.};
  for(auto it = energy_mip_limits.cbegin(); it != energy_mip_limits.cend(); ++it)
    {
      vec3d<float_> compvals;
      compvals = get_values_for_compensation("an_mask/CEE_HEF_HEB", ireg, flr, *it);
      for(auto it2 = this->p.enreg.cbegin(); it2 != this->p.enreg.cend()-1; ++it2)
	{
	  //int_ idx = std::distance(p.enreg.cbegin(), it2);
	  std::string enstr = std::to_string(*it2) + "_" + std::to_string(*(it2+1));
	  TH2D* h = new TH2D(("Cglob_"+enstr).c_str(), ";|C_{glob}|;#DeltaE/E", 20, 0.3, 1.3, 100, -1, 1);
	  delete h;
	}
    }
}

void Calibrator::create_pion_calibration_values(const vec1d<int_>& nq, const std::vector<std::string>& en_calib_str,
						const bool_& pu, const bool_& draw_plot,
						opt<vec1d<CalibrationType>> calib_types_opt,
						opt<vec3d<float_>> calib_vars_opt)
{
  if(pu==true)
    {
      std::cout << "PU calibration is not yet implemented." << std::endl;
      std::exit(0);
    }
  if( (calib_types_opt != std::nullopt and calib_vars_opt == std::nullopt) or
      (calib_types_opt == std::nullopt and calib_vars_opt != std::nullopt) )
    {
      std::cout << "When you specified one of the options in create_photon_calibration_values() you"
	" also have to specify the other one. This is a precaution against non-sense calibrations." << std::endl;
      std::exit(0);
    }

  vec1d<CalibrationType> calib_types;
  if(calib_types_opt == std::nullopt)
    calib_types = {CalibrationType::GenEta, CalibrationType::RecoEn, CalibrationType::GenPhi};
  else
    calib_types = calib_types_opt.value();

  //first dimension: variable type; second dimension: integrated signal region
  vec3d<float_> calib_vars;
  if(calib_vars_opt == std::nullopt)
    {
      for(uint_ j=0; j<calib_types.size(); ++j)
	{
	  vec2d<float_> tmp;
	  calib_vars.emplace_back(tmp);
	}
      for(uint_ j=0; j<p.nreg; ++j)
	{
	  calib_vars[0].emplace_back(this->p.phireg);
	  calib_vars[1].emplace_back(this->p.etareg);
	  calib_vars[2].emplace_back(this->p.enreg);
	}
      assert(calib_types.size() == 3);
    }
  else
    calib_vars = calib_vars_opt.value();
  assert(calib_types.size() == calib_vars.size());
  assert(nq.size() == calib_vars.size());

  umap< CalibrationType, std::tuple< DataRow,std::string, std::string> > type_table = {
    {CalibrationType::GenEta, std::make_tuple(DataRow::Geta, "resVSeta_", ";|#eta|;#DeltaE/E")},
    {CalibrationType::RecoEn, std::make_tuple(DataRow::En, "resVSen_", ";Reco energy [GeV];#DeltaE/E")},
    {CalibrationType::GenPhi, std::make_tuple(DataRow::Gphi, "resVSphi_", ";|#phi|;#DeltaE/E")}};

  //I am assuming the three signal regions have the same limits on all the variables
  umap< CalibrationType, std::pair<float_,float_> > allowed_calibration_window;
  for(auto it_type_window = calib_types.cbegin(); it_type_window != calib_types.cend(); ++it_type_window)
    {
      int_ index_ = std::distance(calib_types.cbegin(), it_type_window);
      allowed_calibration_window[*it_type_window] = std::make_pair(calib_vars[index_][p.nreg-1][0],calib_vars[index_][p.nreg-1].back());
    }
  umap< CalibrationType, std::pair<float_,float_> > allowed_calibration_window_central;
  allowed_calibration_window[CalibrationType::GenEta] = std::make_pair(this->p.etareg_central[0],this->p.etareg_central.back());

  std::string fregression = ( "../HGCalMaskResolutionAna/summaries/summary_mask" + std::to_string(this->p.mask) 
			      + "_central_Pions.root" );
  std::vector<std::string> subdet = {"CEE", "HEF", "HEB"};
  vec4d<float_> x, x_regression;
  for(uint_ ireg=1; ireg<=p.nreg; ++ireg)
    {
      x_regression.clear();
      for(uint_ idet=0; idet<this->p.nsubdets; ++idet)
	{
	  vec3d<float_> vtmp;
	  x_regression.push_back(vtmp);
	  //only one bin; it does not ot matter which variable we choose
	  x_regression[idet] = Calibrator::get_values_for_calibration(fregression, "summary", en_calib_str[ireg-1],
								      allowed_calibration_window_central,
								      CalibrationType::GenEta, this->p.etareg_central);
	}
      //obtain the factors from a 3d linear regression
      vec1d< tup3<float_> > flr = Calibrator::do_pions_linear_regression(x_regression, ireg);
      if(flr.size() != 1)
	{
	  std::cout << "I am using more than one eta bin in the central region. Modifications are required when accessing flr.";
	  std::exit(0);
	}

      for(auto it_type = calib_types.cbegin(); it_type != calib_types.cend(); ++it_type)
	{
	  x.clear();
	  uint_ type_idx = std::distance(calib_types.cbegin(), it_type);
	  for(uint_ idet=0; idet<this->p.nsubdets; ++idet)
	    {
	      vec3d<float_> vtmp;
	      x.push_back(vtmp);
		//TODO:problem in the third argument
	      x[idet] = Calibrator::get_values_for_calibration(this->p.noPUFile, "summary", en_calib_str[idet],
							       allowed_calibration_window,
							       *it_type, calib_vars[type_idx][ireg-1]);
	    }
	  for(auto it = calib_vars[type_idx][ireg-1].cbegin(); it!=calib_vars[type_idx][ireg-1].cend()-1; ++it)
	    {
	      uint_ idx = std::distance(calib_vars[type_idx][ireg-1].cbegin(), it);
	      float_ vartmp1 = *it;
	      float_ vartmp2 = *(it+1);
	      std::string str1 = round_to_string(vartmp1);
	      std::string str2 = round_to_string(vartmp2);
	      std::string idstr = "sr" + std::to_string(ireg) + "from" + str1 + "to" + str2;

	      uint_ ss = x[0][idx].size();
	      double_ xq[ss];
	      double_ probs[nq[type_idx]+1];
	      for(int_ j=0; j<nq[type_idx]+1; ++j)
		probs[j] = static_cast<double_>(j)/nq[type_idx];
	      float_ deltaE;
	      vec1d<float_> genen, geneta, genphi, recen;
	      for(uint_ j=0; j<ss; ++j) {
		/*The calibration should only be done for detectors that have particles that produced completely contained showers
		  Calibrating individual subdetectors would then be a bad idea (see Sec.6.2.11 Calorimetry, Wigmans, 2nd ed.).
		*/
		//sum the reconstructed energy in the three subdetectors
		genen.push_back(x[0][idx][j][DataRow::Gen]);
		geneta.push_back(x[0][idx][j][DataRow::Geta]);
		genphi.push_back(x[0][idx][j][DataRow::Gphi]);
		recen.push_back( std::get<0>(flr[0])*x[0][idx][j][DataRow::En] +
				 std::get<1>(flr[0])*x[1][idx][j][DataRow::En] +
				 std::get<2>(flr[0])*x[2][idx][j][DataRow::En] );
		using ttup = std::tuple<float_, vec1d<float_>>;
		umap<CalibrationType, ttup> var_table = {{CalibrationType::GenEta, std::make_tuple(geneta[j], this->p.etareg)}, 
							 {CalibrationType::RecoEn, std::make_tuple(recen[j], this->p.enreg)},
							 {CalibrationType::GenPhi, std::make_tuple(genphi[j], this->p.phireg)}};
		for(uint_ k=0; k<type_idx; ++k)
		  {
		    CalibrationType type = calib_types[k];

		    //The reco energy has to be corrected not with the bin of the current variable, but with the bin
		    //corresponding to the variable being used for the correction
		    float vartmp1_other, vartmp2_other;
		    vec1d<float_> bins = std::get<1>(var_table[type]);
		    if( std::get<0>(var_table[type]) <= calib_vars[k][ireg-1][0] )
		      {
			vartmp1_other = calib_vars[k][ireg-1][0];
			vartmp2_other = calib_vars[k][ireg-1][1];
		      }
		    else if( std::get<0>(var_table[type]) >= calib_vars[k][ireg-1].back() )
		      {
			vartmp1_other = calib_vars[k][ireg-1].end()[-2];
			vartmp2_other = calib_vars[k][ireg-1].back();
		      }
		    else
		      {
			auto const it_bins = std::lower_bound(bins.begin(), bins.end(), std::get<0>(var_table[type]));
			vartmp1_other = *(it_bins-1);
			vartmp2_other = *it_bins;
		      }
		    std::string str1_other = round_to_string(vartmp1_other);
		    std::string str2_other = round_to_string(vartmp2_other);
		    std::string idstr_other = "sr" + std::to_string(ireg) + "from" + str1_other + "to" + str2_other;

		    //the for loop assures that each additional variable being used for the calibration implies
		    //an additional correction for the reconstructed energy
		    recen.back() /= (this->calibration_values[to_underlying(type)][idstr_other]->Eval( std::get<0>(var_table[type])) + 1.);
		  }

		DataRow drow = std::get<0>(type_table[*it_type]);
		xq[j] = static_cast<double_>( x[0][idx][j][ drow ] );
	      }
	      double quantiles[nq[type_idx]+1];
	      TMath::Quantiles(ss, nq[type_idx]+1, xq, quantiles, probs, kFALSE);

	      std::string hn = idstr + "_" + this->p.samples + "_mask" + std::to_string(this->p.mask);
	      TH2D* htmp = new TH2D( (std::get<1>(type_table[*it_type]) + hn).c_str(), (std::get<2>(type_table[*it_type])).c_str(), 
				     nq[type_idx], quantiles, 250, -1, 3);
	      for(uint_ i=0; i<ss; ++i) {
		umap<CalibrationType, float_> vt = {{CalibrationType::GenEta, geneta[i]},
						    {CalibrationType::RecoEn, recen[i]},
						    {CalibrationType::GenPhi, genphi[i]}};
		deltaE = (recen[i]/genen[i]) - 1.;
		htmp->Fill(vt[*it_type], deltaE);
	      }

	      TGraphAsymmErrors* gae = build_median_profile(htmp);
	      this->calibration_values[to_underlying(*it_type)][idstr] = new TSpline3(("spline_"+idstr).c_str(), gae);
	      if(draw_plot)
		Calibrator::draw_spectrum(htmp, gae, this->calibration_values[to_underlying(*it_type)][idstr], "SR"+std::to_string(ireg), 
					  p.label+"(PU=0)");
	      delete htmp;
	    }
	}

      //Needs further work.
      do_pion_compensation(ireg, flr, false, false);
    }
  std::exit(0);
}

void Calibrator::create_photon_calibration_values(const vec1d<int_>& nq, const std::vector<std::string>& en_calib_str,
						  const bool_& pu, const bool_& draw_plot,
						  opt<std::string> file_opt,
						  opt<vec1d<CalibrationType>> calib_types_opt,
						  opt<vec3d<float_>> calib_vars_opt)
{
  if(pu==true)
    {
      std::cout << "PU calibration is not yet implemented." << std::endl;
      std::exit(0);
    }
  if( (calib_types_opt != std::nullopt and calib_vars_opt == std::nullopt) or
      (calib_types_opt == std::nullopt and calib_vars_opt != std::nullopt) )
    {
      std::cout << "When you specified one of the options in create_photon_calibration_values() you"
	" also have to specify the other one. This is a precaution against non-sense calibrations." << std::endl;
      std::exit(0);
    }
  std::string file = file_opt.value_or(this->p.noPUFile);
  vec1d<CalibrationType> calib_types;
  if(calib_types_opt == std::nullopt)
    calib_types = {CalibrationType::GenPhi, CalibrationType::GenEta, CalibrationType::RecoEn};
  else
    calib_types = calib_types_opt.value();
  vec3d<float_> calib_vars;
  if(calib_vars_opt == std::nullopt)
    {
      for(uint_ j=0; j<calib_types.size(); ++j)
	{
	  vec2d<float_> tmp;
	  calib_vars.push_back(tmp);
	}
      for(uint_ j=0; j<p.nreg; ++j)
	{
	  calib_vars[0].push_back(this->p.phireg);
	  calib_vars[1].push_back(this->p.etareg);
	  calib_vars[2].push_back(this->p.enreg);
	}
      assert(calib_types.size() == 3);
    }
  else
    calib_vars = calib_vars_opt.value();

  assert(calib_types.size() == calib_vars.size());
  assert(nq.size() == calib_vars.size());
  umap<CalibrationType, std::tuple< DataRow, std::string, std::string> > type_table = {
    {CalibrationType::GenEta, std::make_tuple(DataRow::Geta, "resVSeta_", ";|#eta|;#DeltaE/E")},
    {CalibrationType::RecoEn, std::make_tuple(DataRow::En, "resVSen_", ";Reco energy [GeV];#DeltaE/E")},
    {CalibrationType::GenPhi, std::make_tuple(DataRow::Gphi, "resVSphi_", ";|#phi|;#DeltaE/E")}};

  //I am assuming the three signal regions have the same limits on all the variables
  umap< CalibrationType, std::pair<float_,float_> > allowed_calibration_window;
  for(auto it_type_window = calib_types.cbegin(); it_type_window != calib_types.cend(); ++it_type_window)
    {
      int_ index_ = std::distance(calib_types.cbegin(), it_type_window);
      allowed_calibration_window[*it_type_window] = std::make_pair(calib_vars[index_][p.nreg-1][0],calib_vars[index_][p.nreg-1].back());
    }

  for(uint_ ireg=1; ireg<=p.nreg; ++ireg) {
    std::string hstr = "SR" + std::to_string(ireg);
    //loop over calibration variables
    for(auto it_type = calib_types.cbegin(); it_type != calib_types.cend(); ++it_type)
      {
	uint_ type_idx = std::distance(calib_types.cbegin(), it_type);
	vec3d<float_> x = Calibrator::get_values_for_calibration(file, "data", en_calib_str[ireg-1], allowed_calibration_window,
								 *it_type, calib_vars[type_idx][ireg-1]);

	for(auto it = calib_vars[type_idx][ireg-1].cbegin(); it!=calib_vars[type_idx][ireg-1].cend()-1; ++it)
	  {
	    uint_ idx = std::distance(calib_vars[type_idx][ireg-1].cbegin(), it);
	    float_ vartmp1 = *it;
	    float_ vartmp2 = *(it+1);
	    std::string str1 = round_to_string(vartmp1);
	    std::string str2 = round_to_string(vartmp2);
	    std::string idstr = "sr" + std::to_string(ireg) + "from" + str1 + "to" + str2;

	    uint_ ss = x[idx].size();
	    std::cout << ss << ", " << *it << ", " << *(it+1) << std::endl;
	    double_ probs[nq[type_idx]+1];
	    for(int_ j=0; j<nq[type_idx]+1; ++j)
	      probs[j] = static_cast<double_>(j)/nq[type_idx];
	    double_ xq[ss];
	    DataRow drow = std::get<0>(type_table[*it_type]);
	    for(uint_ j=0; j<ss; ++j)
		xq[j] = static_cast<double_>( x[idx][j][ drow ] );
	    double quantiles[nq[type_idx]+1];
	    TMath::Quantiles(ss, nq[type_idx]+1, xq, quantiles, probs, kFALSE);
	    //relative calibration versus eta
	    std::string hn = idstr + "_" + this->p.samples + "_mask" + std::to_string(this->p.mask);
	    TH2D* htmp = new TH2D( (std::get<1>(type_table[*it_type]) + hn).c_str(), std::get<2>(type_table[*it_type]).c_str(), 
				   nq[type_idx], quantiles, 100, -1, 1);
	    float_ genen, recen, deltaE;
	    for(uint_ i=0; i<ss; ++i)
	      {
		genen = x[idx][i][DataRow::Gen];
		recen = x[idx][i][DataRow::En];

		umap<CalibrationType, std::tuple<float_, vec1d<float_>> > var_table;
		for(auto it2_type = calib_types.cbegin(); it2_type != calib_types.cend(); ++it2_type)
		  {
		    int_ index_ = std::distance(calib_types.cbegin(), it2_type);
		    var_table[*it2_type] = std::make_tuple(x[idx][i][std::get<0>(type_table[*it2_type])],
							   calib_vars[index_][ireg-1]);
		  }
		for(uint_ k=0; k<type_idx; ++k)
		  {
		    CalibrationType type = calib_types[k];

		    //The reco energy has to be corrected not with the bin of the current variable, but with the bin
		    //corresponding to the variable being used for the correction
		    float vartmp1_other, vartmp2_other;
		    vec1d<float_> bins = std::get<1>(var_table[type]);
		    if( std::get<0>(var_table[type]) <= calib_vars[k][ireg-1][0] )
		      {
			vartmp1_other = calib_vars[k][ireg-1][0];
			vartmp2_other = calib_vars[k][ireg-1][1];
		      }
		    else if( std::get<0>(var_table[type]) >= calib_vars[k][ireg-1].back() )
		      {
			vartmp1_other = calib_vars[k][ireg-1].end()[-2];
			vartmp2_other = calib_vars[k][ireg-1].back();
		      }
		    else
		      {
			auto const it_bins = std::lower_bound(bins.begin(), bins.end(), std::get<0>(var_table[type]));
			vartmp1_other = *(it_bins-1);
			vartmp2_other = *it_bins;
		      }
		    std::string str1_other = round_to_string(vartmp1_other);
		    std::string str2_other = round_to_string(vartmp2_other);
		    std::string idstr_other = "sr" + std::to_string(ireg) + "from" + str1_other + "to" + str2_other;

		    //the for loop assures that each additional variable being used for the calibration implies
		    //an additional correction for the reconstructed energy
		    recen /= (this->calibration_values[to_underlying(type)][idstr_other]->Eval(std::get<0>(var_table[type]))+1.0);
		  }

		deltaE = (recen/genen) - 1.;
		htmp->Fill(std::get<0>(var_table[*it_type]), deltaE);
	      }

	    TGraphAsymmErrors* gae = build_median_profile(htmp);
	    this->calibration_values[to_underlying(*it_type)][idstr] = new TSpline3(("spline_"+idstr).c_str(), gae);
	    if(draw_plot)
	      Calibrator::draw_spectrum(htmp, gae, this->calibration_values[to_underlying(*it_type)][idstr], "SR" + std::to_string(ireg), 
					p.label+"(PU=0)");
	    delete htmp;
	  }
      }
  }
}

void Calibrator::draw_spectrum(TH2D* h, TGraphAsymmErrors* g, TSpline3* s, const std::string& title, const std::string& proc)
{
  TCanvas c("c", "c", 500, 500);
  c.SetTopMargin(0.05);
  c.SetBottomMargin(0.1);
  c.SetLeftMargin(0.12);
  c.SetRightMargin(0.1);
  h->Draw("colz");
  g->Draw("e2p");
  s->SetLineColor(kRed);
  s->Draw("same");

  TLatex tex;
  tex.SetTextFont(42);
  tex.SetTextSize(0.04);
  tex.SetNDC();
  tex.DrawLatex(0.12, 0.96, "#bf{CMS} #it{simulation preliminary}");
  tex.DrawLatex(0.15, 0.88, title.c_str());
  tex.SetTextAlign(31);
  tex.DrawLatex(0.97, 0.96, proc.c_str());
  c.SaveAs( (this->p.outpath + h->GetName() + "_" + title + ".png").c_str() );
}

bool Calibrator::is_outside_calibration_window(const umap< CalibrationType, std::pair<float_,float_> >& window,
					       const umap< CalibrationType, float_>& var_table) const
{
  for(auto const& [type, ranges]: window)
    {
      if(var_table.at(type) < ranges.first or var_table.at(type) > ranges.second)
	return true;
    }
  return false;
}

vec3d<float_> Calibrator::get_values_for_calibration(const std::string& fname, const std::string& tname,
						     const std::string& en_calib_str,
						     const umap< CalibrationType, std::pair<float_,float_> >& calib_window,
						     const CalibrationType& calib_type, opt<vec1d<float_>> calib_var_opt) const
{
  umap<CalibrationType, vec1d<float_>> type_table = {{CalibrationType::GenEta, this->p.etareg},
						     {CalibrationType::RecoEn, this->p.enreg},
						     {CalibrationType::GenPhi, this->p.phireg}};
  vec1d<float_> calib_var;
  if(calib_var_opt == std::nullopt)
    calib_var = type_table[calib_type];
  else
    calib_var = calib_var_opt.value();
  assert(calib_var.size() != 0);
  vec3d<float_> x;
  //fill energies vector with 2d vectors, one per etaregion considered in the calibration
  for(uint_ i=0; i<calib_var.size()-1; ++i)
    {
      vec2d<float_> var_row;
      x.push_back(var_row);
    }
  
  std::mutex mut;
  auto fill_layer_energies = [&](const float& gen, const float& geta, const float& gphi,
				 const float& en, const float& noi)
    {
      const umap< CalibrationType, float_> var_table = {{CalibrationType::GenEta, geta},
							{CalibrationType::RecoEn, en},
							{CalibrationType::GenPhi, gphi}};
      float_ binning_var = var_table.at(calib_type);
      bool_ region_check = false;
      for(auto it = calib_var.cbegin(); calib_var.cend()-it>1; ++it)
	{
	  if(binning_var <= *it or binning_var > *(it+1))
	    continue;
	  assert(!region_check);

	  region_check = true;
	  if( is_outside_calibration_window(calib_window, var_table) )
	    break;
	  //float_ avgnoise = noi * self.sr_area[ireg-1]/self.sr_area[2]

	  vec1d<float_> row(DataRow::NElements,0.);
	  row[DataRow::Gen] = gen;
	  row[DataRow::Geta] = geta;
	  row[DataRow::En] = en;
	  row[DataRow::Noi] = noi;
	  row[DataRow::Gphi] = gphi;

	  int_ idx = std::distance(calib_var.cbegin(), it);
	  { //mutex lock scope
	    std::lock_guard lock(mut);
	    x.at(idx).push_back(row);
	  }
	}
    };

  ROOT::RDataFrame d(tname, fname);
  d.Foreach(fill_layer_energies, {"genen", "abs_geneta", "compress_phi", en_calib_str.c_str(), "noise_sr3_ROI"});
  return x;
}

/*Iterates over the hits, calibrates them, joins them form three subdetectors
  and returns (Cglob, E_reco) event by event*/
vec3d<float_> Calibrator::get_values_for_compensation(const std::string& tname, const uint_& ireg,
						      const vec1d< tup3<float_> > flr,
						      const float_& energy_mip_limit, const uint_& n_pions_per_event)
{
  /*
  size_t n = p.etareg.size();
  auto convert_str = [](const float_& etatmp)
    {
      return std::to_string(etatmp).replace(1,1,"p").erase(5,10);
    };
  */

  vec3d<float_> values;
  //fill vector with 2d vectors, one per energy region considered in the calibration
  for(uint_ i=0; i<this->p.enreg.size()-1; ++i) {
    vec2d<float_> enreg_row;
    values.push_back(enreg_row);
  }

  std::mutex mut;
  auto calculate_cglob = [&, this](const float& gen, const float& geta,
				   const ROOT::VecOps::RVec<float_>& hits_en, const ROOT::VecOps::RVec<int_>& hits_subdet,
				   const ROOT::VecOps::RVec<int_>& hits_sr, const ROOT::VecOps::RVec<int_>& hits_roi_idx)
    {
      vec1d<float_> raw_recen = this->reconstruct_raw_shower_energy(n_pions_per_event, ireg, hits_en, hits_sr, hits_roi_idx);
      //this->calibrate_shower()
      //this->calculate_cglob()

      /*
      typename std::vector<float_>::const_iterator it;
      for(it = this->p.etareg.cbegin(); it != this->p.etareg.cend()-1; ++it) {
	std::string idstr;
	if (geta < this->p.etareg[0] or geta > this->p.etareg[n-1])
	  {
	    idstr = "sr" + std::to_string(ireg) + "from" + convert_str(p.etareg[0]) + "to" + convert_str(p.etareg[n-1]);
	  }
	else if (geta < *it or geta > *(it+1))
	  continue;
	else
	  {
	    idstr = "sr" + std::to_string(ireg) + "from" + convert_str(*it) + "to" + convert_str(*(it+1));
	  }

	//it should only get here once per event
	float_ f1 = 1., f2 = 0.;
	f1 /= calibration_values[0][idstr]->Eval(geta)+1.0;
	f1 /= calibration_values[1][idstr]->Eval(f1*en[ireg-1])+1.0;
	if (calibration_values[2].size() > 0)
	  f2 = calibration_values[2][idstr]->Eval(noi[2]);

	//all hits are being calibrated, including the noise
	ROOT::VecOps::RVec<float_> hits_en_calib(hits_en.size(), 0.);
	typename ROOT::VecOps::RVec<float_>::const_iterator it_hiten;
	for(it_hiten = hits_en.cbegin(); it_hiten != hits_en.cend(); ++it_hiten)
	  {
	    int_ idx = std::distance(hits_en.cbegin(), it_hiten);
	    for(uint_ idet=0; idet<this->p.nsubdets; ++idet)
	      {
		if(hits_subdet[idx] == idet)
		  {
		    const uint_ idet_const = idet;
		    hits_en_calib[idx] = f1*std::get<idet_const>(flr[0])*hits_en[idx] - f2;
		  }
	      }
	  }
      }
      */

      /*Third part: Cglob calculation*/
      /*
      bool_ region_check = false;
      typename std::vector<float_>::const_iterator it;
      for(it = this->p.enreg.cbegin(); it != this->p.enreg.cend()-1; ++it)
	{
	  if(gen <= *it || gen > *(it+1))
	    continue; 
	  assert(!region_check);
	  region_check = true;

	  vec1d<int_> n_under_mean(n_pions_per_event, 0);
	  vec1d<int_> n_under_limit(n_pions_per_event, 0);
	  vec1d<float_> mean(n_pions_per_event, 0.);
	  vec1d<int_> count(n_pions_per_event, 0);

	  //calculation of the average RecHit energy in all the regions of interest
	  typename ROOT::VecOps::RVec<float_>::const_iterator it_hiten;
	  for(it_hiten = hits_en_calib.cbegin(); it_hiten != hits_en_calib.cend(); ++it_hiten)
	    {
	      uint_ idx1 = std::distance(hits_en_calib.cbegin(), it_hiten);
	      for(uint_ iev=0; iev<n_pions_per_event; ++iev)
		{
		  if(hits_sr.at(idx1) <= static_cast<int_>(ireg) and hits_roi_idx.at(idx1) == static_cast<int_>(iev))
		    {
		      mean.at(iev) += *it_hiten;
		      count.at(iev) += 1;
		    }
		}
	    }
	  for(uint_ iev=0; iev<n_pions_per_event; ++iev)
	    {
	      if(count.at(iev) < 2)
		mean.at(iev) = -1;
	      else
		mean.at(iev) /= static_cast<float_>(count.at(iev));
	    }

	  //calculation of the Cglob factor for all regions of interest
	  typename ROOT::VecOps::RVec<float_>::const_iterator it2_hiten;
	  for(it2_hiten = hits_en_calib.cbegin(); it2_hiten != hits_en_calib.cend(); ++it2_hiten)
	    {
	      uint_ idx2 = std::distance(hits_en_calib.cbegin(), it2_hiten);
	      for(uint_ iev=0; iev<n_pions_per_event; ++iev)
		{
		  if(hits_sr.at(idx2) <= static_cast<int_>(ireg) and hits_roi_idx.at(idx2) == static_cast<int_>(iev))
		    {
		      if(*it2_hiten <= mean[iev])
			n_under_mean.at(iev) += 1;
		      if(*it2_hiten <= energy_mip_limit)
			n_under_limit.at(iev) += 1;
		    }
		}
	    }
	  vec1d<float_> cglob(n_pions_per_event, 0.);
	  int_ idx = std::distance(this->p.enreg.cbegin(), it);
	  for(uint_ iev=0; iev<n_pions_per_event; ++iev)
	    {
	      if(mean.at(iev) == -1) //there was 1 hits or none for this region of interest inside the specified signal region
		continue;
	      cglob.at(iev) = static_cast<float_>(n_under_limit.at(iev)) / n_under_mean.at(iev);
	      vec1d<float_> row = {gen, geta, recen[iev], cglob.at(iev)};
	      { //mutex lock scope
		std::lock_guard lock(mut);
		values.at(idx).push_back(row);
	      }
	    }
	}
      */
    };

  auto cutmin = [](float_ var, float_ cut)
    {
      return var > cut;
    };

  std::string mingenen_str = "static_cast<float>(" + std::to_string(this->p.mingenen) + ")";
  std::string genen_str = "static_cast<float>(ROIs.pt_[0]) * static_cast<float>(TMath::CosH(static_cast<double>(ROIs.eta_[0])))";
  std::string abseta_str = "static_cast<float>(ROIs.eta_[0])"; //always positive
  ROOT::RDataFrame d(tname, this->p.noPUFile_raw);
  d.Define("mingenen", mingenen_str)
    .Define("genen", genen_str)
    .Define("abs_geneta", abseta_str)
    .Filter(cutmin, {"genen", "mingenen"})
    .Foreach(calculate_cglob, {"genen", "abs_geneta", "Hits.en_mip_", "Hits.subdet_", "Hits.signalRegion_", "Hits.assocROIidx_"});
  return values;
}

vec1d<float_> Calibrator::reconstruct_raw_shower_energy(const uint_& n_particles_per_event, const uint_& ireg,
							const ROOT::VecOps::RVec<float_>& hits_en,
							const ROOT::VecOps::RVec<float_>& hits_sr,
							const ROOT::VecOps::RVec<float_>& hits_roi_idx)
{
  vec1d<float_> raw_recen(n_particles_per_event, 0.);
  typename ROOT::VecOps::RVec<float_>::const_iterator it_hiten;
  for(it_hiten = hits_en.cbegin(); it_hiten != hits_en.cend(); ++it_hiten)
    {
      uint_ idx1 = std::distance(hits_en.cbegin(), it_hiten);
      for(uint_ iev=0; iev<n_particles_per_event; ++iev)
	{
	  if(hits_sr.at(idx1) <= static_cast<int_>(ireg) and hits_roi_idx.at(idx1) == static_cast<int_>(iev))
	    raw_recen[iev] += *it_hiten;
	}
    }
  return raw_recen;
}

/*The input parameters are obtained from a CSV file where each row specifies the values of the parameter.
The parameter name is always the first item in each row.
The name of the file has to contain one of the following subtrings: 'photon' or 'pion'.
*/
CalibratorInputParameters::CalibratorInputParameters(const std::string& fname, const vec1d<std::string>& varnames, 
						     const std::string& mask, const std::string& samples, const std::string& method,
						     const std::string& particle):
phireg_fineeta(this->nreg, vec1d<float_>(0)), etareg_fineeta(this->nreg, vec1d<float_>(0)), enreg_fineeta(this->nreg, vec1d<float_>(0))
{
  if(std::stoi(mask) > 6 || std::stoi(mask) < 3)
    {
      std::cout << "The mask has to be 3, 4, 5 or 6." << std::endl;
      std::exit(0);
    }
  if(samples != mDetectorRegion[DetectorRegion::Inner] && samples != mDetectorRegion[DetectorRegion::Outer])
    {
      std::cout << "The samples are either 'inner' or 'outer'." << std::endl;
      std::exit(0);
    }
  this->samples = samples;
  this->mask = std::stoi(mask);
  this->method = method;
  if(fname.find("photon") != std::string::npos)
    this->particle = mParticleType[ParticleType::Photon];
  else if(fname.find("pion") != std::string::npos)
    this->particle = mParticleType[ParticleType::Pion];
  else
    {
      std::cout << "The configuration file must refer to photons or pions!" << std::endl;
      std::exit(0);
    }
  CalibratorInputParameters::define_input_parameters(fname, varnames, mask, samples);
}

void CalibratorInputParameters::define_input_parameters(const std::string& fname, const vec1d<std::string>& varnames,
							const std::string& mask, const std::string& samples)
{
  vec1d<std::string> rows;
  std::ifstream infile(fname);
  for(CSVIterator row_it(infile); row_it != CSVIterator(); ++row_it)
    {
      CSVRow row = *row_it;
      rows.push_back(row[0]);
      if(row[0] == "mingenen")
	{
	  if( row.size()!=2 ) {row.bad_row();}
	  this->mingenen = std::stof(row[1]);
	}
      else if(row[0].find("etareg") != std::string::npos)
	{
	  if( row.size()<3 ) {row.bad_row();}
	  vec1d<double_> etareg_d;
	  if(row[0].find(mMethod[Method::BruteForce]) != std::string::npos)
	    {
	      if(row[0].find("mask"+mask) != std::string::npos)
		{
		  uint_ i = 0;
		  while(i < this->nreg)
		    {
		      if(row[0] == "etareg_"+samples+"_fineeta_sr"+std::to_string(i+1)+"_mask"+mask)
			{
			  etareg_d.clear();
			  for(uint_ j=1; j<row.size(); ++j)
			      etareg_d.push_back(std::stof(row[j]));
			  check_vector_is_sorted(etareg_d);
			  this->etareg_fineeta[i].resize(etareg_d.size());
			  for(uint_ k=0; k<etareg_d.size(); ++k)
			    this->etareg_fineeta[i][k] = static_cast<float_>(etareg_d[k]);
			}
		      ++i;
		    }
		}
	    }
	  else
	    {
	      if(row[0].find(samples) != std::string::npos)
		{
		  if( row.size()!=4 ) {row.bad_row();}
		  etareg_d = linspace(std::stod(row[1]), std::stod(row[2]), std::stoi(row[3]));
		  if(row[0].find("central") != std::string::npos)
		    std::copy(etareg_d.begin(), etareg_d.end(), std::back_inserter(this->etareg_central));
		  else
		    std::copy(etareg_d.begin(), etareg_d.end(), std::back_inserter(this->etareg));
		}
	    }
	}
      else if(row[0].find("enreg") != std::string::npos)
	{
	  if( row[3]<2 ) {row.bad_row();}
	  vec1d<double_> enreg_d;
	  if(row[0].find(mMethod[Method::BruteForce]) != std::string::npos)
	    {
	      if(row[0].find("mask"+mask) != std::string::npos)
		{
		  uint_ i = 0;
		  while(i < this->nreg)
		    {
		      if(row[0] == "enreg_fineeta_sr"+std::to_string(i+1)+"_mask"+mask)
			{
			  enreg_d.clear();
			  for(uint_ j=1; j<row.size(); ++j)
			    enreg_d.push_back(std::stof(row[j]));
			  check_vector_is_sorted(enreg_d);
			  this->enreg_fineeta[i].resize(enreg_d.size());
			  for(uint_ k=0; k<enreg_d.size(); ++k)
			    this->enreg_fineeta[i][k] = static_cast<float_>(enreg_d[k]);
			}
		      ++i;
		    }
		}
	    }
	  else
	    {
	      if( row.size()!=4 ) {row.bad_row();}
	      enreg_d = linspace(std::stod(row[1]), std::stod(row[2]), std::stoi(row[3]));
	      std::copy(enreg_d.begin(), enreg_d.end(), std::back_inserter(this->enreg));
	    }
	}
      else if(row[0].find("phireg") != std::string::npos)
	{
	  if( row.size()!=4 ) {row.bad_row();}
	  if( row[3]<2 ) {row.bad_row();}
	  vec1d<double_> phireg_d;
	  if(row[0].find(mMethod[Method::BruteForce]) != std::string::npos)
	    {
	      if(row[0].find("mask"+mask) != std::string::npos)
		{
		  uint_ i = 0;
		  while(i < this->nreg)
		    {
		      if(row[0] == "phireg_fineeta_sr"+std::to_string(i+1)+"_mask"+mask)
			{
			  phireg_d.clear();
			  phireg_d = linspace(std::stod(row[1]), std::stod(row[2]), std::stoi(row[3]));
			  std::copy(phireg_d.begin(), phireg_d.end(), std::back_inserter(this->phireg_fineeta[i]));
			  check_vector_is_sorted(phireg_d);
			}
		      ++i;
		    }
		}
	    }
	  else
	    {
	      if( row.size()!=4 ) {row.bad_row();}
	      phireg_d = linspace(std::stod(row[1]), std::stod(row[2]), std::stoi(row[3]));
	      std::copy(phireg_d.begin(), phireg_d.end(), std::back_inserter(this->phireg));
	    }
	}
      else if(row[0] == "nreg")
	{
	  if( row.size()!=2 ) {row.bad_row();}
	  this->nreg = std::stoi(row[1]);
	}
      else if(row[0] == "input")
	{
	  if( row.size()!=2 ) {row.bad_row();}
	  this->noPUFile = row[1]+"mask"+mask+"_"+samples+"_"+particle+"s.root";
	}
      else if(row[0] == "input_raw")
	{
	  if( row.size()!=2 ) {row.bad_row();}
	  this->noPUFile_raw = row[1]+"/"+particle+"s/mask"+mask+"_"+samples+"/hadd_mask"+mask+"_"+samples+"_nopu.root";
	}
      else if(row[0] == "output")
	{
	  if( row.size()!=2 ) {row.bad_row();}
	  this->outpath = row[1]+samples+"/mask"+mask+"_"+particle+"s";
	}
      else if(row[0] == "bckgcuts_"+samples)
	{
	  if( row.size()!=4 ) {row.bad_row();}
	  this->bckgcuts = {std:stof(row[1]),std:stof(row[2]),std:stof(row[3])};
	}
      else if(row[0] == "nlayers")
	{
	  if( row.size()!=2 ) {row.bad_row();}
	  this->nlayers = std::stoi(row[1]);
	}
	}
  for(auto it = varnames.begin(); it != varnames.end(); ++it)
    {
      if(std::find(rows.begin(), rows.end(), *it) == rows.end())
	{
	  std::cout << "Variable " + *it + " is not defined in the configuration file " + fname + "." << std::endl;
	  std::exit(0);
	}
    }
}

void CalibratorInputParameters::check_vector_is_sorted(const vec1d<double_>& v)
{
  if (!std::is_sorted(v.begin(), v.end()))
    {
      std::cout << "The input etareg array is not sorted!" << std::endl;
      std::exit(EXIT_FAILURE);
    }
}

