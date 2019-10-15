#include "UserCode/AnalysisCode/interface/photon_analysis.h"

PhotonsPartialWafersAnalysis::PhotonsPartialWafersAnalysis(std::string mask, std::string samples, std::string method) 
{
  vec1d<std::string> varnames;
  if(method == mMethod[Method::BruteForce])
    varnames = {"mingenen", "etareg_"+samples, "nreg", "input"};
  else if(method == mMethod[Method::ShowerLeakage])
    varnames = {"mingenen", "etareg_"+samples, "nreg", "nlayers", 
		"input", "bckgcuts_"+samples};
  else
    {
      std::cout << "The specified methods was not implemented." << std::endl;
      std::exit(0);
    }
  params = std::make_unique<CalibratorInputParameters>("params_photons.csv", varnames, mask, samples, method, 
						       mParticleType[ParticleType::Photon]);
  calibrator = std::make_unique<Calibrator>(*params);

  if(params->samples == mDetectorRegion[DetectorRegion::Inner]) {
    boundaries = {5, 5, 5};
    corr_mode = "left";
    lshift = {.65, .59, .48};
  }
  else if(params->samples == mDetectorRegion[DetectorRegion::Outer]) {
    boundaries = {23, 23, 23};
    corr_mode = "right";
    lshift = {1., 1., 1.};
  }
}

void PhotonsPartialWafersAnalysis::do_photon_calibration(const vec1d<int_>& nq, const bool_& pu, const bool_& draw_plot,
							 opt<std::string> file_opt, 
							 opt<vec1d<CalibrationType>> calib_types_opt,
							 opt<vec2d<float_>> calib_vars_opt)
{
  this->calibration_values.clear();
  if(calib_types_opt != std::nullopt)
    this->calib_types = calib_types_opt.value();
  else
    this->calib_types = {CalibrationType::GenEta, CalibrationType::RecoEn, CalibrationType::GenPhi};

  if(calib_vars_opt != std::nullopt)
    this->calib_vars = calib_vars_opt.value();
  else
    this->calib_vars = {params->etareg,params->enreg,params->phireg};

  assert(this->calib_types.size() == this->calib_vars.size());
  assert(nq.size() == calib_vars.size());

  std::string file = file_opt.value_or(params->noPUFile);
  (this->calibrator)->create_photon_calibration_values(nq, false, true, file, this->calib_types, this->calib_vars);
  this->calibration_values = (this->calibrator)->calibration_values;
}

std::function<vec1d<float_>(const float_&, const float_&, const vec1d<float_>&)> PhotonsPartialWafersAnalysis::define_f1()
{
  auto lambda = [this](const float_& geta, const float_& gphi, const vec1d<float_>& en)
    {
      vec1d<float_> f1(params->nreg, 1.);
      for(uint_ ireg=0; ireg<params->nreg; ++ireg) 
	{
	  for(auto it_type = this->calib_types.cbegin(); it_type != this->calib_types.cend(); ++it_type) 
	    {
	      uint_ type_idx = std::distance(this->calib_types.cbegin(), it_type);
	      umap<CalibrationType, float_> prod_table = {{CalibrationType::GenEta, geta}, 
							  {CalibrationType::RecoEn, en[ireg]}, 
							  {CalibrationType::GenPhi, gphi}};
	      float vartmp1, vartmp2;
	      if( prod_table[*it_type] <= this->calib_vars[type_idx][0] )
		{
		  vartmp1 = calib_vars[type_idx][0];
		  vartmp2 = calib_vars[type_idx][1];
		}
	      else if( prod_table[*it_type] >= this->calib_vars[type_idx].back() )
		{
		  vartmp1 = calib_vars[type_idx].end()[-2];
		  vartmp2 = calib_vars[type_idx].back();
		}
	      else
		{
		  auto const it_bins = std::lower_bound(calib_vars[type_idx].begin(), calib_vars[type_idx].end(), prod_table[*it_type]);
		  vartmp1 = *(it_bins-1);
		  vartmp2 = *it_bins;
		}
	      std::string str1 = round_to_string(vartmp1, 3); 
	      std::string str2 = round_to_string(vartmp2, 3); 
	      std::string idstr = "sr" + std::to_string(ireg+1) + "from" + str1 + "to" + str2;

	      //it should only get here once per event
	      check_key(idstr, this->calibration_values[to_underlying(*it_type)]);
	      f1[ireg] /= (this->calibration_values[to_underlying(*it_type)][idstr]->Eval(f1[ireg]*prod_table[*it_type])+1.0);
	    }
	  //assert(f1[ireg] != 1.);
	}
      return f1;
    };
  return lambda;
}

std::function<vec1d<float_>(const float_&, const vec1d<float_>&)> PhotonsPartialWafersAnalysis::define_f2()
{
  auto lambda = [this](const float_& geta, const vec1d<float_>& noi)
    {
      vec1d<float_> f2(params->nreg, 0.);
      /*
      for(uint_ ireg=0; ireg<params->nreg; ++ireg) {
	for(auto it = this->etas.cbegin(); it!=(this->etas.cend()-1); ++it) {
	  std::string idstr;
	  if (geta < this->etas[0] or geta > this->etas[this->netas-1]) 
	    {
	      idstr = "sr" + std::to_string(ireg+1) + "from" + 
		etastr(std::to_string(this->etas[0])) + "to" + etastr(std::to_string(this->etas[this->netas-1]));
	    }
	  else if (geta < *it or geta > *(it+1))
	    continue;
	  else 
	    {
	      idstr = "sr" + std::to_string(ireg+1) + "from" + 
		etastr(std::to_string(*it)) + "to" + etastr(std::to_string(*(it+1)));
	    }

	  //it should only get here once per event
	  if (calibration_values[to_underlying(CalibrationType::PU)].size() > 0) 
	    {
	      check_key(idstr, calibration_values[to_underlying(CalibrationType::PU)]);
	      f2[ireg-1] = calibration_values[to_underlying(CalibrationType::PU)][idstr]->Eval(noi[2]);
	    }
	}
      }
      */
      return f2;
    };
  return lambda;
}

std::function<vec1d<float_>(const vec1d<float_>&, const vec1d<float_>&, const vec1d<float_>&)> PhotonsPartialWafersAnalysis::define_calibrated_energy() 
{
  auto lambda = [this](const vec1d<float_>& en, const vec1d<float_>& f1, const vec1d<float_>& f2)
    {
      vec1d<float_> encalib(params->nreg, 0.);
      for(uint_ ireg=0; ireg<params->nreg; ++ireg)
	encalib[ireg] = f1[ireg]*en[ireg] - f2[ireg];
      return encalib;
    };
  return lambda;
}

std::function<bool_(const float_&, const vec1d<float_>&)> PhotonsPartialWafersAnalysis::complete_showers_filter()
{
  auto lambda = [this](const float_& geta, const vec1d<float_>& encalib)
    {
      bool_ selection = ( ((params->samples == mDetectorRegion[DetectorRegion::Inner] && geta < params->etareg[0]+0.05) or 
			   (params->samples == mDetectorRegion[DetectorRegion::Outer] && geta > params->etareg.back()-0.05))
			  and encalib[0] != 0 );
      return selection;
    };
  return lambda;
}

std::function<vec2d<float_>(const vec2d<float_>&, const vec1d<float_>&, const vec1d<float_>&, const vec1d<float_>&)> PhotonsPartialWafersAnalysis::complete_showers_energy_profile() 
{
  auto lambda = [this](const vec2d<float_>& encalib_layer, const vec1d<float_>& encalib, const vec1d<float_>& f1, const vec1d<float_>& f2)
    {
      vec2d<float_> frac_en(params->nreg, vec1d<float_>(params->nlayers, 0.));
      for(uint_ ireg=0; ireg<params->nreg; ++ireg) 
	  for(uint_ il=0; il<params->nlayers; ++il) 
	      frac_en[ireg][il] = (encalib_layer[ireg][il] / encalib[ireg]);
      return frac_en;
    };
  return lambda;
}

std::function<vec2d<float_>(const vec2d<float_>&, const vec1d<float_>&, const vec1d<float_>&)> PhotonsPartialWafersAnalysis::define_calibrated_energy_per_layer() 
{
  auto lambda = [this](const vec2d<float_>& en_layer, const vec1d<float_>& f1, const vec1d<float_>& f2)
    {
      vec2d<float_> en_calib_layer(params->nreg, vec1d<float_>(params->nlayers, 0.));
      for(uint_ ireg=0; ireg<params->nreg; ++ireg) 
	  for(uint_ il=0; il<params->nlayers; ++il) 
	      en_calib_layer[ireg][il] = f1[ireg]*en_layer[ireg][il] - f2[ireg];
      return en_calib_layer;
    };
  return lambda;
}

std::function<vec1d<int_>(const vec1d<float_>&, const vec2d<float_>&, const vec2d<float_>&)> PhotonsPartialWafersAnalysis::define_shower_index(const vec2d<float_>& frac_en) 
{
  auto lambda = [this, &frac_en](const vec1d<float_>& encalib, const vec2d<float_>& encalib_layer, const vec2d<float_>& ROI_en)
    {
      vec1d<int_> sid(params->nreg,-1);
      for(uint_ ireg=0; ireg<params->nreg; ++ireg) 
	sid[ireg] = diff_ed(ROI_en[ireg], frac_en[ireg], params->bckgcuts, 0.05);
      return sid;
    };
  return lambda;
}

std::function<vec2d<float_>(const vec1d<float_>&, const vec2d<float_>&)> PhotonsPartialWafersAnalysis::define_shower_energy_profile() 
{
  auto lambda = [this](const vec1d<float_>& encalib, const vec2d<float_>& encalib_layer)
    {
      vec2d<float_> ROI_en(params->nreg, vec1d<float_>(params->nlayers, 0.));
      for(uint_ ireg=0; ireg<params->nreg; ++ireg) 
	{
	  for(uint_ il=0; il<params->nlayers; ++il) 
	    {
	      if(encalib[ireg]==0)
		ROI_en[ireg][il] = 0.;
	      else
		ROI_en[ireg][il] = encalib_layer[ireg][il]/encalib[ireg];
	    }
	}
      return ROI_en;
    };
  return lambda;
}

std::function<vec1d<float_>(const vec1d<float_>&, const vec2d<float_>&, const vec1d<int_>&)> PhotonsPartialWafersAnalysis::calculate_weight_corrected_energy(SoftwareCorrector& sc) 
{
  vector3d<float_> weights = sc.weights;
  auto lambda = [this, weights](const vec1d<float_>& encalib, const vec2d<float_>& encalib_layer, const vec1d<int_>& sid) 
    {
      vec1d<float_> encalib_corr(params->nreg, 0.);
      for(uint_ ireg=0; ireg<params->nreg; ++ireg) 
	{
	  assert(params->bckgcuts.size()==lshift.size());

	  bool_ weight_limit;
	  for (uint_ il=1; il<=params->nlayers; ++il) 
	    {
	      if(params->samples == mDetectorRegion[DetectorRegion::Inner])
		weight_limit = il > boundaries[ireg];
	      else if(params->samples == mDetectorRegion[DetectorRegion::Outer])
		weight_limit = il < boundaries[ireg];

	      if(sid[ireg]==0)
		encalib_corr[ireg] += encalib_layer[ireg][il-1];
	      else if(sid[ireg]>0 && sid[ireg]<=static_cast<int_>(params->bckgcuts.size())) 
		{
		  int_ w = sid[ireg]-1;
		  if(weights(ireg,w,il-1) != 0 && weight_limit) 
		    {
		      int_ il_shift = static_cast<int_>(std::round(il*lshift[w])-1);
		      if(il_shift == -1) 
			il_shift = 0;
		      encalib_corr[ireg] += encalib_layer[ireg][il-1]/weights(ireg,w,il_shift);
		    }
		}
	      else 
		{
		  std::cout << "There is a problem with the showerid." << std::endl;
		  std::exit(0);
		}
	    }
	}      
      return encalib_corr;
    };
  return lambda;
}

std::function<vec1d<float_>(const float_&, const vec1d<float_>&)> PhotonsPartialWafersAnalysis::calculate_response()
{
  auto lambda = [this](const float_& gen, const vec1d<float_>& encalib) 
    {
      vec1d<float_> deltaE(params->nreg, -2.);
      for(uint_ ireg=0; ireg<params->nreg; ++ireg)
	deltaE[ireg] = encalib[ireg]/gen - 1.;
      return deltaE;
    };
  return lambda;
}

std::function<vec1d<float_>(const vec1d<float_>&, const vec1d<int_>&)> PhotonsPartialWafersAnalysis::apply_low_stats_factor(SoftwareCorrector& sc)
{
  vec1d<float_> low_stats_factor = sc.low_stats_factor(this->boundaries, this->corr_mode);
  auto lambda = [this, low_stats_factor](const vec1d<float_>& encorr, const vec1d<int_>& sid) 
     {
       vec1d<float_> encalib_corr = encorr;
       for(uint_ ireg=0; ireg<params->nreg; ++ireg) 
	 {
	   if(sid[ireg] > 0) 
	     {
	       encalib_corr[ireg] *= ( 1 / (1-low_stats_factor[ireg]) );
	       if(params->samples == mDetectorRegion[DetectorRegion::Inner]) 
		 encalib_corr[ireg] *= 1/0.09;
	       else if(params->samples == mDetectorRegion[DetectorRegion::Outer]) 
		 encalib_corr[ireg] *= 1/0.08;
	     }
	 }
       return encalib_corr;
     };
   return lambda;
}


const CalibratorInputParameters& PhotonsPartialWafersAnalysis::get_calibrator_parameters() const
{
  return *(this->params);
}

const vec1d<float_>& PhotonsPartialWafersAnalysis::get_software_correction_shift() const
{
  
  return lshift;
}
