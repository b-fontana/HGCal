#include "UserCode/AnalysisCode/interface/photon_analysis.h"

PhotonsPartialWafersAnalysis::PhotonsPartialWafersAnalysis(std::string mask, std::string samples, std::string method) 
{
  vec1d<std::string> varnames;
  if(method == "fineeta")
    varnames = {"mingenen", "etareg_"+samples, "nreg", "input"};
  else if(method == "ed")
    varnames = {"mingenen", "etareg_"+samples, "nreg", "nlayers", 
		"input", "bckgcuts_"+samples};
  else
    {
      std::cout << "The specified methods was not implemented." << std::endl;
      std::exit(0);
    }
  params =  std::make_unique<CalibratorInputParameters>("params_photons.csv", varnames, mask, samples, "Photon");
  calibrator = std::make_unique<Calibrator>(*params);
  ROOT::EnableImplicitMT(std::thread::hardware_concurrency());
  netas = params->etareg.size();
}

void PhotonsPartialWafersAnalysis::create_photon_calibration_values(const int& nq, const bool_& pu, const bool_& draw_plot)
{
  (this->calibrator)->create_photon_calibration_values(6, false, true);
  this->calibration_values = (this->calibrator)->calibration_values;
}

std::function<vec1d<float_>(const float_&, const vec1d<float_>&)> PhotonsPartialWafersAnalysis::define_f1()
{
  std::function<vec1d<float_>(const float_&, const vec1d<float_>&)> lambda = [this](const float_& geta, const vec1d<float_>& en)
    {
      vec1d<float_> f1(params->nreg, 1.);
      for(uint_ ireg=0; ireg<params->nreg; ++ireg) {
	for(auto it = params->etareg.cbegin(); it!= params->etareg.cend()-1; ++it) {
	  std::string idstr;
	  if (geta < params->etareg[0] or geta > params->etareg[this->netas-1]) 
	    {
	      idstr = "sr" + std::to_string(ireg+1) + "from" + 
		etastr(std::to_string(params->etareg[0])) + "to" + etastr(std::to_string(params->etareg[this->netas-1]));
	    }
	  else if (geta<*it or geta>*(it+1))
	    continue;
	  else 
	    {
	      idstr = "sr" + std::to_string(ireg+1) + "from" + 
		etastr(std::to_string(*it)) + "to" + etastr(std::to_string(*(it+1)));
	    }

	  //it should only get here once per event
	  if (calibration_values[0].size() > 0) 
	    {
	      check_key(idstr, calibration_values[0]);
	      f1[ireg] /= calibration_values[0][idstr]->Eval(geta)+1.0;
	      if (calibration_values[1].size() > 0) 
		{
		  check_key(idstr, calibration_values[1]);
		  f1[ireg] /= calibration_values[1][idstr]->Eval(f1[ireg]*en[ireg])+1.0;
		}
	    }
	}
	assert(f1[ireg] != 1.);
      }
      return f1;
    };
  return lambda;
}

std::function<vec1d<float_>(const float_&, const vec1d<float_>&)> PhotonsPartialWafersAnalysis::define_f2()
{
  std::function<vec1d<float_>(const float_&, const vec1d<float_>&)> lambda = [this](const float_& geta, const vec1d<float_>& noi)
    {
      vec1d<float_> f2(params->nreg, 0.);
      for(uint_ ireg=0; ireg<params->nreg; ++ireg) {
	for(auto it = params->etareg.cbegin(); it!=(params->etareg.cend()-1); ++it) {
	  std::string idstr;
	  if (geta < params->etareg[0] or geta > params->etareg[this->netas-1]) 
	    {
	      idstr = "sr" + std::to_string(ireg+1) + "from" + 
		etastr(std::to_string(params->etareg[0])) + "to" + etastr(std::to_string(params->etareg[this->netas-1]));
	    }
	  else if (geta < *it or geta > *(it+1))
	    continue;
	  else 
	    {
	      idstr = "sr" + std::to_string(ireg+1) + "from" + 
		etastr(std::to_string(*it)) + "to" + etastr(std::to_string(*(it+1)));
	    }

	  //it should only get here once per event
	  if (calibration_values[2].size() > 0) 
	    {
	      check_key(idstr, calibration_values[2]);
	      f2[ireg-1] = calibration_values[2][idstr]->Eval(noi[2]);
	    }
	}
      }
      return f2;
    };
  return lambda;
}

std::function<vec1d<float_>(const vec1d<float_>&, const vec1d<float_>&, const vec1d<float_>&)> PhotonsPartialWafersAnalysis::define_calibrated_energy() 
{
  std::function<vec1d<float_>(const vec1d<float_>&, const vec1d<float_>&, const vec1d<float_>&)> lambda = [this](const vec1d<float_>& en, const vec1d<float_>& f1, const vec1d<float_>& f2)
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
  std::function<bool_(const float_&, const vec1d<float_>&)> lambda = [this](const float_& geta, const vec1d<float_>& encalib)
    {
      bool_ selection = ( ((params->samples == "inner" && geta < params->etareg[0]+0.05) or 
			   (params->samples == "outer" && geta > params->etareg[params->etareg.size()-1]-0.05))
			  and encalib[0] != 0 );
      return selection;
    };
  return lambda;
}

std::function<vec2d<float_>(const vec2d<float_>&, const vec1d<float_>&, const vec1d<float_>&, const vec1d<float_>&)> PhotonsPartialWafersAnalysis::complete_showers_energy_profile() 
{
  std::function<vec2d<float_>(const vec2d<float_>&, const vec1d<float_>&, const vec1d<float_>&, const vec1d<float_>&)> lambda = [this](const vec2d<float_>& encalib_layer, const vec1d<float_>& encalib, const vec1d<float_>& f1, const vec1d<float_>& f2)
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
  std::function<vec2d<float_>(const vec2d<float_>&, const vec1d<float_>&, const vec1d<float_>&)> lambda = [this](const vec2d<float_>& en_layer, const vec1d<float_>& f1, const vec1d<float_>& f2)
    {
      vec2d<float_> en_calib_layer(params->nreg, vec1d<float_>(params->nlayers, 0.));
      for(uint_ ireg=0; ireg<params->nreg; ++ireg) 
	  for(uint_ il=0; il<params->nlayers; ++il) 
	      en_calib_layer[ireg][il] = f1[ireg]*en_layer[ireg][il] - f2[ireg];
      return en_calib_layer;
    };
  return lambda;
}

CalibratorInputParameters PhotonsPartialWafersAnalysis::get_input_parameters() 
{
  return *(this->params);
}
