/* --------------------------------------------------------------------------------
#
#   disaggregated malaria transmission model ABM
#   bloodmeal
#   Sean L. Wu (slwu89@berkeley.edu)
#   Dec 2020
#
-------------------------------------------------------------------------------- */

#include "bloodmeal.hpp"

void bloodmeal(
  mosquito_pop_uptr& mpop,
  human_pop_uptr& hpop,
  const double t0,
  const double dt,
  const Rcpp::NumericVector& parameters
){

  double a = parameters["a"];
  double b = parameters["b"];
  double c = parameters["c"];
  double Nh = hpop->N;

  // compute M2H bloodmeals (human inf: spz enter bloodstream)
  for(auto& hh : hpop->pop){
    if(hh->snow == 'S'){
      // get my risk interval, and the Iv trace over that interval
      boost::icl::interval_map<double, int> Iv_trace(mpop->Iv_trace);
      boost::icl::continuous_interval<double> risk_int = getrisk_human(hh,t0,dt);
      Iv_trace &= risk_int;
      // calculate piecewise risk for my interval
      for(auto it = Iv_trace.begin(); it != Iv_trace.end(); it++){
        double t0 = boost::icl::lower(it->first);
        double t1 = boost::icl::upper(it->first);
        int Iv = it->second;
        double haz = a*b*Iv*(1./Nh)*(t1-t0);
        // sample if bite occured, if so, when?
        if(R::runif(0.,1.) < std::exp(-haz)){
          double btime = R::runif(t0,t1);
          push_M2H_bite(hh,btime);
          break;
        }
      }
    }
  }

  // compute H2M bloodmeals (mosy inf: gametocytes enter mosquito)
  for(auto& mm : mpop->pop){
    if(mm->atrisk){
      // get my risk interval, and the Ih trace over that interval
      boost::icl::interval_map<double, int> Ih_trace(hpop->Ih_trace);
      boost::icl::continuous_interval<double> risk_int = getrisk_mosy(mm,t0,dt);
      Ih_trace &= risk_int;
      // calculate piecewise risk for my interval
      for(auto it = Ih_trace.begin(); it != Ih_trace.end(); it++){
        double t0 = boost::icl::lower(it->first);
        double t1 = boost::icl::upper(it->first);
        int Ih = it->second;
        double haz = a*c*(Ih/Nh)*(t1-t0);
        // sample if bite occured, if so, when?
        if(R::runif(0.,1.) < std::exp(-haz)){
          double btime = R::runif(t0,t1);
          push_H2M_bite(mm,btime);
          break;
        }
      }
    }
  }

  // clear traces
  mpop->Iv_trace.clear();
  hpop->Ih_trace.clear();
};
