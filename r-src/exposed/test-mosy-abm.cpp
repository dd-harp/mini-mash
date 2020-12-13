/* --------------------------------------------------------------------------------
#
#   disaggregated malaria transmission model ABM
#   test mosy with constant value for Ih
#   Sean L. Wu (slwu89@berkeley.edu)
#   December 2020
#
-------------------------------------------------------------------------------- */

#include <algorithm>
#include <string>

#include <Rcpp.h>


// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp14)]]


#include "abm/mosquito.hpp"

// [[Rcpp::export]]
Rcpp::NumericMatrix test_mosy(
  const int SV, const int IV, const double IH, const double NH,
  const Rcpp::NumericVector& parameters,
  const double dt, const double tmax
){

  mosquito_pop_uptr mpop = std::make_unique<mosquito_pop>();

  double a = parameters["a"];
  double c = parameters["c"];

  mpop->g = parameters["g"];
  mpop->EIP = parameters["EIP"];
  mpop->lambda = parameters["lambda"];

  int i;
  for(i=0; i<SV; i++){
    mpop->pop.emplace_back(make_mosquito(mpop,0.,'S'));
  }
  for(i=0; i<IV; i++){
    mpop->pop.emplace_back(make_mosquito(mpop,0.,'I'));
  }

  double clock{0.0};
  int tsteps = tmax/dt;
  int j{0};

  Rcpp::Rcout << "--- begin simulation ---\n";
  while(clock < tmax){

    run_mosypop(mpop,clock,dt);

    // compute bloodmeals (H2M bites)
    for(auto& it : mpop->pop){
      if(it->atrisk){
        double tend{clock+dt};
        double risk_t0, risk_t1;
        if(it->shist.back() == 'D'){
          risk_t0 = std::max(clock, it->thist.end()[-2]);
          risk_t1 = it->thist.end()[-1];
        } else if(it->shist.back() == 'S'){
          risk_t0 = std::max(clock,it->thist.back());
          risk_t1 = tend;
        } else {
          printf("illegal state: %c \n", it->shist.back());
          Rcpp::stop("");
        }

        double riskt = risk_t1 - risk_t0;
        assert(riskt <= dt);
        double myrate = a*c*(IH/NH)*riskt;
        double myprob = 1. - exp(-myrate);
        if(R::runif(0.,1.) < myprob){
          double btime = R::runif(risk_t0,risk_t1);
          push_H2M_bite(it,btime);
        }
      }

    }

    mpop->Iv_trace.clear();

    j++;
    if(j % 100 == 0){
      Rcpp::Rcout << " completed TWICE step " << j << " of " << tsteps << "\n";
    }
    clock += dt;

  };

  Rcpp::Rcout << "--- end simulation ---\n";

  mpop->pop.clear();


  // process the actual history
  std::sort(
    mpop->hist.begin(),
    mpop->hist.end(),
    [](const hist_elem& a, const hist_elem& b){
      return a.t < b.t;
    }
  );

  int nn =mpop->hist.size();
  Rcpp::NumericMatrix outagg(nn,4);
  for(int i=0; i<nn; i++){
    outagg.at(i,0) = mpop->hist[i].t;
    outagg.at(i,1) = mpop->hist[i].dS;
    outagg.at(i,2) = mpop->hist[i].dE;
    outagg.at(i,3) = mpop->hist[i].dI;
  };
  return outagg;
};
