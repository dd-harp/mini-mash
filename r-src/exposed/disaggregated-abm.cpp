/* --------------------------------------------------------------------------------
#
#   disaggregated malaria transmission model (ABM)
#   Sean L. Wu (slwu89@berkeley.edu)
#   December 2020
#
-------------------------------------------------------------------------------- */

// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp14)]]

#include "abm/human.hpp"
#include "abm/mosquito.hpp"
#include "abm/bloodmeal.hpp"
#include "abm/shared.hpp"

// [[Rcpp::export]]
Rcpp::List run_abm(
  const int SV, const int IV,
  const int SH, const int IH,
  const Rcpp::NumericVector& parameters,
  const double dt, const double tmax
){

  // make humans
  human_pop_uptr hpop = std::make_unique<human_pop>();

  hpop->N = SH + IH;
  hpop->r = parameters["r"];
  hpop->LEP = parameters["LEP"];

  int i;
  for(i=0; i<SH; i++){
    hpop->pop.emplace_back(make_human(hpop,'S'));
  }
  for(i=0; i<IH; i++){
    hpop->pop.emplace_back(make_human(hpop,'I'));
  }

  // make mosquitoes
  mosquito_pop_uptr mpop = std::make_unique<mosquito_pop>();

  mpop->g = parameters["g"];
  mpop->EIP = parameters["EIP"];
  mpop->lambda = parameters["lambda"];

  for(i=0; i<SV; i++){
    mpop->pop.emplace_back(make_mosquito(mpop,0.,'S'));
  }
  for(i=0; i<IV; i++){
    mpop->pop.emplace_back(make_mosquito(mpop,0.,'I'));
  }

  // run simulation
  double clock{0.0};
  int tsteps = tmax/dt;
  int j{0};

  Rcpp::Rcout << "--- begin simulation ---\n";
  while(clock < tmax){

    run_humanpop(hpop,clock,dt);
    run_mosypop(mpop,clock,dt);

    bloodmeal(mpop,hpop,clock,dt,parameters);

    j++;
    if(j % 100 == 0){
      Rcpp::Rcout << " completed TWICE step " << j << " of " << tsteps << "\n";
    }
    clock += dt;

  };

  Rcpp::Rcout << "--- end simulation ---\n";

  // prepare output
  mpop->pop.clear(); // need to call mosy dtor to output state trajectory

  // human hist
  // process the actual history
  std::sort(
    hpop->hist.begin(),
    hpop->hist.end(),
    [](const hist_elem& a, const hist_elem& b){
      return a.t < b.t;
    }
  );

  int nh = hpop->hist.size();
  Rcpp::NumericMatrix hhist(nh,4);
  for(int i=0; i<nh; i++){
    hhist.at(i,0) = hpop->hist[i].t;
    hhist.at(i,1) = hpop->hist[i].dS;
    hhist.at(i,2) = hpop->hist[i].dE;
    hhist.at(i,3) = hpop->hist[i].dI;
  };

  // mosy hist
  std::sort(
    mpop->hist.begin(),
    mpop->hist.end(),
    [](const hist_elem& a, const hist_elem& b){
      return a.t < b.t;
  }
  );

  int nm = mpop->hist.size();
  Rcpp::NumericMatrix mhist(nm,4);
  for(int i=0; i<nm; i++){
    mhist.at(i,0) = mpop->hist[i].t;
    mhist.at(i,1) = mpop->hist[i].dS;
    mhist.at(i,2) = mpop->hist[i].dE;
    mhist.at(i,3) = mpop->hist[i].dI;
  };

  return Rcpp::List::create(
    Rcpp::Named("human") = hhist,
    Rcpp::Named("mosy") = mhist
  );
}
