/* --------------------------------------------------------------------------------
#
#   disaggregated malaria transmission model (ABM)
#   Sean L. Wu (slwu89@berkeley.edu)
#   December 2020
#
-------------------------------------------------------------------------------- */

#include <algorithm>
#include <array>

// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp14)]]

#include "abm/human.hpp"
#include "abm/mosquito.hpp"
#include "abm/bloodmeal.hpp"
#include "abm/shared.hpp"

// [[Rcpp::export]]
Rcpp::List run_abm_sumout(
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
  std::sort(
    hpop->hist.begin(),
    hpop->hist.end(),
    [](const hist_elem& a, const hist_elem& b){
      return a.t < b.t;
    }
  );

  auto it_h = std::find_if(
    hpop->hist.begin(),
    hpop->hist.end(),
    [](const hist_elem& elem){
      return elem.t > std::numeric_limits<double>::epsilon();
    }
  );

  int nh = std::distance(it_h,hpop->hist.end());
  Rcpp::NumericMatrix hhist(nh+1,4);
  int h_ix{0};
  hhist.at(h_ix,0) = 0.;
  hhist.at(h_ix,1) = SH;
  hhist.at(h_ix,2) = 0.;
  hhist.at(h_ix,3) = IH;
  h_ix++;

  while(it_h != hpop->hist.end()){

    hhist.at(h_ix,0) = it_h->t;
    hhist.at(h_ix,1) = hhist.at(h_ix-1,1) + it_h->dS;
    hhist.at(h_ix,2) = hhist.at(h_ix-1,2) + it_h->dE;
    hhist.at(h_ix,3) = hhist.at(h_ix-1,3) + it_h->dI;

    h_ix++;
    it_h++;
  }

  // mosy hist
  std::sort(
    mpop->hist.begin(),
    mpop->hist.end(),
    [](const hist_elem& a, const hist_elem& b){
      return a.t < b.t;
  }
  );

  auto it_m = std::find_if(
    mpop->hist.begin(),
    mpop->hist.end(),
    [](const hist_elem& elem){
      return elem.t > std::numeric_limits<double>::epsilon();
    }
  );

  int nm = std::distance(it_m,mpop->hist.end());
  Rcpp::NumericMatrix mhist(nm+1,4);
  int m_ix{0};
  mhist.at(m_ix,0) = 0.;
  mhist.at(m_ix,1) = SV;
  mhist.at(m_ix,2) = 0.;
  mhist.at(m_ix,3) = IV;
  m_ix++;

  while(it_m != mpop->hist.end()){

    mhist.at(m_ix,0) = it_m->t;
    mhist.at(m_ix,1) = mhist.at(m_ix-1,1) + it_m->dS;
    mhist.at(m_ix,2) = mhist.at(m_ix-1,2) + it_m->dE;
    mhist.at(m_ix,3) = mhist.at(m_ix-1,3) + it_m->dI;

    m_ix++;
    it_m++;
  }

  return Rcpp::List::create(
    Rcpp::Named("human") = hhist,
    Rcpp::Named("mosy") = mhist
  );
}
