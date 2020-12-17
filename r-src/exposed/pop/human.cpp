/* --------------------------------------------------------------------------------
#
#   disaggregated malaria transmission model: aggregated population
#   human
#   Sean L. Wu (slwu89@berkeley.edu)
#   November 2020
#
-------------------------------------------------------------------------------- */

#include "human.hpp"


/* --------------------------------------------------------------------------------
#   creation and history extraction
-------------------------------------------------------------------------------- */

human_pop_ptr make_human_pop(
  const int SH,
  const int IH,
  const Rcpp::NumericVector& parameters
){

  human_pop_ptr hpop = std::make_unique<human_pop>();

  hpop->S = SH;
  hpop->E = 0;
  hpop->I = IH;

  hpop->N = SH + IH;

  hpop->tnow = 0.0;

  hpop->r = parameters["r"];
  hpop->LEP = parameters["LEP"];

  // calculate initial propensities
  hpop->ak = hpop->r * static_cast<double>(hpop->I);

  // draw initial firing times
  hpop->Pk = log(1. / R::runif(0.,1.));

  hpop->H_inf.emplace(infinity);

  return hpop;
};


Rcpp::NumericMatrix gethist_human_pop(
  human_pop_ptr& hpop,
  const int SH,
  const int IH
){

  std::sort(
    hpop->hist.begin(),
    hpop->hist.end(),
    [](const hist_elem& a, const hist_elem& b){
      return a.t < b.t;
  }
  );

  auto it = std::find_if(
    hpop->hist.begin(),
    hpop->hist.end(),
    [](const hist_elem& elem){
      return elem.t > std::numeric_limits<double>::epsilon();
    }
  );

  int nm = std::distance(it,hpop->hist.end());
  Rcpp::NumericMatrix hhist(nm+1,4);
  Rcpp::colnames(hhist) = Rcpp::CharacterVector::create("time","SH","EH","IH");
  int ix{0};
  hhist.at(ix,0) = 0.;
  hhist.at(ix,1) = SH;
  hhist.at(ix,2) = 0.;
  hhist.at(ix,3) = IH;
  ix++;

  while(it != hpop->hist.end()){

    hhist.at(ix,0) = it->t;
    hhist.at(ix,1) = hhist.at(ix-1,1) + it->dS;
    hhist.at(ix,2) = hhist.at(ix-1,2) + it->dE;
    hhist.at(ix,3) = hhist.at(ix-1,3) + it->dI;

    ix++;
    it++;
  }

  return hhist;
};


/* --------------------------------------------------------------------------------
#   simulate over a time interval
-------------------------------------------------------------------------------- */

void run_human_pop(
  human_pop_ptr& hpop,
  const double t0,
  const double dt
);
