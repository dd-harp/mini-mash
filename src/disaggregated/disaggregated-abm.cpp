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
Rcpp::List run_miniMASH_abm(
  const int SV, const int EV, const int IV,
  const int SH, const int EH, const int IH,
  const Rcpp::NumericVector& parameters,
  const double dt, const double tmax
){

  // make populations
  human_pop_uptr hpop = make_humanpop(SH,EH,IH,parameters);
  mosquito_pop_uptr mpop = make_mosypop(SV,EV,IV,parameters);

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

  Rcpp::NumericMatrix hhist = gethist_humanpop(SH,EH,IH,hpop);
  Rcpp::NumericMatrix mhist = gethist_mosypop(SV,EV,IV,mpop);

  return Rcpp::List::create(
    Rcpp::Named("human") = hhist,
    Rcpp::Named("mosy") = mhist
  );
}
