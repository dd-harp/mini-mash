/* --------------------------------------------------------------------------------
#
#   disaggregated malaria transmission model (population representation)
#   Sean L. Wu (slwu89@berkeley.edu)
#   December 2020
#
-------------------------------------------------------------------------------- */

// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp14)]]

#include "pop/human.hpp"
#include "pop/mosquito.hpp"
#include "pop/bloodmeal.hpp"
#include "pop/shared.hpp"

// [[Rcpp::export]]
Rcpp::List run_miniMASH_pop(
  const int SV, const int IV,
  const int SH, const int IH,
  const Rcpp::NumericVector& parameters,
  const double dt, const double tmax
){

  // make populations
  human_pop_ptr hpop = make_human_pop(SH,IH,parameters);
  mosy_pop_ptr mpop = make_mosy_pop(SV,IV,parameters);

  // run simulation
  double clock{0.0};
  int tsteps = tmax/dt;
  int j{0};

  Rcpp::Rcout << "--- begin simulation ---\n";
  while(clock < tmax){

    run_human_pop(hpop,clock,dt);
    run_mosy_pop(mpop,clock,dt);

    bloodmeal(mpop,hpop,clock,dt,parameters);

    j++;
    if(j % 100 == 0){
      Rcpp::Rcout << " completed TWICE step " << j << " of " << tsteps << "\n";
    }
    clock += dt;

  };

  Rcpp::Rcout << "--- end simulation ---\n";

  Rcpp::NumericMatrix hhist = gethist_human_pop(hpop,SH,IH);
  Rcpp::NumericMatrix mhist = gethist_mosy_pop(mpop,SV,IV);

  return Rcpp::List::create(
    Rcpp::Named("human") = hhist,
    Rcpp::Named("mosy") = mhist
  );
}
