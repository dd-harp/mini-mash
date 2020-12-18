/* --------------------------------------------------------------------------------
#
#   disaggregated malaria transmission model: population level description
#   bloodmeal
#   Sean L. Wu (slwu89@berkeley.edu)
#   Dec 2020
#
-------------------------------------------------------------------------------- */

#ifndef BLOODMEAL_HPP
#define BLOODMEAL_HPP

#include "human.hpp"
#include "mosquito.hpp"
#include "shared.hpp"


void bloodmeal(
  mosy_pop_ptr& mpop,
  human_pop_ptr& hpop,
  const double t0,
  const double dt,
  const Rcpp::NumericVector& parameters
);

void bloodmeal_old(
  mosy_pop_ptr& mpop,
  human_pop_ptr& hpop,
  const Rcpp::NumericVector& parameters
);

#endif
