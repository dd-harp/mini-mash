/* --------------------------------------------------------------------------------
#
#   disaggregated malaria transmission model ABM
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
  mosquito_pop_uptr& mpop,
  human_pop_uptr& hpop,
  const double t0,
  const double dt,
  const Rcpp::NumericVector& parameters
);

#endif
