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

// call this from the main loop
void bloodmeal(
  mosy_pop_ptr& mpop,
  human_pop_ptr& hpop,
  const double t0,
  const double dt,
  const Rcpp::NumericVector& parameters
);

// sample S->E events in mosquitoes
void bloodmeal_H2M(
  mosy_pop_ptr& mpop,
  human_pop_ptr& hpop,
  const double tmax,
  const double a,
  const double c,
  const double NH
);

// sample S->E events in humans
void bloodmeal_M2H(
  mosy_pop_ptr& mpop,
  human_pop_ptr& hpop,
  const double tmax,
  const double a,
  const double b,
  const double NH
);

#endif
