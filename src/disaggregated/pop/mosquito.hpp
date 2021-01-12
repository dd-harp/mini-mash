/* --------------------------------------------------------------------------------
#
#   disaggregated malaria transmission model: aggregated population
#   mosquito
#   Sean L. Wu (slwu89@berkeley.edu)
#   Dec 2020
#
-------------------------------------------------------------------------------- */


#ifndef MOSQUITO_HPP
#define MOSQUITO_HPP

#include <vector>
#include <array>
#include <memory>
#include <algorithm>

#include <Rcpp.h>

#include "shared.hpp"


/* --------------------------------------------------------------------------------
#   mosquito population
-------------------------------------------------------------------------------- */

typedef struct mosy_pop {

  // states
  int S;
  int E;
  int I;

  // history tracking
  std::vector<hist_elem> hist;

  double tnow;

  // bites when H->M transmission occured
  queue H2M_bites;

  // when the resulting infection in M will manifest (E->I)
  queue M_inf;

  // parameters
  double g;
  double EIP;
  double lambda;

  // elements for MNRM of internal dynamics
  std::array<double,5> delta_t{infinity};
  std::array<double,4> Pk{0.};
  std::array<double,4> Tk{0.};
  std::array<double,4> ak{0.};

  // state transitions for bloodmeal over [t0,t0+dt)
  queue_trace Sv_trace;
  queue_trace Iv_trace;

  mosy_pop(){};
  ~mosy_pop(){};

} mosy_pop;

using mosy_pop_ptr = std::unique_ptr<mosy_pop>;


/* --------------------------------------------------------------------------------
#   creation and history extraction
-------------------------------------------------------------------------------- */

mosy_pop_ptr make_mosy_pop(
  const int SV,
  const int IV,
  const Rcpp::NumericVector& parameters
);

Rcpp::NumericMatrix gethist_mosy_pop(
  mosy_pop_ptr& mpop,
  const int SV,
  const int IV
);


/* --------------------------------------------------------------------------------
#   simulate over a time interval
-------------------------------------------------------------------------------- */

void run_mosy_pop(
  mosy_pop_ptr& mpop,
  const double t0,
  const double dt
);


#endif
