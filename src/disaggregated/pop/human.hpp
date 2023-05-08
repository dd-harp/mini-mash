/* --------------------------------------------------------------------------------
#
#   disaggregated malaria transmission model: aggregated population
#   human
#   Sean L. Wu (slwu89@berkeley.edu)
#   November 2020
#
-------------------------------------------------------------------------------- */

#ifndef HUMAN_HPP
#define HUMAN_HPP

#include <vector>
#include <array>
#include <memory>
#include <algorithm>

#include <Rcpp.h>

#include "shared.hpp"


/* --------------------------------------------------------------------------------
#   human population
-------------------------------------------------------------------------------- */

typedef struct human_pop {

  // states
  int S;
  int E;
  int I;

  // history tracking
  std::vector<hist_elem> hist;

  double tnow;

  // bites when M->H transmission occurred
  queue M2H_bites;

  // when the resulting infection in M will manifest (E->I)
  queue H_inf;

  // parameters
  double r;
  double LEP;
  double N;

  // elements for MNRM of internal dynamics
  std::array<double,2> delta_t{infinity};
  double Pk{0.};
  double Tk{0.};
  double ak{0.};

  // state transitions for bloodmeal over [t0,t0+dt)
  queue_trace Sh_trace;
  queue_trace Ih_trace;

  human_pop(){};
  ~human_pop(){};

} human_pop;

using human_pop_ptr = std::unique_ptr<human_pop>;


/* --------------------------------------------------------------------------------
#   creation and history extraction
-------------------------------------------------------------------------------- */

human_pop_ptr make_human_pop(
  const int SH,
  const int EH,
  const int IH,
  const Rcpp::NumericVector& parameters
);

Rcpp::NumericMatrix gethist_human_pop(
  human_pop_ptr& hpop,
  const int SH,
  const int EH,
  const int IH
);


/* --------------------------------------------------------------------------------
#   simulate over a time interval
-------------------------------------------------------------------------------- */

void run_human_pop(
  human_pop_ptr& hpop,
  const double t0,
  const double dt
);


#endif
