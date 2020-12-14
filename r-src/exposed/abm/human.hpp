/* --------------------------------------------------------------------------------
#
#   disaggregated malaria transmission model ABM
#   human
#   Sean L. Wu (slwu89@berkeley.edu)
#   December 2020
#
-------------------------------------------------------------------------------- */

#ifndef HUMAN_HPP
#define HUMAN_HPP

#include <list>
#include <vector>
#include <memory>
#include <functional>
#include <unordered_map>

#include <boost/icl/interval_map.hpp>

#include <Rcpp.h>

#include "shared.hpp"

// forward declare the mosy pop struct
typedef struct human_pop human_pop;


/* --------------------------------------------------------------------------------
#   individual human
-------------------------------------------------------------------------------- */

typedef struct human {

  double tnow;
  char snow;

  double tnext;
  char snext;

  double I2S;

  human_pop* pop; // pointer to population struct

  human(){};
  ~human(){};

} human;

using human_uptr = std::unique_ptr<human>;


/* --------------------------------------------------------------------------------
#   population
-------------------------------------------------------------------------------- */

struct human_pop {

  // the population of living humans
  std::list<human_uptr> pop;

  double N;

  // parameters
  double r;
  double LEP;

  // Iv trace
  boost::icl::interval_map<double, int> Ih_trace;

  // history
  std::vector<hist_elem> hist;

  human_pop(){};
  ~human_pop(){};

};

using human_pop_uptr = std::unique_ptr<human_pop>;


/* --------------------------------------------------------------------------------
#   population functions
-------------------------------------------------------------------------------- */

human_pop_uptr make_humanpop(const int SH, const int IH, const Rcpp::NumericVector& parameters);

void run_humanpop(human_pop_uptr& hpop, const double t0, const double dt);

Rcpp::NumericMatrix gethist_humanpop(const int SH, const int IH, human_pop_uptr& hpop);


/* --------------------------------------------------------------------------------
#   individual functions
-------------------------------------------------------------------------------- */

human_uptr make_human(human_pop_uptr& hpop, const char state);

void sim_human(human_uptr& hh, const double t0, const double dt);

void sim_human_S(human_uptr& hh, const double t0, const double dt);

void sim_human_E(human_uptr& hh, const double t0, const double dt);

void sim_human_I(human_uptr& hh, const double t0, const double dt);

// calculate the risk interval for this human
boost::icl::continuous_interval<double> getrisk_human(
  human_uptr& hh,
  const double t0,
  const double dt
);

// call this from bloodmeal
void push_M2H_bite(human_uptr& hh, const double btime);

#endif
