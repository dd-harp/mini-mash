/* --------------------------------------------------------------------------------
#
#   disaggregated malaria transmission model ABM
#   mosquito
#   Sean L. Wu (slwu89@berkeley.edu)
#   November 2020
#
-------------------------------------------------------------------------------- */

#ifndef MOSQUITO_HPP
#define MOSQUITO_HPP

#include <list>
#include <vector>
#include <memory>
#include <functional>
#include <unordered_map>

#include <boost/icl/interval_map.hpp>

#include <Rcpp.h>

#include "shared.hpp"

// forward declare the mosy pop struct
typedef struct mosquito_pop mosquito_pop;


/* --------------------------------------------------------------------------------
#   syringe with wings
-------------------------------------------------------------------------------- */

typedef struct mosquito {

  std::vector<double> thist;
  std::vector<char>   shist;

  bool atrisk; // does the mosy need to accum haz?

  double tnext;
  char snext;

  double tdie; // time I will die

  mosquito_pop* pop; // pointer to population struct

  ~mosquito();

} mosquito;

using mosquito_uptr = std::unique_ptr<mosquito>;


/* --------------------------------------------------------------------------------
#   all the skeeters
-------------------------------------------------------------------------------- */

struct mosquito_pop {

  // the population of living mosquitoes
  std::list<mosquito_uptr> pop;

  // parameters
  double g;
  double EIP;
  double lambda;

  // for bloodmeal
  double a;
  double b;
  double c;

  // Iv trace
  boost::icl::interval_map<double, int> Iv_trace;

  // history
  std::vector<hist_elem> hist;

};

using mosquito_pop_uptr = std::unique_ptr<mosquito_pop>;


/* --------------------------------------------------------------------------------
#   population functions
-------------------------------------------------------------------------------- */

mosquito_pop_uptr make_mosypop(const int SV, const int IV, const Rcpp::NumericVector& parameters);

void run_mosypop(mosquito_pop_uptr& mpop, const double t0, const double dt);

Rcpp::NumericMatrix gethist_mosypop(const int SV, const int IV, mosquito_pop_uptr& mpop);


/* --------------------------------------------------------------------------------
#   individual functions
-------------------------------------------------------------------------------- */

// make a skeeter
mosquito_uptr make_mosquito(mosquito_pop_uptr& mpop, const double bday, const char state);

// simulate a skeeter
void sim_mosquito(mosquito_uptr& mosy, const double t0, const double dt);

void sim_mosquito_E(mosquito_uptr& mosy, const double t0, const double dt);

void sim_mosquito_I(mosquito_uptr& mosy, const double t0, const double dt);

void sim_mosquito_D(mosquito_uptr& mosy, const double t0, const double dt);

// calculate the risk interval for this mosquito
boost::icl::continuous_interval<double> getrisk_mosy(
  mosquito_uptr& mosy,
  const double t0,
  const double dt
);

// this is what the bloodmeal module will use to push bites onto Sv
void push_H2M_bite(mosquito_uptr& mosy, const double btime);

#endif
