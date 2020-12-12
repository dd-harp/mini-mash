/* --------------------------------------------------------------------------------
#
#   disaggregated malaria transmission model ABM
#   human
#   Sean L. Wu (slwu89@berkeley.edu)
#   December 2020
#
-------------------------------------------------------------------------------- */

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

// human object MUST have a shorter lifespan than the human_pop one.
typedef struct human {

  double tnow;
  char snow;

  double tnext;
  char snext;

  double I2S;

  human_pop* pop; // pointer to population struct

  human();
  ~human();

  // JUST FOR DEBUGGING
  std::vector<double> thist;
  std::vector<char> shist;

} human;

// they're gonna be unique pointers
using human_uptr = std::unique_ptr<human>;




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

  human_pop();
  ~human_pop();

};

// pointer to the pop
using human_pop_uptr = std::unique_ptr<human_pop>;


void run_humanpop(human_pop_uptr& hpop, const double t0, const double dt);

void sim_human(human_uptr& hh, const double t0, const double dt);

// make a person
human_uptr make_human(human_pop_uptr& hpop, const char state);

void sim_human_S(human_uptr& hh, const double t0, const double dt);

void sim_human_E(human_uptr& hh, const double t0, const double dt);

void sim_human_I(human_uptr& hh, const double t0, const double dt);

// call this from bloodmeal
void push_M2H_bites(human_uptr& hh, const double btime);
