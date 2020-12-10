/* --------------------------------------------------------------------------------
#
#   disaggregated malaria transmission model ABM
#   mosquito
#   Sean L. Wu (slwu89@berkeley.edu)
#   November 2020
#
-------------------------------------------------------------------------------- */

#include <list>
#include <vector>
#include <memory>
#include <functional>
#include <unordered_map>

#include <Rcpp.h>

#include "shared.hpp"

// forward declare the mosy pop struct
typedef struct mosquito_pop mosquito_pop;

// a single skeeter (syringe with wings)
// mosquito object MUST have a shorter lifespan than the mosquito_pop one.
typedef struct mosquito {

  int id;
  std::vector<double> thist;
  std::vector<char>   shist;

  bool atrisk; // does the mosy need to accum haz?

  double tnext;
  char snext;

  double tdie; // time I will die

  mosquito_pop* pop; // pointer to population struct

} mosquito;

// they're gonna be unique pointers
using mosquito_uptr = std::unique_ptr<mosquito>;




struct mosquito_pop {

  // the population of living mosquitoes
  std::list<mosquito_uptr> pop;

  // parameters
  double g;
  double EIP;
  double lambda;

  // Iv trace
  double Iv_time;

  // history
  std::vector<hist_elem> hist;

};

// pointer to the pop
using mosquito_pop_uptr = std::unique_ptr<mosquito_pop>;


// make a skeeter
mosquito_uptr make_mosquito(mosquito_pop_uptr& mpop, const double tnow);

void sim_mosquito_pop(mosquito_pop_uptr& mpop, const double t0, const double dt);

// simulate a skeeter
void sim_mosquito(mosquito_uptr& mosy, const double t0, const double dt);

// this is what the bloodmeal module will use to push bites onto Sv
void push_H2M_bite(mosquito_uptr& mosy, double btime);

// state machine
void sim_mosquito_E(mosquito_uptr& mosy, const double t0, const double dt);

void sim_mosquito_I(mosquito_uptr& mosy, const double t0, const double dt);

void sim_mosquito_D(mosquito_uptr& mosy, const double t0, const double dt);
