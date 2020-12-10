/* --------------------------------------------------------------------------------
#
#   disaggregated malaria transmission model ABM
#   shared data types
#   Sean L. Wu (slwu89@berkeley.edu)
#   December 2020
#
-------------------------------------------------------------------------------- */

#include <limits>

double infinity = std::numeric_limits<double>::infinity();

// history data struct
struct hist_elem {
  double t;
  int S;
  int E;
  int I;
  hist_elem(double t_, int S_, int E_, int I_) : t(t_), S(S_), E(E_), I(I_) {};
};
