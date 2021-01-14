/* --------------------------------------------------------------------------------
#
#   disaggregated malaria transmission model ABM
#   shared data types
#   Sean L. Wu (slwu89@berkeley.edu)
#   December 2020
#
-------------------------------------------------------------------------------- */

#ifndef SHARED_HPP
#define SHARED_HPP

#include <limits>

static double infinity = std::numeric_limits<double>::infinity();
static double epsilon = std::numeric_limits<double>::epsilon();

// history data struct
struct hist_elem {
  double t;
  int dS;
  int dE;
  int dI;
  hist_elem(double t_, int dS_, int dE_, int dI_) : t(t_), dS(dS_), dE(dE_), dI(dI_) {};
};

#endif
