/* --------------------------------------------------------------------------------
#
#   disaggregated malaria transmission model: aggregated population
#   shared data types
#   Sean L. Wu (slwu89@berkeley.edu)
#   December 2020
#
-------------------------------------------------------------------------------- */

#ifndef SHARED_HPP
#define SHARED_HPP

#include <limits>
#include <set>

#include <boost/icl/interval_map.hpp>

// Inf
static double infinity = std::numeric_limits<double>::infinity();

// store events (earliest/smallest value on top)
using queue = std::set<double,std::less<double>>;

// interval map for piecewise trajectories on [t0,t0+dt)
using interval_map = boost::icl::interval_map<double, int>;

// history data struct
struct hist_elem {
  double t;
  int dS;
  int dE;
  int dI;
  hist_elem(double t_, int dS_, int dE_, int dI_) : t(t_), dS(dS_), dE(dE_), dI(dI_) {};
};

#endif
