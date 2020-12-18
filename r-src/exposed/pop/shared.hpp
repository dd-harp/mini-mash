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

// Inf
static double infinity = std::numeric_limits<double>::infinity();

// store events (earliest/smallest value on top)
using queue = std::set<double,std::less<double>>;

// data structure to pass traces
using queue_tuple = std::tuple<double,double>; // 1st element is state, 2nd is time

struct queue_comp {
  bool operator() (const queue_tuple& a, const queue_tuple& b) const {
    return std::get<1>(a) < std::get<1>(b);
  }
};

using queue_trace = std::set<queue_tuple,queue_comp>;


// history data struct
struct hist_elem {
  double t;
  int dS;
  int dE;
  int dI;
  hist_elem(double t_, int dS_, int dE_, int dI_) : t(t_), dS(dS_), dE(dE_), dI(dI_) {};
};

#endif
