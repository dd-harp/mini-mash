/* --------------------------------------------------------------------------------
#
#   disaggregated malaria transmission model ABM
#   human
#   Sean L. Wu (slwu89@berkeley.edu)
#   December 2020
#
-------------------------------------------------------------------------------- */

#include "human.hpp"

// make a person
human_uptr make_human(human_pop_uptr& hpop, const char state){

  human_uptr hh = std::make_unique<human>();

  hh->tnow = 0.;
  hh->snow = state;

  if(state == 'S'){

    hh->tnext = 0.;
    hh->snext = state;

  } else if(state == 'I'){

    double trec = R::rexp(1./hpop->r);
    hh->tnext = trec;
    hh->snext = 'S';

  } else {
    Rcpp::stop("invalid initial state given to 'make_human'");
  }

  hh->pop = hpop.get();

  return hh;
};




void sim_human_S(human_uptr& hh, const double t0, const double dt){

  double tmax{t0+dt};
  hh->tnow = hh->tnext;

  // I2S transition
  if(hh->snow == 'I'){
    hh->snow = 'S';
    hh->pop->hist.emplace_back(hist_elem{hh->tnow,1,0,-1});
  }

  // queue S2S (so that we can accumulate infection hazard on each step properly)
  hh->tnext = tmax + 2E-8;

};


void sim_human_E(human_uptr& hh, const double t0, const double dt){

  // S2E transition
  hh->tnow = hh->tnext;
  hh->snow = 'E';
  hh->pop->hist.emplace_back(hist_elem{hh->tnow,-1,1,0});

  // queue next transition (E2I)
  double E2I = hh->tnow + hh->pop->LEP;
  hh->tnext = E2I;
  hh->snext = 'I';

};

void sim_human_I(human_uptr& hh, const double t0, const double dt){

  // E2I transition
  hh->tnow = hh->tnext;
  hh->snow = 'I';
  hh->pop->hist.emplace_back(hist_elem{hh->tnow,0,-1,1});

  // queue next transition (I2S)
  double I2S = hh->tnow + R::rexp(1./hh->pop->r);
  hh->tnext = I2S;
  hh->snext = 'S';

};

// call this from bloodmeal
void push_M2H_bites(human_uptr& hh, const double btime){
  hh->tnext = btime;
  hh->snext = 'E';
};
