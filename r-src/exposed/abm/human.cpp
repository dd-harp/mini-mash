/* --------------------------------------------------------------------------------
#
#   disaggregated malaria transmission model ABM
#   human
#   Sean L. Wu (slwu89@berkeley.edu)
#   December 2020
#
-------------------------------------------------------------------------------- */

#include "human.hpp"

// ctor and dtor
human::human(){
  // Rcpp::Rcout << "calling human ctor at " << this << "\n";
}

human::~human(){
  // Rcpp::Rcout << "calling human dtor at " << this << "\n";
}

// ctor and dtor
human_pop::human_pop(){
  Rcpp::Rcout << "calling human_pop ctor at " << this << "\n";
}

human_pop::~human_pop(){
  Rcpp::Rcout << "calling human_pop dtor at " << this << "\n";
}

// simulate the population
void run_humanpop(human_pop_uptr& hpop, const double t0, const double dt){
  // simulate all humans over [t0,t0+dt)
  for(auto& hh : hpop->pop){
    sim_human(hh,t0,dt);
  }

};

// simulate 1 human
void sim_human(human_uptr& hh, const double t0, const double dt){

  double tmax{t0+dt};

  while(hh->tnext < tmax){
    if(hh->snext == 'S'){
      sim_human_S(hh,t0,dt);
    } else if(hh->snext == 'E'){
      sim_human_E(hh,t0,dt);
    } else if(hh->snext == 'I'){
      sim_human_I(hh,t0,dt);
    } else {
      Rcpp::stop("error: illegal human state detected from sim_human");
    }
  }
};


// make a person
human_uptr make_human(human_pop_uptr& hpop, const char state){

  human_uptr hh = std::make_unique<human>();

  hh->tnow = 0.;
  hh->snow = state;

  hh->pop = hpop.get();

  if(state == 'S'){

    hh->I2S = std::nan("");
    hh->tnext = 0.;
    hh->snext = state;

    hh->pop->hist.emplace_back(hist_elem{hh->tnow,1,0,0});

  } else if(state == 'I'){

    hh->I2S = R::rexp(1./hpop->r);
    hh->tnext = 0.;
    hh->snext = state;

    hh->pop->hist.emplace_back(hist_elem{hh->tnow,0,0,1});

  } else {
    Rcpp::stop("invalid initial state given to 'make_human'");
  }

  // JUST FOR DEBUGGING
  hh->thist.push_back(hh->tnow);
  hh->shist.push_back(hh->snow);

  return hh;
};


void sim_human_S(human_uptr& hh, const double t0, const double dt){

  double tmax{t0+dt};
  hh->tnow = hh->tnext;

  // I2S transition
  if(hh->snow == 'I'){
    hh->snow = 'S';
    hh->pop->hist.emplace_back(hist_elem{hh->tnow,1,0,-1});

    // JUST FOR DEBUGGING
    hh->thist.push_back(hh->tnow);
    hh->shist.push_back(hh->snow);
  }

  // queue S2S (so that we can accumulate infection hazard on each step properly)
  hh->tnext = tmax + 2E-8;
  hh->snext = 'S';

};


void sim_human_E(human_uptr& hh, const double t0, const double dt){

  // S2E transition
  hh->tnow = hh->tnext;
  hh->snow = 'E';
  hh->pop->hist.emplace_back(hist_elem{hh->tnow,-1,1,0});

  // JUST FOR DEBUGGING
  hh->thist.push_back(hh->tnow);
  hh->shist.push_back(hh->snow);

  // queue next transition (E2I)
  double E2I = hh->tnow + hh->pop->LEP;
  hh->tnext = E2I;
  hh->snext = 'I';

};

void sim_human_I(human_uptr& hh, const double t0, const double dt){

  double tmax{t0+dt};
  hh->tnow = hh->tnext;
  double risk_t0, risk_t1;

  // E2I transition
  if(hh->snow == 'E'){

    hh->snow = 'I';

    hh->pop->hist.emplace_back(hist_elem{hh->tnow,0,-1,1});
    risk_t0 = hh->tnow;

    hh->I2S = hh->tnow + R::rexp(1./hh->pop->r);

    // JUST FOR DEBUGGING
    hh->thist.push_back(hh->tnow);
    hh->shist.push_back(hh->snow);

  } else {
    risk_t0 = t0;
  }

  // check if we make self-loop (I2I) or not (I2S)
  if(hh->I2S < tmax){
    hh->tnext = hh->I2S;
    hh->snext = 'S';
    risk_t1 = hh->I2S;
  } else {
    hh->tnext = tmax + 2E-8;
    hh->snext = 'I';
    risk_t1 = tmax;
  }

  // push the interval this person contributes to mosquito risk to the trace
  hh->pop->Ih_trace += std::make_pair(
    boost::icl::continuous_interval<double>(risk_t0,risk_t1),
    1
  );

};

// call this from bloodmeal
void push_M2H_bites(human_uptr& hh, const double btime){
  assert(hh->snow == 'S');
  hh->tnext = btime;
  hh->snext = 'E';
};
