/* --------------------------------------------------------------------------------
#
#   disaggregated malaria transmission model ABM
#   human
#   Sean L. Wu (slwu89@berkeley.edu)
#   December 2020
#
-------------------------------------------------------------------------------- */

#include "human.hpp"


/* --------------------------------------------------------------------------------
#   population functions
-------------------------------------------------------------------------------- */

human_pop_uptr make_humanpop(
  const int SH,
  const int IH,
  const Rcpp::NumericVector& parameters
){

  human_pop_uptr hpop = std::make_unique<human_pop>();

  hpop->N = SH + IH;
  hpop->r = parameters["r"];
  hpop->LEP = parameters["LEP"];

  int i;
  for(i=0; i<SH; i++){
    hpop->pop.emplace_back(make_human(hpop,'S'));
  }
  for(i=0; i<IH; i++){
    hpop->pop.emplace_back(make_human(hpop,'I'));
  }

  return hpop;
};


void run_humanpop(
  human_pop_uptr& hpop,
  const double t0,
  const double dt
){

  // simulate all humans over [t0,t0+dt)
  for(auto& hh : hpop->pop){
    sim_human(hh,t0,dt);
  }

};

Rcpp::NumericMatrix gethist_humanpop(
  const int SH,
  const int IH,
  human_pop_uptr& hpop
){

  std::sort(
    hpop->hist.begin(),
    hpop->hist.end(),
    [](const hist_elem& a, const hist_elem& b){
      return a.t < b.t;
    }
  );

  auto it_h = std::find_if(
    hpop->hist.begin(),
    hpop->hist.end(),
    [](const hist_elem& elem){
      return elem.t > std::numeric_limits<double>::epsilon();
    }
  );

  int nh = std::distance(it_h,hpop->hist.end());
  Rcpp::NumericMatrix hhist(nh+1,4);
  Rcpp::colnames(hhist) = Rcpp::CharacterVector::create("time","SH","EH","IH");
  int h_ix{0};
  hhist.at(h_ix,0) = 0.;
  hhist.at(h_ix,1) = SH;
  hhist.at(h_ix,2) = 0.;
  hhist.at(h_ix,3) = IH;
  h_ix++;

  while(it_h != hpop->hist.end()){

    hhist.at(h_ix,0) = it_h->t;
    hhist.at(h_ix,1) = hhist.at(h_ix-1,1) + it_h->dS;
    hhist.at(h_ix,2) = hhist.at(h_ix-1,2) + it_h->dE;
    hhist.at(h_ix,3) = hhist.at(h_ix-1,3) + it_h->dI;

    h_ix++;
    it_h++;
  }

  return hhist;
};


/* --------------------------------------------------------------------------------
#   individual functions
-------------------------------------------------------------------------------- */

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
  hh->snext = 'S';

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

  double tmax{t0+dt};
  hh->tnow = hh->tnext;
  double risk_t0, risk_t1;

  // E2I transition
  if(hh->snow == 'E'){

    hh->snow = 'I';

    hh->pop->hist.emplace_back(hist_elem{hh->tnow,0,-1,1});
    risk_t0 = hh->tnow;

    hh->I2S = hh->tnow + R::rexp(1./hh->pop->r);

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
void push_M2H_bite(human_uptr& hh, const double btime){
  assert(hh->snow == 'S');
  hh->tnext = btime;
  hh->snext = 'E';
};
