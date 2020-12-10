/* --------------------------------------------------------------------------------
#
#   disaggregated malaria transmission model ABM
#   mosquito
#   Sean L. Wu (slwu89@berkeley.edu)
#   November 2020
#
-------------------------------------------------------------------------------- */

#include "mosquito.hpp"

// put the functions in a hash table
static const std::unordered_map<char, std::function<void(mosquito_uptr&, const double, const double)> > mosquito_functions  = {
  {'E',sim_mosquito_E},
  {'I',sim_mosquito_I},
  {'D',sim_mosquito_D}
};

// make a skeeter
mosquito_uptr make_mosquito(mosquito_pop_uptr& mpop, const double tnow){

  mosquito_uptr mosy = std::make_unique<mosquito>();

  // at most, the trace will have 4 elements
  mosy->thist.reserve(4);
  mosy->shist.reserve(4);

  mosy->thist.emplace_back(tnow);
  mosy->shist.emplace_back('S');

  mosy->atrisk = true;

  // set pointer
  mosy->pop = mpop.get();

  // sample time of death now
  mosy->tdie = tnow + R::rexp(1./mpop->g);
  mosy->tnext = mosy->tdie;
  mosy->snext = 'D';

  return mosy;
};




void sim_mosquito_pop(mosquito_pop_uptr& mpop, const double t0, const double dt){

  double tmax{t0+dt};

  // remove dead mosquitoes
  mpop->pop.remove_if(
    [](const mosquito_uptr& mm) -> bool {
      return mm->shist.back() == 'D';
    }
  );

  // mosquitoes emerge
  int nborn = R::rpois(mpop->lambda) * dt;
  for(int i=0; i<nborn; i++){
    double bday = R::runif(t0,tmax);
    mpop->pop.emplace_back(make_mosquito(mpop,bday));
  }

  // simulate all mosquitoes over [t0,t0+dt)
  for(auto& mm : mpop->pop){
    sim_mosquito(mm,t0,dt);
  }

};



// simulate a skeeter
void sim_mosquito(mosquito_uptr& mosy, const double t0, const double dt){

  double tmax{t0+dt};

  // update this mosquito over [t0,t0+dt)
  while(mosy->tnext < tmax){
    mosquito_functions.at(mosy->snext)(mosy,t0,dt);
  }

}

// call this from bloodmeal
void push_H2M_bite(mosquito_uptr& mosy, double btime){

  // infection event on a living mosquito
  if(mosy->shist.back() == 'S'){
    mosy->tnext = btime;
    mosy->snext = 'E';

  // infection event on a mosquito that died *after* infection
  } else if(mosy->shist.back() == 'D'){
    auto it_t = mosy->thist.end();
    it_t--;
    mosy->thist.insert(it_t, btime);

    auto it_s = mosy->shist.end();
    it_s--;
    mosy->shist.insert(it_s, 'E');
  } else {
    Rcpp::stop("push_H2M_bite error: mosquito in illegal state E or I");
  }

}








// state machine
void sim_mosquito_E(mosquito_uptr& mosy, const double t0, const double dt){

  // change state
  double tnow = mosy->tnext;
  mosy->thist.emplace_back(mosy->tnext);
  mosy->shist.emplace_back('E');

  mosy->atrisk = false;

  double tinf = tnow + mosy->pop->EIP;
  if(mosy->tdie < tinf){
    mosy->tnext = mosy->tdie;
    mosy->snext = 'D';
  } else {
    mosy->tnext = tinf;
    mosy->snext = 'I';
  }

};

void sim_mosquito_I(mosquito_uptr& mosy, const double t0, const double dt){

  double tmax{t0+dt};

  // if this is the first time entering the I state (E2I) transition
  if(mosy->shist.back() == 'E'){

    // change state
    double tnow = mosy->tnext;
    mosy->thist.emplace_back(tnow);
    mosy->shist.emplace_back('I');

  }

  // next state is death, but we assign I so that this function will get called again
  if(mosy->tdie < tmax){
    mosy->tnext = mosy->tdie;
    mosy->snext = 'D';
  } else {
    mosy->tnext = tmax;
    mosy->snext = 'I';
  }

  // time the mosquito begins contributing to risk on humans
  double risk_begin = std::max(t0,mosy->thist.back());
  double risk_end = std::min(mosy->tdie,tmax);

  
  // accumluated Iv mosquito time to send to bloodmeal
  // double Iv_time = std::min(mosy->tdie,tmax) - std::max(t0,mosy->thist.back());
  // push Iv_time to the mosquito_pop

  // check
  // assert(Iv_time <= dt);

};

void sim_mosquito_D(mosquito_uptr& mosy, const double t0, const double dt){

  if(mosy->shist.back() != 'S'){
    mosy->atrisk = false;
  }

  // change state
  double tnow = mosy->tnext;
  mosy->thist.emplace_back(tnow);
  mosy->shist.emplace_back('D');

  mosy->tnext = infinity;

  // output the history

};
