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

// mosquito dtor
mosquito::~mosquito(){

  assert(thist.size() >= 2); // at minimum, have {S,D}
  assert(thist.size() == shist.size());

  int i{0};
  if(shist[i] == 'S'){
    pop->hist.emplace_back(hist_elem{thist[i],1,0,0});
  } else {
    pop->hist.emplace_back(hist_elem{thist[i],0,0,1});
  }

  for(i=1; i<thist.size(); i++){

    if(shist[i-1] == 'S'){
      // S2D
      if(shist[i] == 'D'){
        pop->hist.emplace_back(hist_elem{thist[i],-1,0,0});
      // S2E
      } else if(shist[i] == 'E'){
        pop->hist.emplace_back(hist_elem{thist[i],-1,1,0});
      }
    } else if(shist[i-1] == 'E'){
      // E2D
      if(shist[i] == 'D'){
        pop->hist.emplace_back(hist_elem{thist[i],0,-1,0});
      // E2I
      } else if(shist[i] == 'I'){
        pop->hist.emplace_back(hist_elem{thist[i],0,-1,1});
      }
    } else if(shist[i-1] == 'I'){
      assert(shist[i] == 'D');
      pop->hist.emplace_back(hist_elem{thist[i],0,0,-1});
    } else {
      Rcpp::stop("illegal state detected, called from 'mosquito::~mosquito'");
    }
  }
  // Rcpp::Rcout << "DONE printing mosquito trajectory -- \n";
};

// make a skeeter
mosquito_uptr make_mosquito(mosquito_pop_uptr& mpop, const double bday, const char state){

  if(state == 'E'){
    Rcpp::stop("invalid initial state given to 'make_mosquito'");
  }

  mosquito_uptr mosy = std::make_unique<mosquito>();

  // at most, the trace will have 4 elements
  mosy->thist.reserve(4);
  mosy->shist.reserve(4);

  mosy->thist.emplace_back(bday);
  mosy->shist.emplace_back(state);

  // tells us that the bloodmeal needs to sample S2E events for this mosquito
  if(state == 'S'){
    mosy->atrisk = true;
  } else {
    mosy->atrisk = false;
  }

  // set pointer
  mosy->pop = mpop.get();

  // sample time of death now
  mosy->tdie = bday + R::rexp(1./mpop->g);
  mosy->tnext = mosy->tdie;
  mosy->snext = 'D';

  return mosy;
};




void run_mosypop(mosquito_pop_uptr& mpop, const double t0, const double dt){

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
    mpop->pop.emplace_back(make_mosquito(mpop,bday,'S'));
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
void push_H2M_bite(mosquito_uptr& mosy, const double btime){

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
    mosy->tnext = tmax + 2E-8;
    mosy->snext = 'I';
  }

  // push the interval this mosquito contributes to human risk to the trace
  double risk_t0 = std::max(t0,mosy->thist.back());
  double risk_t1 = std::min(mosy->tdie,tmax);
  mosy->pop->Iv_trace += std::make_pair(
    boost::icl::continuous_interval<double>(risk_t0,risk_t1),
    1
  );

};

void sim_mosquito_D(mosquito_uptr& mosy, const double t0, const double dt){

  // S2D
  if(mosy->shist.back() == 'S'){
    mosy->atrisk = true;
  // other 2 D
  } else {
    mosy->atrisk = false;
  }

  // change state
  double tnow = mosy->tnext;
  mosy->thist.emplace_back(tnow);
  mosy->shist.emplace_back('D');

  mosy->tnext = infinity;

};


// calculate the risk interval for this mosquito
boost::icl::continuous_interval<double> getrisk_mosy(
  mosquito_uptr& mosy,
  const double t0,
  const double dt
){
  double tend{t0+dt};
  double risk_t0, risk_t1;
  if(mosy->shist.back() == 'D'){
    risk_t0 = std::max(t0, mosy->thist.end()[-2]);
    risk_t1 = mosy->thist.end()[-1];
  } else if(mosy->shist.back() == 'S'){
    risk_t0 = std::max(t0,mosy->thist.back());
    risk_t1 = tend;
  } else {
    printf("illegal state: %c \n", mosy->shist.back());
    Rcpp::stop("");
  }
  return boost::icl::continuous_interval<double>::right_open(risk_t0,risk_t1);
};
