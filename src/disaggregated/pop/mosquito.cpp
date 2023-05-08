/* --------------------------------------------------------------------------------
#
#   disaggregated malaria transmission model: aggregated population
#   mosquito
#   Sean L. Wu (slwu89@berkeley.edu)
#   Dec 2020
#
-------------------------------------------------------------------------------- */

#include "mosquito.hpp"


/* --------------------------------------------------------------------------------
#   creation and history extraction
-------------------------------------------------------------------------------- */

mosy_pop_ptr make_mosy_pop(
  const int SV,
  const int EV,
  const int IV,
  const Rcpp::NumericVector& parameters
){

  mosy_pop_ptr mpop = std::make_unique<mosy_pop>();

  mpop->S = SV;
  mpop->E = EV;
  mpop->I = IV;

  mpop->tnow = 0.0;

  mpop->g = parameters["g"];
  mpop->lambda = parameters["lambda"];
  mpop->EIP = parameters["EIP"];

  // calculate initial propensities
  mpop->ak[0] = mpop->lambda;
  mpop->ak[1] = mpop->g * static_cast<double>(mpop->S);
  mpop->ak[2] = mpop->g * static_cast<double>(mpop->E);
  mpop->ak[3] = mpop->g * static_cast<double>(mpop->I);

  // draw initial firing times
  for(int j=0; j<4; j++){
    mpop->Pk[j] = log(1. / R::runif(0.,1.));
  }
  
  // draw initial EV->IV times for delayed individuals
  mpop->M_inf.emplace(infinity);
  for (int i=0; i<EV; i++) {
    double t_inf = R::runif(-mpop->EIP,0.0);
    mpop->M_inf.emplace(t_inf + mpop->EIP);
  }

  return mpop;
};

Rcpp::NumericMatrix gethist_mosy_pop(
  mosy_pop_ptr& mpop,
  const int SV,
  const int EV,
  const int IV
){

  std::sort(
    mpop->hist.begin(),
    mpop->hist.end(),
    [](const hist_elem& a, const hist_elem& b){
      return a.t < b.t;
    }
  );

  int n = mpop->hist.size();

  Rcpp::NumericMatrix histout(n+1,4);
  Rcpp::colnames(histout) = Rcpp::CharacterVector::create("time","SV","EV","IV");
  int ix{0};
  histout.at(ix,0) = 0.;
  histout.at(ix,1) = SV;
  histout.at(ix,2) = EV;
  histout.at(ix,3) = IV;
  ix++;

  for(auto it = mpop->hist.begin(); it != mpop->hist.end(); it++){
    histout.at(ix,0) = it->t;
    histout.at(ix,1) = histout.at(ix-1,1) + it->dS;
    histout.at(ix,2) = histout.at(ix-1,2) + it->dE;
    histout.at(ix,3) = histout.at(ix-1,3) + it->dI;
    ix++;
  }

  return histout;
};


/* --------------------------------------------------------------------------------
#   simulate over a time interval
-------------------------------------------------------------------------------- */

void run_mosy_pop(
  mosy_pop_ptr& mpop,
  const double t0,
  const double dt
){

  double tmax{t0+dt};

  // get bites from the previous time step [t0-dt,t0)
  while(!mpop->H2M_bites.empty()){

    // time of that bite (S->E)
    double btime = *mpop->H2M_bites.begin();

    // find out how much hazard they accumulate between [btime,t0)
    double haz = (t0 - btime) * mpop->g;

    // did they survive until t0? if not, they died between [btime,t0)
    if(R::runif(0.,1.) < std::exp(-haz)){

      // S->E
      mpop->S -= 1;
      mpop->E += 1;

      // queue the future EIP completion event (E->I), which will happen on this time step (or next)
      mpop->M_inf.emplace(btime + mpop->EIP);

      // track state change
      mpop->hist.emplace_back(hist_elem{btime,-1,1,0});
    }

    // pop that bite off
    mpop->H2M_bites.erase(mpop->H2M_bites.begin());
  }

  assert(mpop->E == mpop->M_inf.size()-1); // -1 because always has Inf

  // push initial state into traces
  mpop->Sv_trace.emplace(queue_tuple(mpop->S,t0));
  mpop->Iv_trace.emplace(queue_tuple(mpop->I,t0));

  // recalculate propentities valid at t0
  mpop->ak[1] = mpop->g * static_cast<double>(mpop->S);
  mpop->ak[2] = mpop->g * static_cast<double>(mpop->E);
  mpop->ak[3] = mpop->g * static_cast<double>(mpop->I);

  // simulate dynamics over this time step [t0,tmax)
  while(mpop->tnow < tmax){

    // absolute next firing times
    for(int j=0; j<4; j++){
      mpop->delta_t[j] = (mpop->Pk[j] - mpop->Tk[j]) / mpop->ak[j];
    }
    mpop->delta_t[4] = *mpop->M_inf.begin() - mpop->tnow;

    // find minimum
    auto min_elem = std::min_element(mpop->delta_t.begin(), mpop->delta_t.end());
    int mu = std::distance(mpop->delta_t.begin(), min_elem);

    // update putative time
    double delta = mpop->delta_t[mu];
    double tnext = mpop->tnow + delta;

    // if overrun, advance time to tmax and adjust the Poisson process for the fast forward
    if(tnext > tmax){
      double remaining = tmax - mpop->tnow;
      for(int j=0; j<4; j++){
        mpop->Tk[j] += mpop->ak[j] * remaining;
      }
      mpop->tnow = tmax;
      break;
    }

    mpop->tnow = tnext;

    // S emerges
    if(mu == 0){
      mpop->S += 1;
      // update history and trace
      mpop->hist.emplace_back(hist_elem{mpop->tnow,1,0,0});
      mpop->Sv_trace.emplace(queue_tuple(mpop->S,mpop->tnow));
    // S dies
    } else if(mu == 1){
      mpop->S -= 1;
      // update history and trace
      mpop->hist.emplace_back(hist_elem{mpop->tnow,-1,0,0});
      mpop->Sv_trace.emplace(queue_tuple(mpop->S,mpop->tnow));
    // E dies
    } else if(mu == 2){
      mpop->E -= 1;
      // delete a random element from the EIP set (not the end, which is Inf)
      int relem = R::runif(0.0, mpop->M_inf.size()-1);
      auto it = mpop->M_inf.begin();
      std::advance(it, relem);
      mpop->M_inf.erase(it);
      // update trace
      mpop->hist.emplace_back(hist_elem{mpop->tnow,0,-1,0});
    // I dies
    } else if(mu == 3){
      mpop->I -= 1;
      // update history and trace
      mpop->hist.emplace_back(hist_elem{mpop->tnow,0,0,-1});
      mpop->Iv_trace.emplace(queue_tuple(mpop->I,mpop->tnow));
    // E->I completion
    } else if(mu == 4){
      mpop->E -= 1;
      mpop->I += 1;
      mpop->M_inf.erase(mpop->M_inf.begin());
      // update history and trace
      mpop->hist.emplace_back(hist_elem{mpop->tnow,0,-1,1});
      mpop->Iv_trace.emplace(queue_tuple(mpop->I,mpop->tnow));
    } else {
      std::string msg("invalid minimum element in 'run_mosy_pop': " + std::to_string(mu));
      Rcpp::stop(msg);
    }

    // update Tk
    for(int j=0; j<4; j++){
      mpop->Tk[j] += mpop->ak[j] * delta;
    }

    // update Pk[mu]
    if(mu != 4){
      mpop->Pk[mu] += log(1. / R::runif(0.,1.));
    }

    // recalculate propentities
    mpop->ak[1] = mpop->g * static_cast<double>(mpop->S);
    mpop->ak[2] = mpop->g * static_cast<double>(mpop->E);
    mpop->ak[3] = mpop->g * static_cast<double>(mpop->I);

  }

};
