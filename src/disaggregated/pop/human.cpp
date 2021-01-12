/* --------------------------------------------------------------------------------
#
#   disaggregated malaria transmission model: aggregated population
#   human
#   Sean L. Wu (slwu89@berkeley.edu)
#   November 2020
#
-------------------------------------------------------------------------------- */

#include "human.hpp"


/* --------------------------------------------------------------------------------
#   creation and history extraction
-------------------------------------------------------------------------------- */

human_pop_ptr make_human_pop(
  const int SH,
  const int IH,
  const Rcpp::NumericVector& parameters
){

  human_pop_ptr hpop = std::make_unique<human_pop>();

  hpop->S = SH;
  hpop->E = 0;
  hpop->I = IH;

  hpop->N = SH + IH;

  hpop->tnow = 0.0;

  hpop->r = parameters["r"];
  hpop->LEP = parameters["LEP"];

  // calculate initial propensities
  hpop->ak = hpop->r * static_cast<double>(hpop->I);

  // draw initial firing times
  hpop->Pk = log(1. / R::runif(0.,1.));

  hpop->H_inf.emplace(infinity);

  return hpop;
};


Rcpp::NumericMatrix gethist_human_pop(
  human_pop_ptr& hpop,
  const int SH,
  const int IH
){

  std::sort(
    hpop->hist.begin(),
    hpop->hist.end(),
    [](const hist_elem& a, const hist_elem& b){
      return a.t < b.t;
    }
  );

  int n = hpop->hist.size();

  Rcpp::NumericMatrix histout(n+1,4);
  Rcpp::colnames(histout) = Rcpp::CharacterVector::create("time","SH","EH","IH");
  int ix{0};
  histout.at(ix,0) = 0.;
  histout.at(ix,1) = SH;
  histout.at(ix,2) = 0.;
  histout.at(ix,3) = IH;
  ix++;

  for(auto it = hpop->hist.begin(); it != hpop->hist.end(); it++){
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

void run_human_pop(
  human_pop_ptr& hpop,
  const double t0,
  const double dt
){

  double tmax{t0+dt};

  // get bites from the previous time step [t0-dt,t0)
  while(!hpop->M2H_bites.empty()){

    // take out one S person (they move to E)
    hpop->S -= 1;
    hpop->E += 1;

    // time of the S->E event
    double btime = *hpop->M2H_bites.begin();

    // push the future infection (E->I) and clear the bite (S->E)
    hpop->H_inf.emplace(btime + hpop->LEP);
    hpop->M2H_bites.erase(hpop->M2H_bites.begin());

    // track state change
    hpop->hist.emplace_back(hist_elem{btime,-1,1,0});

  }

  assert(hpop->E == hpop->H_inf.size()-1);

  // push initial state into trace
  hpop->Sh_trace.emplace(queue_tuple(hpop->S,t0));
  hpop->Ih_trace.emplace(queue_tuple(hpop->I,t0));

  // recalculate propensity of I->R transition valid at t0
  hpop->ak = hpop->r * static_cast<double>(hpop->I);

  // simulate stochastic dynamics over [t0,t0+dt)
  while(hpop->tnow < tmax){

    // absolute next event times
    hpop->delta_t[0] = (hpop->Pk - hpop->Tk) / hpop->ak;
    hpop->delta_t[1] = *hpop->H_inf.begin() - hpop->tnow;

    // find minimum
    int mu{1};
    if(hpop->delta_t[0] < hpop->delta_t[1]){
      mu = 0;
    }

    // update putative time
    double delta = hpop->delta_t[mu];
    double tnext = hpop->tnow + delta;

    // if overrun, advance time to tmax and adjust the Poisson process for the fast forward
    if(tnext > tmax){
      double remaining = tmax - hpop->tnow;
      hpop->Tk += hpop->ak * remaining;
      hpop->tnow = tmax;
      break;
    }

    hpop->tnow = tnext;

    // I->S (recovery)
    if(mu == 0){
      hpop->I -= 1;
      hpop->S += 1;
      // track state change
      hpop->hist.emplace_back(hist_elem{hpop->tnow,1,0,-1});
    // E->I (completion of latent period)
    } else {
      hpop->E -= 1;
      hpop->I += 1;
      hpop->H_inf.erase(hpop->H_inf.begin());
      // track state change
      hpop->hist.emplace_back(hist_elem{hpop->tnow,0,-1,1});
    }

    // update Tk
    hpop->Tk += hpop->ak * delta;

    // update Pk
    if(mu == 0){
      hpop->Pk += log(1. / R::runif(0.,1.));
    }

    // recalculate propensity
    hpop->ak = hpop->r * static_cast<double>(hpop->I);

    // push to trace
    hpop->Sh_trace.emplace(queue_tuple(hpop->S,hpop->tnow));
    hpop->Ih_trace.emplace(queue_tuple(hpop->I,hpop->tnow));

  }

};
