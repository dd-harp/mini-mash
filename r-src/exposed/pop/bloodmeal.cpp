/* --------------------------------------------------------------------------------
#
#   disaggregated malaria transmission model: population level description
#   bloodmeal
#   Sean L. Wu (slwu89@berkeley.edu)
#   Dec 2020
#
-------------------------------------------------------------------------------- */

#include "bloodmeal.hpp"

// call this from the main loop
void bloodmeal(
  mosy_pop_ptr& mpop,
  human_pop_ptr& hpop,
  const double t0,
  const double dt,
  const Rcpp::NumericVector& parameters
){

  assert(mpop->H2M_bites.empty());
  assert(hpop->M2H_bites.empty());

  double tmax(t0+dt);

  // extract parameters
  double a = parameters["a"];
  double c = parameters["c"];
  double b = parameters["b"];

  double NH{hpop->N};

  // infections in mosquitoes
  bloodmeal_H2M(mpop,hpop,tmax,a,c,NH);

  // infections in humans
  bloodmeal_M2H(mpop,hpop,tmax,a,b,NH);

};

// sample S->E events in mosquitoes
void bloodmeal_H2M(
  mosy_pop_ptr& mpop,
  human_pop_ptr& hpop,
  const double tmax,
  const double a,
  const double c,
  const double NH
){

  // H2M transmission (S->E) events in mosquitoes
  double H2M_lambda{0.};
  int H2M_count{0};
  int H2M_tot{0};
  int mu;

  double t0,t1,delta; // bounds of current interval
  double SV,IH; // intensity is: (acX)SV, X = IH/NH

  std::array<double,2> tnext; // next interval

  // left limit of interval
  queue_tuple Sv_0 = *mpop->Sv_trace.begin();
  queue_tuple Ih_0 = *hpop->Ih_trace.begin();
  mpop->Sv_trace.erase(mpop->Sv_trace.begin());
  hpop->Ih_trace.erase(hpop->Ih_trace.begin());

  t0 = std::get<1>(Sv_0);

  // right limit of interval
  queue_tuple Sv_1 = *mpop->Sv_trace.begin();
  queue_tuple Ih_1 = *hpop->Ih_trace.begin();
  mpop->Sv_trace.erase(mpop->Sv_trace.begin());
  hpop->Ih_trace.erase(hpop->Ih_trace.begin());

  tnext[0] = std::get<1>(Sv_1);
  tnext[1] = std::get<1>(Ih_1);

  // compute H2M point process over [t0,t0+dt)
  while(true){

    // the last interval (to tmax)
    if( std::all_of(tnext.begin(),tnext.end(), [](const double y){return y == infinity;}) ){
      mu = -1; // sign to break
      t1 = tmax;
    } else {
      // find which trace changes next
      auto min_elem = std::min_element(tnext.begin(), tnext.end());
      mu = std::distance(tnext.begin(), min_elem);
      // compute length of this piecewise interval [t0,t1)
      t1 = tnext[mu];
    }

    delta = t1 - t0;

    // trajectory values at t0
    SV = std::get<0>(Sv_0);
    IH = std::get<0>(Ih_0);

    // intensity over [t0,t1)
    H2M_lambda = a * c * (IH/NH) * (SV - (double)H2M_tot) * delta;

    // sample bites
    H2M_count = R::rpois(H2M_lambda);
    H2M_tot += H2M_count;

    // push output
    for(int k=0; k<H2M_count; k++){
      double btime = R::runif(t0,t1);
      mpop->H2M_bites.emplace(btime);
    }

    // next change in Sv
    if(mu == 0){

      Sv_0 = Sv_1;

      if(!mpop->Sv_trace.empty()){

        Sv_1 = *mpop->Sv_trace.begin();
        mpop->Sv_trace.erase(mpop->Sv_trace.begin());
        tnext[0] = std::get<1>(Sv_1);

      } else {
        tnext[0] = infinity;
      }

    // next change in Ih
    } else if(mu == 1){

      Ih_0 = Ih_1;

      if(!hpop->Ih_trace.empty()){

        Ih_1 = *hpop->Ih_trace.begin();
        hpop->Ih_trace.erase(hpop->Ih_trace.begin());
        tnext[1] = std::get<1>(Ih_1);

      } else {
        tnext[1] = infinity;
      }

    // final interval
    } else {
      break;
    }

    // last t1 becomes new t0
    t0 = t1;

  }

};

// sample S->E events in humans
void bloodmeal_M2H(
  mosy_pop_ptr& mpop,
  human_pop_ptr& hpop,
  const double tmax,
  const double a,
  const double b,
  const double NH
){

  // M2H transmission (S->E) events in humans
  double M2H_lambda{0.};
  int M2H_count{0};
  int M2H_tot{0};
  int mu;

  double t0,t1,delta; // bounds of current interval
  double SH,IV; // intensity is: (abZ)IV, Z = SH/NH

  std::array<double,2> tnext; // next interval

  // left limit of interval
  queue_tuple Iv_0 = *mpop->Iv_trace.begin();
  queue_tuple Sh_0 = *hpop->Sh_trace.begin();
  mpop->Iv_trace.erase(mpop->Iv_trace.begin());
  hpop->Sh_trace.erase(hpop->Sh_trace.begin());

  t0 = std::get<1>(Iv_0);

  // right limit of interval
  queue_tuple Iv_1 = *mpop->Iv_trace.begin();
  queue_tuple Sh_1 = *hpop->Sh_trace.begin();
  mpop->Iv_trace.erase(mpop->Iv_trace.begin());
  hpop->Sh_trace.erase(hpop->Sh_trace.begin());

  tnext[0] = std::get<1>(Iv_1);
  tnext[1] = std::get<1>(Sh_1);

  // compute M2H point process over [t0,t0+dt)
  while(true){

    // the last interval (to tmax)
    if( std::all_of(tnext.begin(),tnext.end(), [](const double y){return y == infinity;}) ){
      mu = -1; // sign to break
      t1 = tmax;
    } else {
      // find which trace changes next
      auto min_elem = std::min_element(tnext.begin(), tnext.end());
      mu = std::distance(tnext.begin(), min_elem);
      // compute length of this piecewise interval [t0,t1)
      t1 = tnext[mu];
    }

    delta = t1 - t0;

    // trajectory values at t0
    IV = std::get<0>(Iv_0);
    SH = std::get<0>(Sh_0);

    // intensity over [t0,t1)
    M2H_lambda = a * b * ((SH - (double)M2H_tot)/NH) * IV * delta;

    // sample bites
    M2H_count = R::rpois(M2H_lambda);
    M2H_tot += M2H_count;

    // push output
    for(int k=0; k<M2H_count; k++){
      double btime = R::runif(t0,t1);
      hpop->M2H_bites.emplace(btime);
    }

    // next change in Iv
    if(mu == 0){

      Iv_0 = Iv_1;

      if(!mpop->Iv_trace.empty()){

        Iv_1 = *mpop->Iv_trace.begin();
        mpop->Iv_trace.erase(mpop->Iv_trace.begin());
        tnext[0] = std::get<1>(Iv_1);

      } else {
        tnext[0] = infinity;
      }

    // next change in Sh
    } else if(mu == 1){

      Sh_0 = Sh_1;

      if(!hpop->Sh_trace.empty()){

        Sh_1 = *hpop->Sh_trace.begin();
        hpop->Sh_trace.erase(hpop->Sh_trace.begin());
        tnext[1] = std::get<1>(Sh_1);

      } else {
        tnext[1] = infinity;
      }

    // final interval
    } else {
      break;
    }

    // last t1 becomes new t0
    t0 = t1;

  }

};
