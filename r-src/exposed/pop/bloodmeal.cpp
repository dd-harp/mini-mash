/* --------------------------------------------------------------------------------
#
#   disaggregated malaria transmission model: population level description
#   bloodmeal
#   Sean L. Wu (slwu89@berkeley.edu)
#   Dec 2020
#
-------------------------------------------------------------------------------- */

#include "bloodmeal.hpp"


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

  // H2M transmission (S->E) events in mosquitoes
  double H2M_lambda{0.};
  int H2M_count{0};
  int H2M_tot{0};

  double tl,tr; // bounds of current interval
  double SV,IH; // intensity is: (acX)SV

  std::array<double,2> tnext; // next interval

  // left limit of interval
  queue_tuple Sv_0 = *mpop->Sv_trace.begin();
  queue_tuple Ih_0 = *hpop->Ih_trace.begin();
  mpop->Sv_trace.erase(mpop->Sv_trace.begin());
  hpop->Ih_trace.erase(hpop->Ih_trace.begin());

  tl = std::get<1>(Sv_0);

  // right limit of interval
  queue_tuple Sv_1 = *mpop->Sv_trace.begin();
  queue_tuple Ih_1 = *hpop->Ih_trace.begin();
  mpop->Sv_trace.erase(mpop->Sv_trace.begin());
  hpop->Ih_trace.erase(hpop->Ih_trace.begin());

  tnext[0] = std::get<1>(Sv_1);
  tnext[1] = std::get<1>(Ih_1);

  // compute H2M point process over [t0,t0+dt)
  while( std::all_of(tnext.begin(),tnext.end(), [](const double y){return y < infinity;}) ){

    // find which trace changes next
    auto min_elem = std::min_element(tnext.begin(), tnext.end());
    int mu = std::distance(tnext.begin(), min_elem);

    // compute length of this piecewise interval [tl,tr)
    tr = tnext[mu];
    double delta = tr - tl;

    // trajectory values at tl
    SV = std::get<0>(Sv_0);
    IH = std::get<0>(Ih_0);

    // intensity over [tl,tr)
    H2M_lambda = a * c * (IH/NH) * (SV - (double)H2M_tot) * delta;

    // sample bites
    H2M_count = R::rpois(H2M_lambda);
    H2M_tot += H2M_count;

    for(int k=0; k<H2M_count; k++){
      double btime = R::runif(tl,tr);
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

    // next change in Iv
    } else {

      Ih_0 = Ih_1;

      if(!hpop->Ih_trace.empty()){

        Ih_1 = *hpop->Ih_trace.begin();
        hpop->Ih_trace.erase(hpop->Ih_trace.begin());
        tnext[1] = std::get<1>(Ih_1);

      } else {
        tnext[1] = infinity;
      }

    }

    // last tr becomes new tl
    tl = tr;

  }

};





void bloodmeal_old(
  mosy_pop_ptr& mpop,
  human_pop_ptr& hpop,
  const Rcpp::NumericVector& parameters
){

  assert(mpop->H2M_bites.empty());
  assert(hpop->M2H_bites.empty());

  double a = parameters["a"];
  double c = parameters["c"];
  double b = parameters["b"];

  double NH{hpop->N};

  // mosy S->E events occur with intensity (acX)S_{V}
  // human S->E events occur with intensity (abZ)I_{V}
  double mosyinf_intensity, humaninf_intensity;
  int mosyinf_samp, humaninf_samp;  // number of events that occur
  // cumulative number of events
  int mosyinf_tot{0};
  int humaninf_tot{0};

  // length of intervals
  double t1,t0,dt;
  // double S,I,X,Z;
  double SV,IV,SH,IH;

  // next state change
  std::array<double,4> t1_array{0.};

  // beginning of piecewise trajectories
  queue_tuple t0_SV = *mpop->Sv_trace.begin();
  queue_tuple t0_IV = *mpop->Iv_trace.begin();
  queue_tuple t0_SH = *hpop->Sh_trace.begin();
  queue_tuple t0_IH = *hpop->Ih_trace.begin();

  mpop->Sv_trace.erase(mpop->Sv_trace.begin());
  mpop->Iv_trace.erase(mpop->Iv_trace.begin());
  hpop->Sh_trace.erase(hpop->Sh_trace.begin());
  hpop->Ih_trace.erase(hpop->Ih_trace.begin());

  t0 = std::get<1>(t0_SV);

  // end of the piecewise trajectories
  queue_tuple t1_SV = *mpop->Sv_trace.begin();
  queue_tuple t1_IV = *mpop->Iv_trace.begin();
  queue_tuple t1_SH = *hpop->Sh_trace.begin();
  queue_tuple t1_IH = *hpop->Ih_trace.begin();

  mpop->Sv_trace.erase(mpop->Sv_trace.begin());
  mpop->Iv_trace.erase(mpop->Iv_trace.begin());
  hpop->Sh_trace.erase(hpop->Sh_trace.begin());
  hpop->Ih_trace.erase(hpop->Ih_trace.begin());

  t1_array[0] = std::get<1>(t1_SV);
  t1_array[1] = std::get<1>(t1_IV);
  t1_array[2] = std::get<1>(t1_SH);
  t1_array[3] = std::get<1>(t1_IH);

  // compute over the TWICE step
  while( std::all_of(t1_array.begin(), t1_array.end(), [](const double y){return y < infinity;}) ){

    mosyinf_intensity = 0.; // a c X S_v
    humaninf_intensity = 0.; // a b Z I_v

    // find which trajectory changes next
    auto min_elem = std::min_element(t1_array.begin(), t1_array.end());
    int mu = std::distance(t1_array.begin(), min_elem);

    // compute length of this piecewise interval [t0,d1)
    t1 = t1_array[mu];
    dt = t1 - t0;

    // state values at t0
    SV = std::get<0>(t0_SV);
    IV = std::get<0>(t0_IV);
    SH = std::get<0>(t0_SH);
    IH = std::get<0>(t0_IH);

    // intensity over the interval
    double X = IH / NH;
    double Z = (SH - static_cast<double>(humaninf_tot)) / NH;
    mosyinf_intensity = a * c * X * (SV - static_cast<double>(mosyinf_tot)) * dt;
    humaninf_intensity = a * b * Z * IV * dt;

    // // intensity over the interval
    // double X = IH / NH;
    // double Z = SH / NH;
    // mosyinf_intensity = a * c * X * SV * dt;
    // humaninf_intensity = a * b * Z * IV * dt;

    // sample values
    mosyinf_samp = R::rpois(mosyinf_intensity);
    humaninf_samp = R::rpois(humaninf_intensity);

    // cumulative S->E events
    mosyinf_tot += mosyinf_samp;
    humaninf_tot += humaninf_samp;

    // add the bites to the mosquito for next time
    for(int k=0; k<mosyinf_samp; k++){
      double btime = R::runif(t0,t1);
      mpop->H2M_bites.emplace(btime);
    }

    // add the bites to the human for next time
    for(int k=0; k<humaninf_samp; k++){
      double btime = R::runif(t0,t1);
      hpop->M2H_bites.emplace(btime);
    }

    // state change that happened first becomes new starting point
    if(mu==0){
      t0_SV = t1_SV;
      if(!mpop->Sv_trace.empty()){
        t1_SV = *mpop->Sv_trace.begin();
        mpop->Sv_trace.erase(mpop->Sv_trace.begin());
        t1_array[0] = std::get<1>(t1_SV);
      } else {
        t1_array[0] = infinity;
      }
    } else if(mu==1){
      t0_IV = t1_IV;
      if(!mpop->Iv_trace.empty()){
        t1_IV = *mpop->Iv_trace.begin();
        mpop->Iv_trace.erase(mpop->Iv_trace.begin());
        t1_array[1] = std::get<1>(t1_IV);
      } else {
        t1_array[1] = infinity;
      }
    } else if(mu==2){
      t0_SH = t1_SH;
      if(!hpop->Sh_trace.empty()){
        t1_SH = *hpop->Sh_trace.begin();
        hpop->Sh_trace.erase(hpop->Sh_trace.begin());
        t1_array[2] = std::get<1>(t1_SH);
      } else {
        t1_array[2] = infinity;
      }
    } else if(mu==3){
      t0_IH = t1_IH;
      if(!hpop->Ih_trace.empty()){
        t1_IH = *hpop->Ih_trace.begin();
        hpop->Ih_trace.erase(hpop->Ih_trace.begin());
        t1_array[3] = std::get<1>(t1_IH);
      } else {
        t1_array[3] = infinity;
      }
    } else {
      std::string msg("invalid minimum element in 'bloodmeal': " + std::to_string(mu));
      Rcpp::stop(msg);
    }

    // last end time becomes next beginning time
    t0 = t1;
  }


};
