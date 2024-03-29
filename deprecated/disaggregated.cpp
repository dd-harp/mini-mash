/* --------------------------------------------------------------------------------
#
#   disaggregated malaria transmission model
#   Sean L. Wu (slwu89@berkeley.edu)
#   November 2020
#
-------------------------------------------------------------------------------- */

// [[Rcpp::plugins(cpp14)]]

#include <vector>
#include <array>
#include <string>
#include <limits>
#include <tuple>
#include <set>
#include <memory>
#include <algorithm>
// #include <queue>

#include <Rcpp.h>


/* --------------------------------------------------------------------------------
#   Stuff the whole simulation needs
-------------------------------------------------------------------------------- */

// const static double eps = 1.E-9;
const static double infinity = std::numeric_limits<double>::infinity();

// data structure
using queue = std::set<double,std::less<double>>;

// data structure to pass traces
using queue_tuple = std::tuple<double,double>; // 1st element is state, 2nd is time

struct queue_comp {
    bool operator() (const queue_tuple& a, const queue_tuple& b) const {
      return std::get<1>(a) < std::get<1>(b);
    }
};

using queue_trace = std::set<queue_tuple,queue_comp>;


/* --------------------------------------------------------------------------------
#   syringes with wings: the mosquitoes
#   make the struct
-------------------------------------------------------------------------------- */

// mosquito population
typedef struct mosypop_str {

  // states
  int S;
  int E;
  int I;

  // history tracking
  // hist_set mosy_hist;
  std::vector<double> t_hist;
  std::vector<double> S_hist;
  std::vector<double> E_hist;
  std::vector<double> I_hist;

  double tnow;

  // bites when H->M transmission occured
  queue H2M_bites;

  // when the resulting infection in M will manifest
  queue M_inf;

  // parameters
  double g;
  double EIP;
  double lambda;

  // elements for MNRM of internal mosy dynamics
  std::array<double,5> delta_t{infinity};
  std::array<double,4> Pk{0.};
  std::array<double,4> Tk{0.};
  std::array<double,4> ak{0.};

  // integrated susceptible (Sv) biting over [t0,t0+dt)
  // integrated infectious (Iv) biting over [t0,t0+dt)
  queue_trace Sv_trace;
  queue_trace Iv_trace;

  mosypop_str(){Rcpp::Rcout << "'mosypop_str' ctor called at " << this << std::endl;};
  ~mosypop_str(){Rcpp::Rcout << "'mosypop_str' dtor called at " << this << std::endl;};

} mosypop_str;


// use this pointer
using mosypop_ptr = std::unique_ptr<mosypop_str>;

// make the mosquitoes
mosypop_ptr make_mosypop(const Rcpp::NumericVector parameters, const int SV, const int IV, const int outsize){

  mosypop_ptr mosypop = std::make_unique<mosypop_str>();

  mosypop->S = SV;
  mosypop->E = 0;
  mosypop->I = IV;

  mosypop->t_hist.reserve(outsize);
  mosypop->S_hist.reserve(outsize);
  mosypop->E_hist.reserve(outsize);
  mosypop->I_hist.reserve(outsize);

  mosypop->tnow = 0.0;

  mosypop->M_inf.emplace(infinity);

  mosypop->g = parameters["g"];
  mosypop->lambda = parameters["lambda"];
  mosypop->EIP = parameters["EIP"];

  // calculate initial propensities
  mosypop->ak[0] = mosypop->lambda;
  mosypop->ak[1] = mosypop->g * static_cast<double>(mosypop->S);
  mosypop->ak[2] = mosypop->g * static_cast<double>(mosypop->E);
  mosypop->ak[3] = mosypop->g * static_cast<double>(mosypop->I);

  // draw initial firing times
  for(int j=0; j<4; j++){
    mosypop->Pk[j] = log(1. / R::runif(0.,1.));
  }

  return mosypop;
};


/* --------------------------------------------------------------------------------
#   Simulate the mosquitoes over a TWICE step
-------------------------------------------------------------------------------- */

// iterate the mosquitoes
void run_mosypop(mosypop_ptr& mosypop, double t0, double dt){

  double tmax{t0+dt};

  // traces will give relevant state changes to bloodmeal, over [t0,tmax)
  if(!mosypop->Sv_trace.empty() || !mosypop->Iv_trace.empty()){
    Rcpp::stop("'Sv_trace' and 'Iv_trace' should always be empty at start of time step \n");
  }

  // get bites from the previous time step [t0-dt,t0)
  while(!mosypop->H2M_bites.empty()){

    // take out one S mosquito and move to E
    // mosypop->S -= 1;

    // time of that bite (S->E)
    double btime = *mosypop->H2M_bites.begin();
    assert(btime < t0);

    // find out how much hazard they accumulate between [btime,t0)
    double haz = (t0 - btime) * mosypop->g;
    assert(haz > 0);

    // did they survive until t0? if not, they died between [btime,t0)
    if(R::runif(0.,1.) < exp(-haz)){
      // queue the future EIP completion event, which will happen on this time step (or next)
      mosypop->M_inf.emplace(btime + mosypop->EIP);
      mosypop->S -= 1;
      mosypop->E += 1;
    }

    // pop that bite off
    mosypop->H2M_bites.erase(mosypop->H2M_bites.begin());

  }
  // end looping over bites from last step

  mosypop->t_hist.push_back(mosypop->tnow);
  mosypop->S_hist.push_back(mosypop->S);
  mosypop->E_hist.push_back(mosypop->E);
  mosypop->I_hist.push_back(mosypop->I);


  assert(mosypop->E == mosypop->M_inf.size() - 1); // -1 because always has Inf

  // recalculate propentities valid at t0
  mosypop->ak[1] = mosypop->g * static_cast<double>(mosypop->S);
  mosypop->ak[2] = mosypop->g * static_cast<double>(mosypop->E);
  mosypop->ak[3] = mosypop->g * static_cast<double>(mosypop->I);

  // push initial state into traces
  mosypop->Sv_trace.emplace(queue_tuple(mosypop->S,t0));
  mosypop->Iv_trace.emplace(queue_tuple(mosypop->I,t0));

  // simulate dynamics over this time step [t0,tmax)
  while(mosypop->tnow < tmax){

    // absolute next firing times
    for(int j=0; j<4; j++){
      mosypop->delta_t[j] = (mosypop->Pk[j] - mosypop->Tk[j]) / mosypop->ak[j];
    }
    mosypop->delta_t[4] = *mosypop->M_inf.begin() - mosypop->tnow;

    // find minimum
    auto min_elem = std::min_element(mosypop->delta_t.begin(), mosypop->delta_t.end());
    int mu = std::distance(mosypop->delta_t.begin(), min_elem);

    // update putative time
    double delta = mosypop->delta_t[mu];
    double tsamp = mosypop->tnow + delta;

    // if tsamp > tmax, advance time to tmax and adjust the Poisson processes for the fast-forward
    if(tsamp > tmax){
      double remaining = tmax - mosypop->tnow;
      for(int j=0; j<4; j++){
        mosypop->Tk[j] += mosypop->ak[j] * remaining;
      }
      mosypop->Sv_trace.emplace(queue_tuple(mosypop->S,tmax));
      mosypop->Iv_trace.emplace(queue_tuple(mosypop->I,tmax));
      mosypop->tnow = tmax;
      break;
    }
    mosypop->tnow = tsamp;

    // update system

    // S emerges
    if(mu == 0){
      mosypop->S += 1;
      mosypop->Sv_trace.emplace(queue_tuple(mosypop->S,mosypop->tnow));
    // S dies
    } else if(mu == 1){
      mosypop->S -= 1;
      mosypop->Sv_trace.emplace(queue_tuple(mosypop->S,mosypop->tnow));
    // E dies
    } else if(mu == 2){
      mosypop->E -= 1;
      // delete a random element from the EIP set (not the end, which is Inf)
      int relem = R::runif(0.0, mosypop->M_inf.size()-1);
      auto it = mosypop->M_inf.begin();
      std::advance(it, relem);
      mosypop->M_inf.erase(it);
    // I dies
    } else if(mu == 3){
      mosypop->I -= 1;
      mosypop->Iv_trace.emplace(queue_tuple(mosypop->I,mosypop->tnow));
    // E->I completion
    } else if(mu == 4){
      mosypop->E -= 1;
      mosypop->I += 1;
      mosypop->M_inf.erase(mosypop->M_inf.begin());
      mosypop->Iv_trace.emplace(queue_tuple(mosypop->I,mosypop->tnow));
    } else {
      std::string msg("invalid minimum element in 'run_mosypop': " + std::to_string(mu));
      Rcpp::stop(msg);
    }

    // update Tk
    for(int j=0; j<4; j++){
      mosypop->Tk[j] += mosypop->ak[j] * delta;
    }

    // update Pk[mu]
    if(mu != 4){
      mosypop->Pk[mu] += log(1. / R::runif(0.,1.));
    }

    // recalculate propentities
    mosypop->ak[1] = mosypop->g * static_cast<double>(mosypop->S);
    mosypop->ak[2] = mosypop->g * static_cast<double>(mosypop->E);
    mosypop->ak[3] = mosypop->g * static_cast<double>(mosypop->I);

  }

}


/* --------------------------------------------------------------------------------
#   human population
-------------------------------------------------------------------------------- */

typedef struct humanpop_str {

  // states
  int S;
  int E;
  int I;

  // history tracking
  // hist_set human_hist;
  std::vector<double> t_hist;
  std::vector<double> S_hist;
  std::vector<double> E_hist;
  std::vector<double> I_hist;

  double tnow;

  // bites when M->H transmission occured
  queue M2H_bites;

  // when the resulting infection in M will manifest
  queue H_inf;

  // parameters
  double r;
  double LEP;

  // elements for MNRM of internal dynamics
  std::array<double,2> delta_t{infinity};
  double Pk{0.};
  double Tk{0.};
  double ak{0.};

  // state transitions for bloodmeal over [t0,t0+dt)
  queue_trace Sh_trace;
  queue_trace Ih_trace;

  humanpop_str(){Rcpp::Rcout << "'humanpop_str' ctor called at " << this << std::endl;};
  ~humanpop_str(){Rcpp::Rcout << "'humanpop_str' dtor called at " << this << std::endl;};

} humanpop_str;


// use this pointer
using humanpop_ptr = std::unique_ptr<humanpop_str>;

// make the humans
humanpop_ptr make_humanpop(const Rcpp::NumericVector parameters, const int SH, const int IH, const int outsize){

  humanpop_ptr humanpop = std::make_unique<humanpop_str>();

  humanpop->S = SH;
  humanpop->E = 0;
  humanpop->I = IH;

  humanpop->t_hist.reserve(outsize);
  humanpop->S_hist.reserve(outsize);
  humanpop->E_hist.reserve(outsize);
  humanpop->I_hist.reserve(outsize);

  humanpop->tnow = 0.0;

  humanpop->H_inf.emplace(infinity);

  humanpop->r = parameters["r"];
  humanpop->LEP = parameters["LEP"];

  // calculate initial propensities
  humanpop->ak = humanpop->r * static_cast<double>(humanpop->I);

  // draw initial firing times
  humanpop->Pk = log(1. / R::runif(0.,1.));

  return humanpop;
};


/* --------------------------------------------------------------------------------
#   simulate the humans over a TWICE step
-------------------------------------------------------------------------------- */

// iterate the humans
void run_humanpop(humanpop_ptr& humanpop, double t0, double dt){

  // push initial state into traces
  if(!humanpop->Sh_trace.empty() || !humanpop->Ih_trace.empty()){
    Rcpp::stop("'S_trace' and 'I_trace' should always be empty at start of time step \n");
  }

  double tmax{t0+dt};

  // get bites from the previous time step [t0-dt,t0)
  while(!humanpop->M2H_bites.empty()){

    // take out one S person (they move to E)
    humanpop->S -= 1;
    humanpop->E += 1;

    // time of the S->E event
    double btime = *humanpop->M2H_bites.begin();

    // push the future infection
    humanpop->H_inf.emplace(btime + humanpop->LEP);

    humanpop->M2H_bites.erase(humanpop->M2H_bites.begin());

  }
  // end loop over previous bites

  humanpop->t_hist.push_back(humanpop->tnow);
  humanpop->S_hist.push_back(humanpop->S);
  humanpop->E_hist.push_back(humanpop->E);
  humanpop->I_hist.push_back(humanpop->I);

  assert(humanpop->E == humanpop->H_inf.size() - 1);

  // recalculate propensity valid at t0
  humanpop->ak = humanpop->r * static_cast<double>(humanpop->I);

  // push initial state into trace
  humanpop->Sh_trace.emplace(queue_tuple(humanpop->S,t0));
  humanpop->Ih_trace.emplace(queue_tuple(humanpop->I,t0));

  // simulate dynamics over this time step [t0,tmax)
  while(humanpop->tnow < tmax){

    // absolute next firing times
    humanpop->delta_t[0] = (humanpop->Pk - humanpop->Tk) / humanpop->ak;
    humanpop->delta_t[1] = *humanpop->H_inf.begin() - humanpop->tnow;

    // find minimum
    int mu;
    if(humanpop->delta_t[0] < humanpop->delta_t[1]){
      mu = 0;
    } else {
      mu = 1;
    }

    // update putative time
    double delta = humanpop->delta_t[mu];
    double tsamp = humanpop->tnow + delta;

    // if overrun, advance time to tmax and adjust the Poisson process for the fast forward
    if(tsamp > tmax){
      double remaining = tmax - humanpop->tnow;
      humanpop->Tk += humanpop->ak * remaining;
      humanpop->Sh_trace.emplace(queue_tuple(humanpop->S,tmax));
      humanpop->Ih_trace.emplace(queue_tuple(humanpop->I,tmax));
      humanpop->tnow = tmax;
      break;
    }
    humanpop->tnow = tsamp;

    // update system
    // I->S recovery
    if(mu == 0){
      humanpop->I -= 1;
      humanpop->S += 1;
    // E->I
    } else if(mu == 1){
      humanpop->E -= 1;
      humanpop->I += 1;
      humanpop->H_inf.erase(humanpop->H_inf.begin());
    } else {
      std::string msg("invalid minimum element in 'run_humanpop': " + std::to_string(mu));
      Rcpp::stop(msg);
    }

    humanpop->Sh_trace.emplace(queue_tuple(humanpop->S,humanpop->tnow));
    humanpop->Ih_trace.emplace(queue_tuple(humanpop->I,humanpop->tnow));

    // update Tk
    humanpop->Tk += humanpop->ak * delta;

    // update Pk
    if(mu == 0){
      humanpop->Pk += log(1. / R::runif(0.,1.));
    }

    // recalculate propensity
    humanpop->ak = humanpop->r * static_cast<double>(humanpop->I);

  }

}


/* --------------------------------------------------------------------------------
#   bloodmeal
-------------------------------------------------------------------------------- */

void bloodmeal(const Rcpp::NumericVector parameters, humanpop_ptr& humanpop, mosypop_ptr& mosypop){

  assert(mosypop->H2M_bites.size() == 0);
  assert(humanpop->M2H_bites.size() == 0);

  double a = parameters["a"];
  double c = parameters["c"];
  double b = parameters["b"];

  // human population size is constant
  double NH = static_cast<double>(humanpop->S) + static_cast<double>(humanpop->E) + static_cast<double>(humanpop->I);

  // mosy S->E events occur with intensity (acX)S_{V}
  // human S->E events occur with intensity (abZ)I_{V}
  double mosyinf_intensity, humaninf_intensity;
  int mosyinf_samp, humaninf_samp;  // number of events that occur
  // cumulative number of events
  int mosyinf_cum{0};
  int humaninf_cum{0};

  // length of intervals
  double t1,t0,dt;
  // double S,I,X,Z;
  double SV,IV,SH,IH;

  // next state change
  std::array<double,4> t1_array{0.};

  // beginning of piecewise trajectories
  // Rcpp::Rcout << "accessing t0\n";
  queue_tuple t0_SV = *mosypop->Sv_trace.begin();
  queue_tuple t0_IV = *mosypop->Iv_trace.begin();
  queue_tuple t0_SH = *humanpop->Sh_trace.begin();
  queue_tuple t0_IH = *humanpop->Ih_trace.begin();

  mosypop->Sv_trace.erase(mosypop->Sv_trace.begin());
  mosypop->Iv_trace.erase(mosypop->Iv_trace.begin());
  humanpop->Sh_trace.erase(humanpop->Sh_trace.begin());
  humanpop->Ih_trace.erase(humanpop->Ih_trace.begin());

  t0 = std::get<1>(t0_SV);

  // end of the piecewise trajectories
  // Rcpp::Rcout << "accessing t1\n";
  queue_tuple t1_SV = *mosypop->Sv_trace.begin();
  queue_tuple t1_IV = *mosypop->Iv_trace.begin();
  queue_tuple t1_SH = *humanpop->Sh_trace.begin();
  queue_tuple t1_IH = *humanpop->Ih_trace.begin();

  mosypop->Sv_trace.erase(mosypop->Sv_trace.begin());
  mosypop->Iv_trace.erase(mosypop->Iv_trace.begin());
  humanpop->Sh_trace.erase(humanpop->Sh_trace.begin());
  humanpop->Ih_trace.erase(humanpop->Ih_trace.begin());

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
    double Z = (SH - static_cast<double>(humaninf_cum)) / NH;
    mosyinf_intensity = a * c * X * (SV - static_cast<double>(mosyinf_cum)) * dt;
    humaninf_intensity = a * b * Z * IV * dt;

    // // intensity over the interval
    // double X = IH / NH;
    // double Z = SH / NH;
    // mosyinf_intensity = a * c * X * SV * dt;
    // humaninf_intensity = a * b * Z * IV * dt;

    // sample values
    mosyinf_samp = R::rpois(mosyinf_intensity);
    humaninf_samp = R::rpois(humaninf_intensity);

    // // cumulative S->E events
    mosyinf_cum += mosyinf_samp;
    humaninf_cum += humaninf_samp;

    // add the bites to the mosquito for next time
    if(mosyinf_samp > 0){
      for(int k=0; k<mosyinf_samp; k++){
        double btime = R::runif(t0,t1);
        mosypop->H2M_bites.emplace(btime);
      }
    }

    // add the bites to the human for next time
    if(humaninf_samp > 0){
      for(int k=0; k<humaninf_samp; k++){
        double btime = R::runif(t0,t1);
        humanpop->M2H_bites.emplace(btime);
      }
    }

    // state change that happened first becomes new starting point
    if(mu==0){
      t0_SV = t1_SV;
      if(!mosypop->Sv_trace.empty()){
        t1_SV = *mosypop->Sv_trace.begin();
        mosypop->Sv_trace.erase(mosypop->Sv_trace.begin());
        t1_array[0] = std::get<1>(t1_SV);
      } else {
        t1_array[0] = infinity;
      }
    } else if(mu==1){
      t0_IV = t1_IV;
      if(!mosypop->Iv_trace.empty()){
        t1_IV = *mosypop->Iv_trace.begin();
        mosypop->Iv_trace.erase(mosypop->Iv_trace.begin());
        t1_array[1] = std::get<1>(t1_IV);
      } else {
        t1_array[1] = infinity;
      }
    } else if(mu==2){
      t0_SH = t1_SH;
      if(!humanpop->Sh_trace.empty()){
        t1_SH = *humanpop->Sh_trace.begin();
        humanpop->Sh_trace.erase(humanpop->Sh_trace.begin());
        t1_array[2] = std::get<1>(t1_SH);
      } else {
        t1_array[2] = infinity;
      }
    } else if(mu==3){
      t0_IH = t1_IH;
      if(!humanpop->Ih_trace.empty()){
        t1_IH = *humanpop->Ih_trace.begin();
        humanpop->Ih_trace.erase(humanpop->Ih_trace.begin());
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

  assert(mosyinf_cum == mosypop->H2M_bites.size());
  assert(humaninf_cum == humanpop->M2H_bites.size());

};


/* --------------------------------------------------------------------------------
#   run (main)
-------------------------------------------------------------------------------- */

// [[Rcpp::export]]
Rcpp::List run_miniMASH_exactbm(const Rcpp::NumericVector parameters, const Rcpp::IntegerVector y0, const double dt, const double tmax, bool verbose){

  int SH = y0["SH"];
  int IH = y0["IH"];
  int SV = y0["SV"];
  int IV = y0["IV"];

  const int outsize = 1E5;
  mosypop_ptr mosypop = make_mosypop(parameters,SV,IV,outsize);
  humanpop_ptr humanpop = make_humanpop(parameters,SH,IH,outsize);

  int tsteps = tmax/dt;
  int i{0};
  double clock{0.0};

  Rcpp::Rcout << "--- begin simulation ---\n";
  while(clock < tmax){

    run_mosypop(mosypop, clock, dt);
    run_humanpop(humanpop, clock, dt);

    bloodmeal(parameters,humanpop,mosypop);

    i++;
    if((i % 100 == 0) && verbose){
      Rcpp::Rcout << " completed TWICE step " << i << " of " << tsteps << "\n";
    }
    clock += dt;
  }
  Rcpp::Rcout << "--- ending simulation ---\n";


  int kh = humanpop->t_hist.size();
  Rcpp::NumericMatrix outh(kh,4);
  Rcpp::colnames(outh) = Rcpp::CharacterVector::create("time","SH","EH","IH");
  for(int j=0; j<kh; j++){
    outh(j,0) = humanpop->t_hist[j];
    outh(j,1) = humanpop->S_hist[j];
    outh(j,2) = humanpop->E_hist[j];
    outh(j,3) = humanpop->I_hist[j];
  }

  int km = mosypop->t_hist.size();
  Rcpp::NumericMatrix outm(km,4);
  Rcpp::colnames(outm) = Rcpp::CharacterVector::create("time","SV","EV","IV");
  for(int j=0; j<km; j++){
    outm(j,0) = mosypop->t_hist[j];
    outm(j,1) = mosypop->S_hist[j];
    outm(j,2) = mosypop->E_hist[j];
    outm(j,3) = mosypop->I_hist[j];
  }

  return Rcpp::List::create(
    Rcpp::Named("mosquito") = outm,
    Rcpp::Named("human") = outh
  );
}
