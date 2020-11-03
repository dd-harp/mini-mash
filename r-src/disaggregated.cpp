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
#include <memory>
#include <queue>

#include <Rcpp.h>


/* --------------------------------------------------------------------------------
#   Stuff the whole simulation needs
-------------------------------------------------------------------------------- */

// const static double eps = 1.E-9;
const static double infinity = std::numeric_limits<double>::infinity();

// data structure
using queue = std::priority_queue<double,std::vector<double>,std::greater<double> >;

// data structure to pass traces
using queue_tuple = std::tuple<double,double>; // 1st element is state, 2nd is time
auto queue_comp = [](const queue_tuple& a, const queue_tuple& b ) -> bool { return std::get<1>(a) > std::get<1>(b); };
using queue_trace = std::priority_queue<queue_tuple,std::vector<queue_tuple>, decltype(queue_comp) >;

// [[Rcpp::export]]
void testqueue(){
  queue_trace sample(queue_comp);
  sample.push(queue_tuple(1,R::rlnorm(0.,1.)));
  sample.push(queue_tuple(2,R::rlnorm(0.,1.)));
  sample.push(queue_tuple(3,R::rlnorm(0.,1.)));
  sample.push(queue_tuple(4,infinity));
  while(!sample.empty()){
    queue_tuple n = sample.top();
    sample.pop();
    Rcpp::Rcout <<  " int val: " << std::get<0>(n) << " double val: " << std::get<1>(n) << " --- \n";
  }
};



/* --------------------------------------------------------------------------------
#   syringes with wings: the mosquitoes
-------------------------------------------------------------------------------- */

// mosquito population
typedef struct mosypop_str {

  // states
  int S;
  int I;

  // history tracking
  std::vector<double> t_hist;
  std::vector<double> S_hist;
  std::vector<double> I_hist;

  double tnow;

  // bites when H->M transmission occured
  queue H2M_bites;

  // when the resulting infection in M will manifest
  queue M_inf;

  // parameters
  double a;
  double g;
  double EIP;
  double lambda;

  // elements for MNRM of internal mosy dynamics
  std::array<double,4> delta_t{infinity};
  std::array<double,3> Pk{0.};
  std::array<double,3> Tk{0.};
  std::array<double,3> ak{0.};

  // integrated susceptible (Sv) biting over [t0,t0+dt)
  // integrated infectious (Iv) biting over [t0,t0+dt)
  queue_trace Sv_trace;
  queue_trace Iv_trace;

  mosypop_str() : Sv_trace(queue_comp), Iv_trace(queue_comp) {};

} mosypop_str;

// use this pointer
using mosypop_ptr = std::unique_ptr<mosypop_str>;

// make the mosquitoes
mosypop_ptr make_mosypop(const Rcpp::NumericVector parameters, const int SV, const int IV, const int outsize){
  mosypop_ptr mosypop = std::make_unique<mosypop_str>();

  mosypop->S = SV;
  mosypop->I = IV;

  mosypop->t_hist.reserve(outsize);
  mosypop->S_hist.reserve(outsize);
  mosypop->I_hist.reserve(outsize);

  mosypop->tnow = 0.0;

  // mosypop->H2M_bites.push(infinity);
  mosypop->M_inf.push(infinity);

  mosypop->a = parameters["a"];
  mosypop->g = parameters["g"];
  mosypop->lambda = parameters["lambda"];
  mosypop->EIP = parameters["EIP"];

  // calculate initial propensities
  mosypop->ak[0] = mosypop->lambda;
  mosypop->ak[1] = mosypop->g * static_cast<double>(mosypop->S);
  mosypop->ak[2] = mosypop->g * static_cast<double>(mosypop->I);

  // draw initial firing times
  for(int j=0; j<3; j++){
    mosypop->Pk[j] = log(1. / R::runif(0.,1.));
  }

  return mosypop;
};

// iterate the mosquitoes
void run_mosypop(mosypop_ptr& mosypop, double t0, double dt){

  // push initial state into traces
  if(!mosypop->Sv_trace.empty() || !mosypop->Iv_trace.empty()){
    Rcpp::stop("'Sv_trace' and 'Iv_trace' should always be empty at start of time step \n");
  }
  mosypop->Sv_trace.push(queue_tuple(mosypop->S,t0));
  mosypop->Iv_trace.push(queue_tuple(mosypop->I,t0));

  double tmax{t0+dt};
  double P = exp(-mosypop->g*mosypop->EIP);

  // Rcpp::Rcout << "getting bites from prev time step \n";

  // get bites from the previous time step
  while(!mosypop->H2M_bites.empty()){
    // take out one S mosquito (they move to E)
    mosypop->S -= 1;
    // if that mosquito would survive the EIP, queue an infection event in the future
    if(R::runif(0.,1.) < P){
      mosypop->M_inf.push(mosypop->H2M_bites.top() + mosypop->EIP);
    }
    mosypop->H2M_bites.pop();
  }

  // simulate dynamics over this time step
  while(mosypop->tnow < tmax){

    // absolute next firing times
    for(int j=0; j<3; j++){
      mosypop->delta_t[j] = (mosypop->Pk[j] - mosypop->Tk[j]) / mosypop->ak[j];
    }
    mosypop->delta_t[3] = mosypop->M_inf.top() - mosypop->tnow;

    // find minimum
    auto min_elem = std::min_element(mosypop->delta_t.begin(), mosypop->delta_t.end());
    int mu = std::distance(mosypop->delta_t.begin(), min_elem);

    // update putative time
    double delta = mosypop->delta_t[mu];
    double tsamp = mosypop->tnow + delta;
    // if overrun, advance time to tmax and adjust the Poisson processes for the fast-forward
    if(tsamp > tmax){
      double remaining = tmax - mosypop->tnow;
      for(int j=0; j<3; j++){
        mosypop->Tk[j] += mosypop->ak[j] * remaining;
      }
      mosypop->Sv_trace.push(queue_tuple(mosypop->S,tmax));
      mosypop->Iv_trace.push(queue_tuple(mosypop->I,tmax));
      mosypop->tnow = tmax;
      break;
    }
    mosypop->tnow = tsamp;

    // update system
    if(mu == 0){
      mosypop->S += 1;
      mosypop->Sv_trace.push(queue_tuple(mosypop->S,mosypop->tnow));
    } else if(mu == 1){
      mosypop->S -= 1;
      mosypop->Sv_trace.push(queue_tuple(mosypop->S,mosypop->tnow));
    } else if(mu == 2){
      mosypop->I -= 1;
      mosypop->Iv_trace.push(queue_tuple(mosypop->I,mosypop->tnow));
    } else if(mu == 3){
      mosypop->M_inf.pop();
      mosypop->I += 1;
      mosypop->Iv_trace.push(queue_tuple(mosypop->I,mosypop->tnow));
    } else {
      std::string msg("invalid minimum element in 'run_mosypop': " + std::to_string(mu));
      Rcpp::stop(msg);
    }

    // update Tk
    for(int j=0; j<3; j++){
      mosypop->Tk[j] += mosypop->ak[j] * delta;
    }

    // update Pk[mu]
    if(mu != 3){
      mosypop->Pk[mu] += log(1. / R::runif(0.,1.));
    }

    // recalculate propentities
    mosypop->ak[1] = mosypop->g * static_cast<double>(mosypop->S);
    mosypop->ak[2] = mosypop->g * static_cast<double>(mosypop->I);

    // track output
    mosypop->t_hist.push_back(mosypop->tnow);
    mosypop->S_hist.push_back(mosypop->S);
    mosypop->I_hist.push_back(mosypop->I);

  }

}


/* --------------------------------------------------------------------------------
#   human population
-------------------------------------------------------------------------------- */

typedef struct humanpop_str {

  // states
  int S;
  int I;
  // double N;

  // history tracking
  std::vector<double> t_hist;
  std::vector<double> S_hist;
  std::vector<double> I_hist;

  double tnow;

  // bites when M->H transmission occured
  queue M2H_bites;

  // when the resulting infection in M will manifest
  queue H_inf;

  // parameters
  double r;
  double LEP;

  // elements for MNRM of internal mosy dynamics
  std::array<double,2> delta_t{infinity};
  double Pk{0.};
  double Tk{0.};
  double ak{0.};

  // integrated prevalence over [t0,t0+dt)
  queue_trace X_trace;

  humanpop_str() : X_trace(queue_comp) {};

} humanpop_str;

// use this pointer
using humanpop_ptr = std::unique_ptr<humanpop_str>;

// make the humans
humanpop_ptr make_humanpop(const Rcpp::NumericVector parameters, const int SH, const int IH, const int outsize){
  humanpop_ptr humanpop = std::make_unique<humanpop_str>();

  humanpop->S = SH;
  humanpop->I = IH;
  // humanpop->N = SH+IH;

  humanpop->t_hist.reserve(outsize);
  humanpop->S_hist.reserve(outsize);
  humanpop->I_hist.reserve(outsize);

  humanpop->tnow = 0.0;

  humanpop->H_inf.push(infinity);

  humanpop->r = parameters["r"];
  humanpop->LEP = parameters["LEP"];

  // calculate initial propensities
  humanpop->ak = humanpop->r * static_cast<double>(humanpop->I);

  // draw initial firing times
  humanpop->Pk = log(1. / R::runif(0.,1.));

  return humanpop;
};

// iterate the humans
void run_humanpop(humanpop_ptr& humanpop, double t0, double dt){

  // push initial state into traces
  if(!humanpop->X_trace.empty()){
    Rcpp::stop("'X_trace' should always be empty at start of time step \n");
  }

  double tmax{t0+dt};

  double X = static_cast<double>(humanpop->I) / (static_cast<double>(humanpop->S) + static_cast<double>(humanpop->I));
  humanpop->X_trace.push(queue_tuple(X,t0));

  // get bites from the previous time step
  while(!humanpop->M2H_bites.empty()){
    // take out one S person (they move to E)
    humanpop->S -= 1;
    // push the future infection
    humanpop->H_inf.push(humanpop->M2H_bites.top() + humanpop->LEP);
    humanpop->M2H_bites.pop();
  }

  // simulate dynamics over this time step
  while(humanpop->tnow < tmax){

    // absolute next firing times
    humanpop->delta_t[0] = (humanpop->Pk - humanpop->Tk) / humanpop->ak;
    humanpop->delta_t[1] = humanpop->H_inf.top() - humanpop->tnow;

    // find minimum
    int mu{1};
    if(humanpop->delta_t[0] < humanpop->delta_t[1]){
      mu = 0;
    }

    // update putative time
    double delta = humanpop->delta_t[mu];
    double tsamp = humanpop->tnow + delta;
    // if overrun, advance time to tmax and adjust the Poisson process for the fast forward
    if(tsamp > tmax){
      double remaining = tmax - humanpop->tnow;
      humanpop->Tk += humanpop->ak * remaining;
      X = static_cast<double>(humanpop->I) / (static_cast<double>(humanpop->S) + static_cast<double>(humanpop->I));
      humanpop->X_trace.push(queue_tuple(X,tmax));
      humanpop->tnow = tmax;
      break;
    }
    humanpop->tnow = tsamp;

    // update system
    if(mu == 0){
      humanpop->I -= 1;
      humanpop->S += 1;
    } else if(mu == 1){
      humanpop->H_inf.pop();
      humanpop->I += 1;
    } else {
      std::string msg("invalid minimum element in 'run_mosypop': " + std::to_string(mu));
      Rcpp::stop(msg);
    }

    X = static_cast<double>(humanpop->I) / (static_cast<double>(humanpop->S) + static_cast<double>(humanpop->I));
    humanpop->X_trace.push(queue_tuple(X,humanpop->tnow));

    // update Tk
    humanpop->Tk += humanpop->ak * delta;

    // update Pk
    if(mu == 0){
      humanpop->Pk += log(1. / R::runif(0.,1.));
    }

    // recalculate propensity
    humanpop->ak = humanpop->r * static_cast<double>(humanpop->I);

    // track output
    humanpop->t_hist.push_back(humanpop->tnow);
    humanpop->S_hist.push_back(humanpop->S);
    humanpop->I_hist.push_back(humanpop->I);

  }

}


/* --------------------------------------------------------------------------------
#   bloodmeal
-------------------------------------------------------------------------------- */

void bloodmeal(const Rcpp::NumericVector parameters, humanpop_ptr& humanpop, mosypop_ptr& mosypop){

  double a = parameters["a"];
  double c = parameters["c"];
  double b = parameters["b"];

  // H2M bites occur with intensity (acX)S_{V}
  // M2H bites occur with intensity ab(1-X)I_{V}
  double H2M_intensity, M2H_intensity;
  int H2M_bites, M2H_bites;

  // length of intervals
  double t1,t0,dt;
  double S,I,X;

  // next state change
  std::array<double,3> t1_array{0.};

  // beginning of piecewise trajectories
  queue_tuple t0_SV = mosypop->Sv_trace.top();
  queue_tuple t0_IV = mosypop->Iv_trace.top();
  queue_tuple t0_X = humanpop->X_trace.top();

  mosypop->Sv_trace.pop();
  mosypop->Iv_trace.pop();
  humanpop->X_trace.pop();

  t0 = std::get<1>(t0_SV);

  // end of the piecewise trajectories
  queue_tuple t1_SV;
  queue_tuple t1_IV;
  queue_tuple t1_X;

  // compute over the TWICE step
  while(!mosypop->Sv_trace.empty() || !mosypop->Iv_trace.empty() || !humanpop->X_trace.empty()){

    H2M_intensity = 0.;
    M2H_intensity = 0.;

    // if not empty, get the next state change in each trajectory
    std::fill(t1_array.begin(),t1_array.end(),infinity);
    if(!mosypop->Sv_trace.empty()){
      t1_SV = mosypop->Sv_trace.top();
      mosypop->Sv_trace.pop();
      t1_array[0] = std::get<1>(t1_SV);
    }
    if(!mosypop->Iv_trace.empty()){
      t1_IV = mosypop->Iv_trace.top();
      mosypop->Iv_trace.pop();
      t1_array[1] = std::get<1>(t1_IV);
    }
    if(!humanpop->X_trace.empty()){
      t1_X = humanpop->X_trace.top();
      humanpop->X_trace.pop();
      t1_array[2] = std::get<1>(t1_X);
    }

    // find which trajectory changes next
    auto min_elem = std::min_element(t1_array.begin(), t1_array.end());
    int mu = std::distance(t1_array.begin(), min_elem);

    // compute length of this piecewise interval [t0,d1)
    t1 = t1_array[mu];
    dt = t1 - t0;

    // state values at t0
    S = std::get<0>(t0_SV);
    I = std::get<0>(t0_IV);
    X = std::get<0>(t0_X);

    // intensity over the interval
    H2M_intensity = a * c * X * S * dt;
    M2H_intensity = a * b * (1. - X) * I * dt;

    // sample values
    H2M_bites = R::rpois(H2M_intensity);
    M2H_bites = R::rpois(M2H_intensity);

    // add the bites to the mosquito for next time
    if(H2M_bites > 0){
      double btime = R::runif(t0,t1);
      mosypop->H2M_bites.push(btime);
    }

    // add the bites to the human for next time
    if(M2H_bites > 0){
      double btime = R::runif(t0,t1);
      humanpop->M2H_bites.push(btime);
    }

    // state change that happened first becomes new starting point
    if(mu==0){
      t0_SV = t1_SV;
    } else if(mu==1){
      t0_IV = t1_IV;
    } else if(mu==2){
      t0_X = t1_X;
    } else{
      std::string msg("invalid minimum element in 'bloodmeal': " + std::to_string(mu));
      Rcpp::stop(msg);
    }

    // last end time becomes next beginning time
    t0 = t1;

  }

};











/* --------------------------------------------------------------------------------
#   run (main)
-------------------------------------------------------------------------------- */

// [[Rcpp::export]]
Rcpp::NumericMatrix test_mosquitoes(const Rcpp::NumericVector parameters, const double X, const int SV, const int IV, const double dt, const double tmax){

  const int outsize = 1E5;
  mosypop_ptr mosypop = make_mosypop(parameters,SV,IV,outsize);

  // some parameters
  double a = parameters["a"];
  double c = parameters["c"];

  // 1st iteration
  double clock{0.0};

  // subsequent iteration
  while(clock < tmax){

    // iterate the mosquitoes
    run_mosypop(mosypop, clock, dt);

    // compute bloodmeal times (the BLOODMEAL MODULE)

    // start of this piecewise trajectory
    queue_tuple t0_block = mosypop->Sv_trace.top();
    mosypop->Sv_trace.pop();
    // end of this piecewise trajectory
    queue_tuple t1_block;
    // iterate over parts
    // double ddt{0.};
    while(!mosypop->Sv_trace.empty()){
      // grab the end
      t1_block = mosypop->Sv_trace.top();
      mosypop->Sv_trace.pop();

      // dt: length of trajectory, S: the susceptible mosy over this traj
      double t1 = std::get<1>(t1_block);
      double t0 = std::get<1>(t0_block);
      double dt = t1 - t0;
      double S = std::get<0>(t0_block);

      // intensity for this piecewise block
      double intensity = a * c * X * S * dt;
      int bites = R::rpois(intensity);

      // add the bites to the mosquito for next time
      if(bites > 0){
        double btime = R::runif(t0,t1);
        mosypop->H2M_bites.push(btime);
      }

      // set endpoint to startpoint for next iteration
      t0_block = t1_block;
    }
    // Rcpp::Rcout << "done queueing bites, total dt accumluated: " << ddt << " \n";

    // clear out unused trace
    while(!mosypop->Iv_trace.empty()){
      mosypop->Iv_trace.pop();
    }

    // increment tmax before next iteration
    clock += dt;
  }

  int k = mosypop->t_hist.size();
  Rcpp::NumericMatrix out(k,3);
  Rcpp::colnames(out) = Rcpp::CharacterVector::create("time","S","I");
  for(int j=0; j<k; j++){
    out(j,0) = mosypop->t_hist[j];
    out(j,1) = mosypop->S_hist[j];
    out(j,2) = mosypop->I_hist[j];
  }
  return out;
};


// [[Rcpp::export]]
Rcpp::NumericMatrix test_humans(const Rcpp::NumericVector parameters, const int SH, const int IH, const double IV, const double dt, const double tmax){

  const int outsize = 1E5;
  humanpop_ptr humanpop = make_humanpop(parameters,SH,IH,outsize);

  // some parameters
  double a = parameters["a"];
  double b = parameters["b"];

  // 1st iteration
  double clock{0.0};

  // subsequent iteration
  while(clock < tmax){

    // iterate the mosquitoes
    run_humanpop(humanpop, clock, dt);

    // compute bloodmeal times (the BLOODMEAL MODULE)

    // start of this piecewise trajectory
    queue_tuple t0_block = humanpop->X_trace.top();
    humanpop->X_trace.pop();
    // end of this piecewise trajectory
    queue_tuple t1_block;
    // iterate over parts
    // double ddt{0.};
    while(!humanpop->X_trace.empty()){
      // grab the end
      t1_block = humanpop->X_trace.top();
      humanpop->X_trace.pop();

      // dt: length of trajectory, S: the susceptible mosy over this traj
      double t1 = std::get<1>(t1_block);
      double t0 = std::get<1>(t0_block);
      double dt = t1 - t0;
      double X = std::get<0>(t0_block);
      // Rcpp::Rcout <<"X: " << X << "---\n";

      // intensity for this piecewise block
      double intensity = a * b * (1. - X) * IV * dt;
      int bites = R::rpois(intensity);

      // add the bites to the mosquito for next time
      if(bites > 0){
        double btime = R::runif(t0,t1);
        humanpop->M2H_bites.push(btime);
      }

      // set endpoint to startpoint for next iteration
      t0_block = t1_block;
    }

    // increment tmax before next iteration
    clock += dt;
  }

  int k = humanpop->t_hist.size();
  Rcpp::NumericMatrix out(k,3);
  Rcpp::colnames(out) = Rcpp::CharacterVector::create("time","S","I");
  for(int j=0; j<k; j++){
    out(j,0) = humanpop->t_hist[j];
    out(j,1) = humanpop->S_hist[j];
    out(j,2) = humanpop->I_hist[j];
  }
  return out;
};
