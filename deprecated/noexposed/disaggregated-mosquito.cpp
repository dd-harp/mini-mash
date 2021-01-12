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
        for(int k=0; k<bites; k++){
          double btime = R::runif(t0,t1);
          mosypop->H2M_bites.push(btime);
        }
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
