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

const static double infinity = std::numeric_limits<double>::infinity();

// data structure
using queue = std::priority_queue<double,std::vector<double>,std::greater<double> >;

// data structure to pass traces
using queue_tuple = std::tuple<int,double>;
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



// mosquito population
typedef struct mosypop_str {

  // states
  int S;
  int E;
  int I;

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

using mosypop_ptr = std::unique_ptr<mosypop_str>;

// make the mosquitoes
mosypop_ptr make_mosypop(const Rcpp::NumericVector parameters, int SV, int IV){
  mosypop_ptr mosypop = std::make_unique<mosypop_str>();

  mosypop->S = SV;
  mosypop->E = 0;
  mosypop->I = IV;

  mosypop->tnow = 0.0;

  mosypop->H2M_bites.push(infinity);
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

  // get bites from the previous time step
  while(mosypop->H2M_bites.top() != infinity){
    // if that mosquito would survive the EIP, queue an infection event in the future
    if(R::runif(0.,1.) < P){
      mosypop->M_inf.push(mosypop->H2M_bites.top() + mosypop->EIP);
      mosypop->H2M_bites.pop();
    }
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
    if(mu == 1){
      mosypop->S += 1;
      mosypop->Sv_trace.push(queue_tuple(mosypop->S,mosypop->tnow));
    } else if(mu == 2){
      mosypop->S -= 1;
      mosypop->Sv_trace.push(queue_tuple(mosypop->S,mosypop->tnow));
    } else if(mu == 3){
      mosypop->I -= 1;
      mosypop->Iv_trace.push(queue_tuple(mosypop->I,mosypop->tnow));
    } else if(mu == 4){
      mosypop->M_inf.pop();
      mosypop->I += 1;
      mosypop->Iv_trace.push(queue_tuple(mosypop->I,mosypop->tnow));
    } else {
      std::string msg("invalid minimum element in 'run_mosypop': " + std::to_string(mu));
      Rcpp::stop(msg);
    }

    // update Tk


  }

  // draw internal firing times


}
