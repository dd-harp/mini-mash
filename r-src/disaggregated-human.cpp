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

  // elements for MNRM of internal dynamics
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
  double X = static_cast<double>(humanpop->I) / (static_cast<double>(humanpop->S) + static_cast<double>(humanpop->I));
  // double X = static_cast<double>(humanpop->S) / (humanpop->N);
  humanpop->X_trace.push(queue_tuple(X,t0));

  double tmax{t0+dt};

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
    auto min_elem = std::min_element(humanpop->delta_t.begin(), humanpop->delta_t.end());
    int mu = std::distance(humanpop->delta_t.begin(), min_elem);

    // update putative time
    double delta = humanpop->delta_t[mu];
    double tsamp = humanpop->tnow + delta;
    // if overrun, advance time to tmax and adjust the Poisson process for the fast forward
    if(tsamp > tmax){
      double remaining = tmax - humanpop->tnow;
      humanpop->Tk += humanpop->ak * remaining;
      X = static_cast<double>(humanpop->I) / (static_cast<double>(humanpop->S) + static_cast<double>(humanpop->I));
      // double X = static_cast<double>(humanpop->S) / (humanpop->N);
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
    // double X = static_cast<double>(humanpop->S) / (humanpop->N);
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
#   run (main)
-------------------------------------------------------------------------------- */


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

      // intensity for this piecewise block
      double intensity = a * b * (1. - X) * IV * dt;
      // double intensity = a * b * X * IV * dt;
      int bites = R::rpois(intensity);

      // add the bites to the human for next time
      if(bites > 0){
        for(int k=0; k<bites; k++){
          double btime = R::runif(t0,t1);
          humanpop->M2H_bites.push(btime);
        }
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
