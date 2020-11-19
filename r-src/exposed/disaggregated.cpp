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

// iterate the mosquitoes
void run_mosypop(mosypop_ptr& mosypop, double t0, double dt){

  // push initial state into traces
  if(!mosypop->Sv_trace.empty() || !mosypop->Iv_trace.empty()){
    Rcpp::stop("'Sv_trace' and 'Iv_trace' should always be empty at start of time step \n");
  }
  mosypop->Sv_trace.emplace(queue_tuple(mosypop->S,t0));
  mosypop->Iv_trace.emplace(queue_tuple(mosypop->I,t0));

  double tmax{t0+dt};

  // get bites from the previous time step
  while(!mosypop->H2M_bites.empty()){
    // take out one S mosquito (they move to E)
    mosypop->S -= 1;
    mosypop->E += 1;
    // queue the future EIP completion event
    mosypop->M_inf.emplace(*mosypop->H2M_bites.begin() + mosypop->EIP);
    mosypop->H2M_bites.erase(mosypop->H2M_bites.begin());
  }

  // simulate dynamics over this time step
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
    // if overrun, advance time to tmax and adjust the Poisson processes for the fast-forward
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

    // track output
    mosypop->t_hist.push_back(mosypop->tnow);
    mosypop->S_hist.push_back(mosypop->S);
    mosypop->E_hist.push_back(mosypop->E);
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

  // elements for MNRM of internal dynamics
  std::array<double,2> delta_t{infinity};
  double Pk{0.};
  double Tk{0.};
  double ak{0.};

  // integrated prevalence over [t0,t0+dt)
  queue_trace X_trace;

  humanpop_str(){Rcpp::Rcout << "'humanpop_str' ctor called at " << this << std::endl;};
  ~humanpop_str(){Rcpp::Rcout << "'humanpop_str' dtor called at " << this << std::endl;};

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

  humanpop->H_inf.emplace(infinity);

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
  humanpop->X_trace.emplace(queue_tuple(X,t0));

  // get bites from the previous time step
  while(!humanpop->M2H_bites.empty()){
    // take out one S person (they move to E)
    humanpop->S -= 1;
    // push the future infection
    humanpop->H_inf.emplace(*humanpop->M2H_bites.begin() + humanpop->LEP);
    humanpop->M2H_bites.erase(humanpop->M2H_bites.begin());
  }

  // simulate dynamics over this time step
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
      X = static_cast<double>(humanpop->I) / (static_cast<double>(humanpop->S) + static_cast<double>(humanpop->I));
      humanpop->X_trace.emplace(queue_tuple(X,tmax));
      humanpop->tnow = tmax;
      break;
    }
    humanpop->tnow = tsamp;

    // update system
    if(mu == 0){
      humanpop->I -= 1;
      humanpop->S += 1;
    } else if(mu == 1){
      humanpop->H_inf.erase(humanpop->H_inf.begin());
      humanpop->I += 1;
    } else {
      std::string msg("invalid minimum element in 'run_mosypop': " + std::to_string(mu));
      Rcpp::stop(msg);
    }

    X = static_cast<double>(humanpop->I) / (static_cast<double>(humanpop->S) + static_cast<double>(humanpop->I));
    humanpop->X_trace.emplace(queue_tuple(X,humanpop->tnow));

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

// void bloodmeal(const Rcpp::NumericVector parameters, humanpop_ptr& humanpop, mosypop_ptr& mosypop){
//
//   double a = parameters["a"];
//   double c = parameters["c"];
//   double b = parameters["b"];
//
//   // H2M bites occur with intensity (acX)S_{V}
//   // M2H bites occur with intensity ab(1-X)I_{V}
//   double H2M_intensity, M2H_intensity;
//   int H2M_bites, M2H_bites;
//
//   // length of intervals
//   double t1,t0,dt;
//   double S,I,X;
//
//   // next state change
//   std::array<double,3> t1_array{0.};
//
//   // beginning of piecewise trajectories
//   queue_tuple t0_SV = *mosypop->Sv_trace.begin();
//   queue_tuple t0_IV = *mosypop->Iv_trace.begin();
//   queue_tuple t0_X = *humanpop->X_trace.begin();
//
//   mosypop->Sv_trace.erase(mosypop->Sv_trace.begin());
//   mosypop->Iv_trace.erase(mosypop->Iv_trace.begin());
//   humanpop->X_trace.erase(humanpop->X_trace.begin());
//
//   t0 = std::get<1>(t0_SV);
//
//   // end of the piecewise trajectories
//   queue_tuple t1_SV;
//   queue_tuple t1_IV;
//   queue_tuple t1_X;
//
//   // compute over the TWICE step
//   while(!mosypop->Sv_trace.empty() || !mosypop->Iv_trace.empty() || !humanpop->X_trace.empty()){
//
//     H2M_intensity = 0.;
//     M2H_intensity = 0.;
//
//     // if not empty, get the next state change in each trajectory
//     std::fill(t1_array.begin(),t1_array.end(),infinity);
//     if(!mosypop->Sv_trace.empty()){
//       t1_SV = *mosypop->Sv_trace.begin();
//       mosypop->Sv_trace.erase(mosypop->Sv_trace.begin());
//       t1_array[0] = std::get<1>(t1_SV);
//     }
//     if(!mosypop->Iv_trace.empty()){
//       t1_IV = *mosypop->Iv_trace.begin();
//       mosypop->Iv_trace.erase(mosypop->Iv_trace.begin());
//       t1_array[1] = std::get<1>(t1_IV);
//     }
//     if(!humanpop->X_trace.empty()){
//       t1_X = *humanpop->X_trace.begin();
//       humanpop->X_trace.erase(humanpop->X_trace.begin());
//       t1_array[2] = std::get<1>(t1_X);
//     }
//
//     // find which trajectory changes next
//     auto min_elem = std::min_element(t1_array.begin(), t1_array.end());
//     int mu = std::distance(t1_array.begin(), min_elem);
//
//     // compute length of this piecewise interval [t0,d1)
//     t1 = t1_array[mu];
//     dt = t1 - t0;
//
//     // state values at t0
//     S = std::get<0>(t0_SV);
//     I = std::get<0>(t0_IV);
//     X = std::get<0>(t0_X);
//
//     // intensity over the interval
//     H2M_intensity = a * c * X * S * dt;
//     M2H_intensity = a * b * (1. - X) * I * dt;
//
//     // sample values
//     H2M_bites = R::rpois(H2M_intensity);
//     M2H_bites = R::rpois(M2H_intensity);
//
//     // add the bites to the mosquito for next time
//     if(H2M_bites > 0){
//       for(int k=0; k<H2M_bites; k++){
//         double btime = R::runif(t0,t1);
//         mosypop->H2M_bites.emplace(btime);
//       }
//     }
//
//     // add the bites to the human for next time
//     if(M2H_bites > 0){
//       for(int k=0; k<M2H_bites; k++){
//         double btime = R::runif(t0,t1);
//         humanpop->M2H_bites.emplace(btime);
//       }
//     }
//
//     // state change that happened first becomes new starting point
//     if(mu==0){
//       t0_SV = t1_SV;
//     } else if(mu==1){
//       t0_IV = t1_IV;
//     } else if(mu==2){
//       t0_X = t1_X;
//     } else{
//       std::string msg("invalid minimum element in 'bloodmeal': " + std::to_string(mu));
//       Rcpp::stop(msg);
//     }
//
//     // last end time becomes next beginning time
//     t0 = t1;
//   }
//
// };

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
  queue_tuple t0_SV = *mosypop->Sv_trace.begin();
  queue_tuple t0_IV = *mosypop->Iv_trace.begin();
  queue_tuple t0_X = *humanpop->X_trace.begin();

  mosypop->Sv_trace.erase(mosypop->Sv_trace.begin());
  mosypop->Iv_trace.erase(mosypop->Iv_trace.begin());
  humanpop->X_trace.erase(humanpop->X_trace.begin());

  t0 = std::get<1>(t0_SV);

  // end of the piecewise trajectories
  queue_tuple t1_SV = *mosypop->Sv_trace.begin();
  queue_tuple t1_IV = *mosypop->Iv_trace.begin();
  queue_tuple t1_X = *humanpop->X_trace.begin();

  mosypop->Sv_trace.erase(mosypop->Sv_trace.begin());
  mosypop->Iv_trace.erase(mosypop->Iv_trace.begin());
  humanpop->X_trace.erase(humanpop->X_trace.begin());

  t1_array[0] = std::get<1>(t1_SV);
  t1_array[1] = std::get<1>(t1_IV);
  t1_array[2] = std::get<1>(t1_X);

  // compute over the TWICE step
  while( std::all_of(t1_array.begin(), t1_array.end(), [](const double y){return y < infinity;}) ){

    H2M_intensity = 0.;
    M2H_intensity = 0.;

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
      for(int k=0; k<H2M_bites; k++){
        double btime = R::runif(t0,t1);
        mosypop->H2M_bites.emplace(btime);
      }
    }

    // add the bites to the human for next time
    if(M2H_bites > 0){
      for(int k=0; k<M2H_bites; k++){
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
      t0_X = t1_X;
      if(!humanpop->X_trace.empty()){
        t1_X = *humanpop->X_trace.begin();
        humanpop->X_trace.erase(humanpop->X_trace.begin());
        t1_array[2] = std::get<1>(t1_X);
      } else {
        t1_array[2] = infinity;
      }
    } else {
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
Rcpp::List run_miniMASH(const Rcpp::NumericVector parameters, const Rcpp::IntegerVector y0, const double dt, const double tmax){

  int SH = y0[0];
  int IH = y0[1];
  int SV = y0[2];
  int IV = y0[3];

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
    if(i % 100 == 0){
      Rcpp::Rcout << " completed TWICE step " << i << " of " << tsteps << "\n";
    }
    clock += dt;
  }
  Rcpp::Rcout << "--- ending simulation ---\n";


  int kh = humanpop->t_hist.size();
  Rcpp::NumericMatrix outh(kh,3);
  Rcpp::colnames(outh) = Rcpp::CharacterVector::create("time","S","I");
  for(int j=0; j<kh; j++){
    outh(j,0) = humanpop->t_hist[j];
    outh(j,1) = humanpop->S_hist[j];
    outh(j,2) = humanpop->I_hist[j];
  }

  int km = mosypop->t_hist.size();
  Rcpp::NumericMatrix outm(km,3);
  Rcpp::colnames(outm) = Rcpp::CharacterVector::create("time","S","I");
  for(int j=0; j<km; j++){
    outm(j,0) = mosypop->t_hist[j];
    outm(j,1) = mosypop->S_hist[j];
    outm(j,2) = mosypop->I_hist[j];
  }

  return Rcpp::List::create(
    Rcpp::Named("mosquito") = outm,
    Rcpp::Named("human") = outh
  );
}
