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
#include <memory>
#include <queue>

#include <Rcpp.h>

const static double infinity = std::numeric_limits<double>::infinity();

using queue = std::priority_queue<double,std::vector<double>,std::greater<double> >;

enum class intensity_type {SV,IV,X};

// if have these, plus the end time of the step, can calculate the integrated intensity
typedef struct integrated_intensity {

  double t0;
  double val;
  intensity_type type;

} integrated_intensity;

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

  // integrated susceptible (Sv) biting
  // integrated infectious (Iv) biting
  double Sv_bite;
  double Iv_bite;

  // mosypop_str(){};
  // ~mosypop_str(){};

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

  double tmax{t0+dt};
  double P = exp(-mosypop->g*mosypop->EIP);

  // integrated susceptible (Sv) biting
  // integrated infectious (Iv) biting
  mosypop->Sv_bite = 0.;
  mosypop->Iv_bite = 0.;

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
    // if overrun, advance time to tmax and adjust the Poisson processes
    if(tsamp > tmax){
      double remaining = tmax - mosypop->tnow;
      for(int j=0; j<3; j++){
        mosypop->Tk[j] += mosypop->ak[j] * remaining;
      }
      mosypop->Sv_bite += mosypop->a * mosypop->S * remaining;
      mosypop->Iv_bite += mosypop->a * mosypop->I * remaining;
      mosypop->tnow = tmax;
      break;
    }

    // update time after updating integrated biting
    mosypop->Sv_bite += mosypop->a * mosypop->S * delta;
    mosypop->Iv_bite += mosypop->a * mosypop->I * delta;
    mosypop->tnow = tsamp;

    // update system
    if(mu == 1){
      mosypop->S += 1;
    } else if(mu == 2){
      mosypop->S -= 1;
    } else if(mu == 3){
      mosypop->I -= 1;
    } else if(mu == 4){
      mosypop->M_inf.pop();
      mosypop->I += 1;
    } else {
      std::string msg("invalid minimum element in 'run_mosypop': " + std::to_string(mu));
      Rcpp::stop(msg);
    }




  }

  // drawk internal firing times


}
