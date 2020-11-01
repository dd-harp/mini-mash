/* --------------------------------------------------------------------------------
#
#   aggregated malaria transmission model
#   Sean L. Wu (slwu89@berkeley.edu)
#   October 2020
#
-------------------------------------------------------------------------------- */

// [[Rcpp::plugins(cpp14)]]

#include <vector>
#include <array>
#include <limits>
#include <queue>
#include <functional>

#include <Rcpp.h>

const static double eps = 1.E-9;
const static double infinity = std::numeric_limits<double>::infinity();


//' Simulate non-Markovian SIR Model via Modified Next Reaction Method (MNRM)
//'
//' In this non-Markovian variant of the SIR model, the infectious period has a Gamma distribution.
//'
//' Sample a trajectory from the Markovian SIR model using the MNRM algorithm presented in:
//'   * Anderson, D. F. (2007). A modified next reaction method for simulating chemical systems with time dependent propensities and delays. Journal of Chemical Physics, 127(21). \url{https://doi.org/10.1063/1.2799998}
//'
//' @param tmax maximum time of simulation
//' @param S initial number of susceptible individuals
//' @param I initial number of infected & infectious individuals
//' @param R initial number of recovered individuals
//' @param beta the product of transmission probability and contact rate
//' @param gamma_shape shape parameter of Gamma distributed infectious period
//' @param gamma_scale scale parameter of Gamma distributed infectious period
//' @param verbose print extra information?
//'
//' @return a matrix
//'
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix pfsim_aggregated(
  const double tmax,
  const int S,
  const int I,
  const int R,
  const double beta,
  const double gamma_shape,
  const double gamma_scale,
  const bool verbose
){

  // system state
  std::array<int,3> X{S,I,R};

  double tnow{0.};
  const int outsize = 1E5;
  int i{1};

  // trajectory of the process
  std::vector<double> t_hist;
  std::vector<double> S_hist;
  std::vector<double> I_hist;
  std::vector<double> R_hist;
  t_hist.reserve(outsize);
  S_hist.reserve(outsize);
  I_hist.reserve(outsize);
  R_hist.reserve(outsize);
  t_hist.push_back(tnow);
  S_hist.push_back(X[0]);
  I_hist.push_back(X[1]);
  R_hist.push_back(X[2]);

  // 1. initialize
  double Pk{0.};        // next internal firing time of Poisson process (Pk > Tk)
  double Tk{0.};        // internal time of Poisson process (integrated propensity)
  double delta_t{0.};   // absolute time to fire
  double ak{0.};        // propensity function
  std::priority_queue<double,std::vector<double>,std::greater<double> > sk; // delayed recovery
  sk.push(infinity);

  double Delta{0.};
  bool   delay{false};

  // 2. calculate propensities
  ak = beta * static_cast<double>(X[0]) * static_cast<double>(X[1]);

  // initial completions of delay recovery for I(0) individuals
  for(int k=0; k<I; k++){
    sk.push(R::rgamma(gamma_shape,gamma_scale));
  }

  // 3-4: draw internal jump times
  Pk = log(1. / R::runif(0.,1.));

  if(verbose){
    Rcpp::Rcout << " --- beginning simulation --- \n";
  }

  while(tnow < tmax){

    if(verbose){
      if(i % 100 == 0){
        Rcpp::Rcout << " --- simulated " << i << " reactions at time " << tnow << " --- \n";
      }
    }

    // 5. set absolute times to fire
    delta_t = (Pk - Tk) / ak;

    // 6. find minimum
    if(delta_t < (sk.top() - tnow)){
      delay = false;
      Delta = delta_t;
    } else {
      delay = true;
      Delta = sk.top() - tnow;
    }

    tnow += Delta;
    // check if we can break early
    if(tnow > tmax){
      break;
    }

    // 7-10: update state
    if(delay){
      // I -> R
      X[1] -= 1;
      X[2] += 1;
      // remove the recovery time
      sk.pop();
    } else {
      // S -> I
      X[0] -= 1;
      X[1] += 1;
      // draw a recovery time
      sk.push(tnow + R::rgamma(gamma_shape,gamma_scale));
    }

    // 11. update Tk
    Tk += ak*Delta;

    // 12. update P_mu
    if(!delay){
      Pk += log(1. / R::runif(0.,1.));
    }

    // 13. recalculate propensities
    ak = beta * static_cast<double>(X[0]) * static_cast<double>(X[1]);

    // store output
    t_hist.push_back(tnow);
    S_hist.push_back(X[0]);
    I_hist.push_back(X[1]);
    R_hist.push_back(X[2]);
    i++;

    // check to see if we can end early
    if((ak < eps) && (sk.top() == infinity)){
      Rcpp::Rcout << " --- all propensities approximately zero; returning output early --- \n";
      int k = t_hist.size();
      Rcpp::NumericMatrix out(k,4);
      Rcpp::colnames(out) = Rcpp::CharacterVector::create("time","S","I","R");
      for(int j=0; j<k; j++){
        out(j,0) = t_hist[j];
        out(j,1) = S_hist[j];
        out(j,2) = I_hist[j];
        out(j,3) = R_hist[j];
      }
      return out;
    }
  }

  int k = t_hist.size();
  Rcpp::NumericMatrix out(k,4);
  Rcpp::colnames(out) = Rcpp::CharacterVector::create("time","S","I","R");
  for(int j=0; j<k; j++){
    out(j,0) = t_hist[j];
    out(j,1) = S_hist[j];
    out(j,2) = I_hist[j];
    out(j,3) = R_hist[j];
  }
  return out;
};
