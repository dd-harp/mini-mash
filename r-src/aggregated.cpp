/* --------------------------------------------------------------------------------
#
#   aggregated malaria transmission model
#   Sean L. Wu (slwu89@berkeley.edu)
#   November 2020
#
-------------------------------------------------------------------------------- */

// [[Rcpp::plugins(cpp14)]]

#include <vector>
#include <array>
#include <string>
#include <limits>
#include <queue>
#include <functional>

#include <Rcpp.h>

enum events {
  emerge = 0,
  SV_die = 1,
  IV_die = 2,
  m_inf = 3,
  h_inf = 4,
  h_rec = 5,
  m_EIP = 6,
  h_LEP = 7
};

enum states {
  SH = 0,
  IH = 1,
  SV = 2,
  IV = 3
};

const static double eps = 1.E-9;
const static double infinity = std::numeric_limits<double>::infinity();



// [[Rcpp::export]]
Rcpp::NumericMatrix pfsim_aggregated(
  const double tmax,
  const int SH,
  const int IH,
  const int SV,
  const int IV,
  const Rcpp::NumericVector parameters,
  const bool verbose
){

  // extract parameters
  double g = parameters["g"];
  double lambda = parameters["lambda"];
  double a = parameters["a"];
  double b = parameters["b"];
  double c = parameters["c"];
  double r = parameters["r"];
  double EIP = parameters["EIP"];
  double LEP = parameters["LEP"];
  double surv_EIP = exp(-g * EIP);

  // state
  std::array<int,4> X{0};
  X[states::SH] = SH;
  X[states::IH] = IH;
  X[states::SV] = SV;
  X[states::IV] = IV;

  // next time each event fires
  std::array<double,8> delta_t{infinity};

  // Pk: next internal firing time of the Poisson process
  std::array<double,6> Pk{0.0};

  // Tk: current internal time of the Poisson process (integrated propensity)
  std::array<double,6> Tk{0.0};

  // ak: propensity functions
  std::array<double,6> ak{0.0};

  // delay effects
  std::priority_queue<double,std::vector<double>,std::greater<double> > m_EIP;
  std::priority_queue<double,std::vector<double>,std::greater<double> > h_LEP;
  m_EIP.push(infinity);
  h_LEP.push(infinity);

  double tnow{0.};
  const int outsize = 1E5;
  int i{1};

  // trajectory of the process
  std::vector<double> t_hist;
  std::vector<double> SH_hist;
  std::vector<double> IH_hist;
  std::vector<double> SV_hist;
  std::vector<double> IV_hist;
  t_hist.reserve(outsize);
  SH_hist.reserve(outsize);
  IH_hist.reserve(outsize);
  SV_hist.reserve(outsize);
  IV_hist.reserve(outsize);
  t_hist.push_back(tnow);
  SH_hist.push_back(X[states::SH]);
  IH_hist.push_back(X[states::IH]);
  SV_hist.push_back(X[states::SV]);
  IV_hist.push_back(X[states::IV]);


  // calculate initial propensities
  ak[events::emerge] = lambda;
  ak[events::SV_die] = g * static_cast<double>(X[states::SV]);
  ak[events::IV_die] = g * static_cast<double>(X[states::IV]);
  double x = static_cast<double>(X[states::IH]) / (static_cast<double>(X[states::IH]) + static_cast<double>(X[states::SH]));
  ak[events::m_inf] = a * c * x * static_cast<double>(X[states::SV]);
  ak[events::h_inf] = a * b * static_cast<double>(X[states::IV]) * (1. - x);
  ak[events::h_rec] = r * static_cast<double>(X[states::IH]);

  // draw internal firing times
  for(int j=0; j<6; j++){
    Pk[j] = log(1. / R::runif(0., 1.));
  }

  if(verbose){
    Rcpp::Rcout << " --- beginning simulation --- \n";
  }

  // main loop
  while(tnow < tmax){

    if(verbose){
      if(i % 10000 == 0){
        Rcpp::Rcout << " --- simulated " << i << " events with remaining time: " << (tmax - tnow) << " --- \n";
      }
    }

    // calculate absolute times to fire (wall clock)
    for(int j=0; j<6; j++){
      delta_t[j] = (Pk[j] - Tk[j]) / ak[j];
    }
    delta_t[events::m_EIP] = m_EIP.top() - tnow;
    delta_t[events::h_LEP] = h_LEP.top() - tnow;

    // find minimum
    auto min_elem = std::min_element(delta_t.begin(), delta_t.end());
    int min_mu = std::distance(delta_t.begin(), min_elem);
    auto mu = static_cast<events>(min_mu);

    // update time
    double delta = delta_t[min_mu];
    tnow += delta;
    if(tnow > tmax){
      break;
    }

    // update system
    if(mu == events::emerge){
      X[states::SV] += 1;
    } else if(mu == events::SV_die){
      X[states::SV] -= 1;
    } else if(mu == events::IV_die){
      X[states::IV] -= 1;
    } else if(mu == events::m_inf){
      X[states::SV] -= 1;
      // does this incubating mosquito survive the EIP?
      if(R::runif(0., 1.) < surv_EIP){
        m_EIP.push(tnow + EIP);
      }
    } else if(mu == events::h_inf){
      X[states::SH] -= 1;
      h_LEP.push(tnow + LEP);
    } else if(mu == events::h_rec){
      X[states::IH] -= 1;
      X[states::SH] += 1;
    } else if(mu == events::m_EIP){
      m_EIP.pop();
      X[states::IV] += 1;
    } else if(mu == events::h_LEP){
      h_LEP.pop();
      X[states::IH] += 1;
    } else {
      std::string msg("invalid minimum element: " + std::to_string(min_mu));
      Rcpp::stop(msg);
    }

    // update Tk
    for(int j=0; j<6; j++){
      Tk[j] += ak[j] * delta;
    }

    // update Pk[mu]
    if(mu != events::h_LEP && mu != events::m_EIP){
      Pk[min_mu] += log(1. / R::runif(0., 1.));
    }

    // recalculate propentities
    ak[events::emerge] = lambda;
    ak[events::SV_die] = g * static_cast<double>(X[states::SV]);
    ak[events::IV_die] = g * static_cast<double>(X[states::IV]);
    double x = static_cast<double>(X[states::IH]) / (static_cast<double>(X[states::IH]) + static_cast<double>(X[states::SH]));
    ak[events::m_inf] = a * c * x * static_cast<double>(X[states::SV]);
    ak[events::h_inf] = a * b * static_cast<double>(X[states::IV]) * (1. - x);
    ak[events::h_rec] = r * static_cast<double>(X[states::IH]);

    // track output
    t_hist.push_back(tnow);
    SH_hist.push_back(X[states::SH]);
    IH_hist.push_back(X[states::IH]);
    SV_hist.push_back(X[states::SV]);
    IV_hist.push_back(X[states::IV]);
    i++;

    // check break conditions
    bool ak_zero = std::all_of(ak.begin(), ak.end(), [](double a){return a < eps;});
    bool lag_inf = m_EIP.top() == infinity && h_LEP.top() == infinity;
    if(ak_zero && lag_inf){
      break;
    }

  }

  int k = t_hist.size();
  Rcpp::NumericMatrix out(k,5);
  Rcpp::colnames(out) = Rcpp::CharacterVector::create("time","SH","IH","SV","IV");
  for(int j=0; j<k; j++){
    out(j,0) = t_hist[j];
    out(j,1) = SH_hist[j];
    out(j,2) = IH_hist[j];
    out(j,3) = SV_hist[j];
    out(j,4) = IV_hist[j];
  }
  return out;

};
