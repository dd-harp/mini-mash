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
#include <limits>
#include <set>
#include <algorithm>


#include <Rcpp.h>

enum events {
  emerge = 0,
  SV_die = 1,
  EV_die = 2,
  IV_die = 3,
  m_inf = 4,
  h_inf = 5,
  h_rec = 6,
  m_EIP = 7,
  h_LEP = 8
};

enum states {
  SH = 0,
  EH = 1,
  IH = 2,
  SV = 3,
  EV = 4,
  IV = 5
};

// const static double eps = 1.E-9;
const static double infinity = std::numeric_limits<double>::infinity();



// [[Rcpp::export]]
Rcpp::NumericMatrix pfsim_aggregated(
  const double tmax,
  const int SH,
  const int EH,
  const int IH,
  const int SV,
  const int EV,
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

  // state
  std::array<int,6> X{0};
  X[states::SH] = SH;
  X[states::EH] = EH;
  X[states::IH] = IH;
  X[states::SV] = SV;
  X[states::EV] = EV;
  X[states::IV] = IV;

  // next time each event fires
  std::array<double,9> delta_t{infinity};

  // Pk: next internal firing time of the Poisson process
  std::array<double,7> Pk{0.0};

  // Tk: current internal time of the Poisson process (integrated propensity)
  std::array<double,7> Tk{0.0};

  // ak: propensity functions
  std::array<double,7> ak{0.0};

  // delay effects
  std::set<double, std::less<double>> m_EIP;
  std::set<double, std::less<double>> h_LEP;
  m_EIP.emplace(infinity);
  h_LEP.emplace(infinity);

  double tnow{0.};
  const int outsize = 1E5;
  int i{1};

  // trajectory of the process
  std::vector<double> t_hist;
  std::vector<double> SH_hist;
  std::vector<double> EH_hist;
  std::vector<double> IH_hist;
  std::vector<double> SV_hist;
  std::vector<double> EV_hist;
  std::vector<double> IV_hist;
  t_hist.reserve(outsize);
  SH_hist.reserve(outsize);
  EH_hist.reserve(outsize);
  IH_hist.reserve(outsize);
  SV_hist.reserve(outsize);
  EV_hist.reserve(outsize);
  IV_hist.reserve(outsize);
  t_hist.push_back(tnow);
  SH_hist.push_back(X[states::SH]);
  EH_hist.push_back(X[states::EH]);
  IH_hist.push_back(X[states::IH]);
  SV_hist.push_back(X[states::SV]);
  EV_hist.push_back(X[states::EV]);
  IV_hist.push_back(X[states::IV]);

  // calculate initial propensities
  double NH = static_cast<double>(X[states::SH]) + static_cast<double>(X[states::EH]) + static_cast<double>(X[states::IH]);
  double x = static_cast<double>(X[states::IH]) / NH;
  double z = static_cast<double>(X[states::SH]) / NH;

  ak[events::emerge] = lambda;
  ak[events::SV_die] = g * static_cast<double>(X[states::SV]);
  ak[events::EV_die] = g * static_cast<double>(X[states::EV]);
  ak[events::IV_die] = g * static_cast<double>(X[states::IV]);
  ak[events::m_inf] = a * c * x * static_cast<double>(X[states::SV]);
  ak[events::h_inf] = a * b * z * static_cast<double>(X[states::IV]);
  ak[events::h_rec] = r * static_cast<double>(X[states::IH]);

  // draw internal firing times
  for(int j=0; j<7; j++){
    Pk[j] = log(1. / R::runif(0., 1.));
  }
  
  // draw initial EH->IH and EV->IV times for delayed individuals
  for (int i=0; i<EH; i++) {
    double t_inf = R::runif(-LEP,0.0);
    h_LEP.emplace(t_inf + LEP);
  }
  
  for (int i=0; i<EV; i++) {
    double t_inf = R::runif(-EIP,0.0);
    m_EIP.emplace(t_inf + EIP);
  }

  // main loop
  if(verbose){
    Rcpp::Rcout << " --- beginning simulation --- \n";
  }
  
  while(tnow < tmax){

    if(verbose){
      if(i % 10000 == 0){
        Rcpp::Rcout << " --- simulated " << i << " events with remaining time: " << (tmax - tnow) << " --- \n";
      }
    }

    // calculate absolute times to fire (wall clock)
    for(int j=0; j<7; j++){
      delta_t[j] = (Pk[j] - Tk[j]) / ak[j];
    }
    delta_t[events::m_EIP] = *(m_EIP.begin()) - tnow;
    delta_t[events::h_LEP] = *(h_LEP.begin()) - tnow;

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
    } else if(mu == events::EV_die){
      // delete a random element from the EIP set (not the end, which is Inf)
      int relem = R::runif(0.0, m_EIP.size()-1);
      auto it = m_EIP.begin();
      std::advance(it, relem);
      m_EIP.erase(it);
      X[states::EV] -= 1;
      assert(m_EIP.size()-1 == X[states::EV]);
    } else if(mu == events::IV_die){
      X[states::IV] -= 1;
    } else if(mu == events::m_inf){
      X[states::SV] -= 1;
      X[states::EV] += 1;
      m_EIP.emplace(tnow + EIP);
      assert(m_EIP.size()-1 == X[states::EV]);
    } else if(mu == events::h_inf){
      X[states::SH] -= 1;
      X[states::EH] += 1;
      h_LEP.emplace(tnow + LEP);
      assert(h_LEP.size()-1 == X[states::EH]);
    } else if(mu == events::h_rec){
      X[states::IH] -= 1;
      X[states::SH] += 1;
    } else if(mu == events::m_EIP){
      m_EIP.erase(m_EIP.begin());
      X[states::EV] -= 1;
      X[states::IV] += 1;
      assert(m_EIP.size()-1 == X[states::EV]);
    } else if(mu == events::h_LEP){
      h_LEP.erase(h_LEP.begin());
      X[states::EH] -= 1;
      X[states::IH] += 1;
      assert(h_LEP.size()-1 == X[states::EH]);
    } else {
      std::string msg("invalid minimum element: " + std::to_string(min_mu));
      Rcpp::stop(msg);
    }

    // update Tk
    for(int j=0; j<7; j++){
      Tk[j] += ak[j] * delta;
    }

    // update Pk[mu]
    if(mu != events::h_LEP && mu != events::m_EIP){
      Pk[min_mu] += log(1. / R::runif(0., 1.));
    }

    // recalculate propentities
    double NH = static_cast<double>(X[states::SH]) + static_cast<double>(X[states::EH]) + static_cast<double>(X[states::IH]);
    double x = static_cast<double>(X[states::IH]) / NH;
    double z = static_cast<double>(X[states::SH]) / NH;

    ak[events::emerge] = lambda;
    ak[events::SV_die] = g * static_cast<double>(X[states::SV]);
    ak[events::EV_die] = g * static_cast<double>(X[states::EV]);
    ak[events::IV_die] = g * static_cast<double>(X[states::IV]);
    ak[events::m_inf] = a * c * x * static_cast<double>(X[states::SV]);
    ak[events::h_inf] = a * b * z * static_cast<double>(X[states::IV]);
    ak[events::h_rec] = r * static_cast<double>(X[states::IH]);

    assert(std::all_of(X.begin(), X.end(), [](int i){return i >= 0;}));

    // track output
    t_hist.push_back(tnow);
    SH_hist.push_back(X[states::SH]);
    EH_hist.push_back(X[states::EH]);
    IH_hist.push_back(X[states::IH]);
    SV_hist.push_back(X[states::SV]);
    EV_hist.push_back(X[states::EV]);
    IV_hist.push_back(X[states::IV]);
    i++;

  }

  int k = t_hist.size();
  Rcpp::NumericMatrix out(k,7);
  Rcpp::colnames(out) = Rcpp::CharacterVector::create("time","SH","EH","IH","SV","EV","IV");
  for(int j=0; j<k; j++){
    out(j,0) = t_hist[j];
    out(j,1) = SH_hist[j];
    out(j,2) = EH_hist[j];
    out(j,3) = IH_hist[j];
    out(j,4) = SV_hist[j];
    out(j,5) = EV_hist[j];
    out(j,6) = IV_hist[j];
  }
  return out;

};
