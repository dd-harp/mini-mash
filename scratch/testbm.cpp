
#include <set>
#include <algorithm>
#include <vector>
#include <array>

// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp14)]]

#include <Rcpp.h>

// Inf
static double infinity = std::numeric_limits<double>::infinity();

// store events (earliest/smallest value on top)
using queue = std::set<double,std::less<double>>;

// data structure to pass traces
using queue_tuple = std::tuple<double,double>; // 1st element is state, 2nd is time

struct queue_comp {
  bool operator() (const queue_tuple& a, const queue_tuple& b) const {
    return std::get<1>(a) < std::get<1>(b);
  }
};

using queue_trace = std::set<queue_tuple,queue_comp>;



// [[Rcpp::export]]
Rcpp::List test_H2M(
  const Rcpp::NumericMatrix Sv_mat,
  const Rcpp::NumericMatrix Ih_mat,
  const double NH,
  const double t0,
  const double dt,
  const Rcpp::NumericVector& parameters
){

  queue_trace Sv_trace;
  queue_trace Ih_trace;

  int i;

  for(i=0; i<Sv_mat.nrow(); i++){
    Sv_trace.emplace(queue_tuple(Sv_mat.at(i,0),Sv_mat.at(i,1)));
  }

  for(i=0; i<Ih_mat.nrow(); i++){
    Ih_trace.emplace(queue_tuple(Ih_mat.at(i,0),Ih_mat.at(i,1)));
  }

  // output
  std::vector<double> intensity;
  std::vector<double> left_lims;
  std::vector<int>    sampled_bites;

  double tmax(t0+dt);

  // extract parameters
  double a = parameters["a"];
  double c = parameters["c"];

  // H2M transmission (S->E) events in mosquitoes
  double H2M_lambda{0.};
  int H2M_count{0};
  int H2M_tot{0};

  double tl,tr; // bounds of current interval
  double SV,IH; // intensity is: (acX)SV

  std::array<double,2> tnext; // next interval

  // left limit of interval
  queue_tuple Sv_0 = *Sv_trace.begin();
  queue_tuple Ih_0 = *Ih_trace.begin();
  Sv_trace.erase(Sv_trace.begin());
  Ih_trace.erase(Ih_trace.begin());

  tl = std::get<1>(Sv_0);

  // right limit of interval
  queue_tuple Sv_1 = *Sv_trace.begin();
  queue_tuple Ih_1 = *Ih_trace.begin();
  Sv_trace.erase(Sv_trace.begin());
  Ih_trace.erase(Ih_trace.begin());

  tnext[0] = std::get<1>(Sv_1);
  tnext[1] = std::get<1>(Ih_1);

  // compute H2M point process over [t0,t0+dt)
  while(true){

    Rcpp::checkUserInterrupt();

    int mu;

    // the last interval (to tmax)
    if( std::all_of(tnext.begin(),tnext.end(), [](const double y){return y == infinity;}) ){

      mu = -1; // sign to break

      tr = tmax;

    } else {

      // find which trace changes next
      auto min_elem = std::min_element(tnext.begin(), tnext.end());
      mu = std::distance(tnext.begin(), min_elem);

      // compute length of this piecewise interval [tl,tr)
      tr = tnext[mu];

    }

    Rcpp::Rcout << " --- evaluating trajectory over [" << tl << "," << tr << "), ";

    double delta = tr - tl;

    // trajectory values at tl
    SV = std::get<0>(Sv_0);
    IH = std::get<0>(Ih_0);

    // intensity over [tl,tr)
    H2M_lambda = a * c * (IH/NH) * (SV - (double)H2M_tot) * delta;

    Rcpp::Rcout << " intensity evaluated as: " << a << " * " << c << " * " << (IH/NH) << " * " << SV << " * " << delta << " --- \n";

    // sample bites
    H2M_count = R::rpois(H2M_lambda);
    H2M_tot += H2M_count;

    // push output
    intensity.push_back(H2M_lambda);
    left_lims.push_back(tl);
    sampled_bites.push_back(H2M_count);

    // next change in Sv
    if(mu == 0){

      Sv_0 = Sv_1;

      if(!Sv_trace.empty()){

        Sv_1 = *Sv_trace.begin();
        Sv_trace.erase(Sv_trace.begin());
        tnext[0] = std::get<1>(Sv_1);

      } else {
        tnext[0] = infinity;
      }

    // next change in Iv
    } else if(mu == 1){

      Ih_0 = Ih_1;

      if(!Ih_trace.empty()){

        Ih_1 = *Ih_trace.begin();
        Ih_trace.erase(Ih_trace.begin());
        tnext[1] = std::get<1>(Ih_1);

      } else {
        tnext[1] = infinity;
      }

    // final interval
    } else {
      break;
    }

    // last tr becomes new tl
    tl = tr;

  }

  return Rcpp::List::create(
    Rcpp::Named("intensity") = intensity,
    Rcpp::Named("left_lims") = left_lims,
    Rcpp::Named("sampled_bites") = sampled_bites
  );
};

/*** R
parameters <- c("a"=0.3,"c"=0.15)
t0 <- 0
dt <- 5
NH <- 100
Sv_mat <- matrix(data = c(
  50,0,
  49,2.3,
  51,4.3
),ncol = 2,byrow = T)
Ih_mat <- matrix(data = c(
  9,0,
  10,3.23
),ncol = 2,byrow = T
)
test_H2M(Sv_mat = Sv_mat,Ih_mat = Ih_mat,NH = NH,t0 = t0,dt = dt,parameters = parameters)
*/
