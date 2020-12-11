#include <vector>
#include <algorithm>
#include <functional>
#include <iostream>

#include <Rcpp.h>

// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp14)]]

#include <boost/icl/interval_map.hpp>


// struct elem {
//   double t;
//   int k;
// };
//
// // [[Rcpp::export]]
// void test(){
//   std::vector<elem> vec;
//   vec.emplace_back(
//     elem{
//       R::runif(0.,10.),
//       0
//     }
//   );
//   vec.emplace_back(
//     elem{
//       R::runif(0.,10.),
//       1
//     }
//   );
//   vec.emplace_back(
//     elem{
//       R::runif(0.,10.),
//       2
//     }
//   );
//   vec.emplace_back(
//     elem{
//       R::runif(0.,10.),
//       3
//     }
//   );
//
//   std::sort(vec.begin(),
//             vec.end(),
//             [](const elem& a, const elem& b) {return a.t < b.t; });
//
//   for(int i=0; i<vec.size(); i++){
//     Rcpp::Rcout << "elem i: " << i << ", t: " << vec[i].t << ", k: " << vec[i].k << " --- \n";
//   }
// }

// [[Rcpp::export]]
void test_icl(const std::vector<double>& input, const std::vector<double>& query){

  // safety
  auto min = std::min_element(input.begin(),input.end());
  auto max = std::max_element(input.begin(),input.end());
  if(query[0] < *(min) || query[1] > *(max)){
    Rcpp::stop("invalid query bounds");
  }

  // the trace! (I_V or I_H)
  boost::icl::interval_map<double, int> mymap;
  for(int i=0; i<input.size()/2; i++){
    mymap += std::make_pair(
      boost::icl::continuous_interval<double>(input[i],input[i+1]),
      1
    );
  };
  Rcpp::Rcout << " --- the map has this many intervals:  " << boost::icl::interval_count(mymap) << " --- \n --- ";
  Rcpp::Rcout <<  mymap << " --- \n";

  // when I was at risk (S_V or S_H)
  boost::icl::continuous_interval<double> query_int = boost::icl::continuous_interval<double>::right_open(query[0],query[1]);
  mymap &= query_int;
  Rcpp::Rcout << "\n --- mymap after restricting to interval has this many intervals: " <<  boost::icl::interval_count(mymap) << " --- \n";
  Rcpp::Rcout << " --- mymap after restricting to interval: " <<  mymap << " --- \n\n";

  // sample a time!
  for(auto it = mymap.begin(); it != mymap.end(); it++){
    double t0 = boost::icl::lower(it->first);
    double t1 = boost::icl::upper(it->first);
    int val = it->second;
    Rcpp::Rcout << "iterating:: dt " <<  t1-t0 << ", val: " << val << " --- \n";
  }
  
  mymap.clear();
  Rcpp::Rcout << "\n --- mymap after clearing has this many intervals: " <<  boost::icl::interval_count(mymap) << ", print the map: " <<  mymap << " --- \n";
}

/*** R
xx <- unlist(replicate(n = 10,expr = {sort(rexp(n = 2))},simplify = FALSE))

yy <- c(0.5,1.2)
test_icl(input = xx,query = yy)
*/
