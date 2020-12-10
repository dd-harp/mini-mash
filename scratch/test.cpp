#include <vector>
#include <algorithm>
#include <functional>
#include <iostream>

#include <Rcpp.h>

// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp14)]]  

#include <boost/icl/interval_map.hpp>


struct elem {
  double t;
  int k;
};

// [[Rcpp::export]]
void test(){
  std::vector<elem> vec;
  vec.emplace_back(
    elem{
      R::runif(0.,10.),
      0
    }
  );
  vec.emplace_back(
    elem{
      R::runif(0.,10.),
      1
    }
  );
  vec.emplace_back(
    elem{
      R::runif(0.,10.),
      2
    }
  );
  vec.emplace_back(
    elem{
      R::runif(0.,10.),
      3
    }
  );
  
  std::sort(vec.begin(), 
            vec.end(), 
            [](const elem& a, const elem& b) {return a.t < b.t; });
  
  for(int i=0; i<vec.size(); i++){
    Rcpp::Rcout << "elem i: " << i << ", t: " << vec[i].t << ", k: " << vec[i].k << " --- \n";
  }
}

// [[Rcpp::export]]
void test_icl(const Rcpp::NumericVector& input){
  boost::icl::interval_map<double, int> mymap;
  for(int i=0; i<input.size()/2; i++){
    mymap += std::make_pair(
      boost::icl::continuous_interval<double>(input[i],input[i+1]),
      1
    );
  }
  Rcpp::Rcout << "the map: " << mymap << " --- \n";
}