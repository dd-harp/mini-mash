#include <vector>
#include <algorithm>
#include <functional>

#include <Rcpp.h>

// [[Rcpp::plugins(cpp14)]]  


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