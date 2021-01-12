#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::NumericMatrix discretise(const Rcpp::NumericMatrix& out, const double dt)
{
  int ncol = out.ncol();
  int events = out.nrow();

  double end = out.at(events-1,0);
  double start = out.at(0,0);

  int len = (int)((end - start) / dt) + 1;
  Rcpp::NumericMatrix x(len,ncol);

  double target{0.};
  int j{0};

  for(int i=0; i<events; i++){
    while(out.at(i,0) >= target){
      if(j>=len){
        Rcpp::Rcout << " --- overrunning matrix --- \n";
        break;
      }
      x.at(j,0) = target;
      for(int k=1; k<ncol; k++){
        x.at(j,k) = out.at(i,k);
      }
      j += 1;
      target += dt;
    }
  }

  Rcpp::CharacterVector outnames = Rcpp::colnames(out);
  Rcpp::colnames(x) = outnames;
  return x;

}
