/* --------------------------------------------------------------------------------
#
#   disaggregated malaria transmission model ABM
#   test humans with constant value for Iv
#   Sean L. Wu (slwu89@berkeley.edu)
#   December 2020
#
-------------------------------------------------------------------------------- */

#include <algorithm>

#include <Rcpp.h>

// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp14)]]


#include "abm/human.hpp"

// [[Rcpp::export]]
Rcpp::List test_humans(
  const int SH, const int IH, const double IV,
  const Rcpp::NumericVector& parameters,
  const double dt, const double tmax
){

  human_pop_uptr hpop = std::make_unique<human_pop>();

  double a = parameters["a"];
  double b = parameters["b"];

  hpop->N = SH + IH;
  hpop->r = parameters["r"];
  hpop->LEP = parameters["LEP"];

  int i;
  for(i=0; i<SH; i++){
    hpop->pop.emplace_back(make_human(hpop,'S'));
  }
  for(i=0; i<IH; i++){
    hpop->pop.emplace_back(make_human(hpop,'I'));
  }

  double clock{0.0};
  int tsteps = tmax/dt;
  int j{0};

  Rcpp::Rcout << "--- begin simulation ---\n";
  while(clock < tmax){

    run_humanpop(hpop,clock,dt);

    // compute bloodmeals (M2H bites)
    for(auto& it : hpop->pop){
      if(it->snow == 'S'){
        assert(it->tnow >= clock);
        double myrate = a*b*IV*(1./hpop->N) * ((clock+dt) - it->tnow);
        double myprob = 1. - exp(-myrate);
        if(R::runif(0.,1.) < myprob){
          double btime = R::runif(it->tnow,(clock+dt));
          push_M2H_bites(it,btime);
        }
      }
    }

    hpop->Ih_trace.clear();

    j++;
    if(j % 100 == 0){
      Rcpp::Rcout << " completed TWICE step " << j << " of " << tsteps << "\n";
    }
    clock += dt;

  };

  // return data
  std::vector<int> ids;
  std::vector<char> states;
  std::vector<double> times;
  int id{0};

  for(auto& it : hpop->pop){
    int njump = it->thist.size();
    states.insert(states.end(),it->shist.begin(),it->shist.end());
    times.insert(times.end(),it->thist.begin(),it->thist.end());
    ids.insert(ids.end(),njump,id);
    id++;
  }

  // process the actual history
  std::sort(
    hpop->hist.begin(),
    hpop->hist.end(),
    [](const hist_elem& a, const hist_elem& b){
      return a.t < b.t;
    }
   );

   int nn =hpop->hist.size();
   Rcpp::NumericMatrix outagg(nn,4);
   for(int i=0; i<nn; i++){
     outagg.at(i,0) = hpop->hist[i].t;
     outagg.at(i,1) = hpop->hist[i].dS;
     outagg.at(i,2) = hpop->hist[i].dE;
     outagg.at(i,3) = hpop->hist[i].dI;
   };

  return Rcpp::List::create(
    Rcpp::Named("ids") = Rcpp::wrap(ids),
    Rcpp::Named("times") = Rcpp::wrap(times),
    Rcpp::Named("states") = Rcpp::wrap(states),
    Rcpp::Named("agghist") =outagg
  );
};
