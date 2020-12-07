/* --------------------------------------------------------------------------------
#
#   disaggregated malaria transmission model ABM
#   mosquito
#   Sean L. Wu (slwu89@berkeley.edu)
#   November 2020
#
-------------------------------------------------------------------------------- */

#include <list>

#include <Rcpp.h>


// a single skeeter (syringe with wings)
typedef struct mosquito {
  double tnow;
  double tnext;
  double tdie;
  char state;
  char statenext;
} mosquito;

// using

typedef struct mosypop_str {

  std::list<mosquito> mosypop;

  // history tracking
  std::vector<double> t_hist;
  std::vector<double> S_hist;
  std::vector<double> E_hist;
  std::vector<double> I_hist;

  double tnow;

  // bites when H->M transmission occured
  queue H2M_bites;

  // when the resulting infection in M will manifest
  queue M_inf;

  // parameters
  double g;
  double EIP;
  double lambda;

  // elements for MNRM of internal mosy dynamics
  std::array<double,5> delta_t{infinity};
  std::array<double,4> Pk{0.};
  std::array<double,4> Tk{0.};
  std::array<double,4> ak{0.};

  // integrated susceptible (Sv) biting over [t0,t0+dt)
  // integrated infectious (Iv) biting over [t0,t0+dt)
  queue_trace Sv_trace;
  queue_trace Iv_trace;

  mosypop_str(){Rcpp::Rcout << "'mosypop_str' ctor called at " << this << std::endl;};
  ~mosypop_str(){Rcpp::Rcout << "'mosypop_str' dtor called at " << this << std::endl;};

} mosypop_str;
