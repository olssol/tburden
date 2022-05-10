#include <Rcpp.h>
#include <Rmath.h>

using namespace Rcpp;

// [[Rcpp::plugins("cpp11")]]

// [[Rcpp::init]]
void my_package_init(DllInfo *dll) {
  // initialization code here
  R_useDynamicSymbols(dll, TRUE);
}


// Test Rcpp function
//
//
// @param test test parameter
//
// @export
// [[Rcpp::export]]
double crtTest(double test) {
  return pow(test,2);
}

// Get utility value by survival function
//
// @export
//
// [[Rcpp::export]]
double c_uti_surv(NumericMatrix surv_f) {

  int    i;
  double l, p, rst = 0;

  for (i = 0; i < surv_f.nrow() - 1; i++) {
    l   =  surv_f(i + 1, 0) - surv_f(i, 0);
    p   =  1 - surv_f(i, 1);
    rst += l * p;
  }


  // return
  return(rst);
}


// Get Pseudo Overall Response
//
// @export
//
// [[Rcpp::export]]
NumericVector c_pseudo_response(NumericVector pchg, double thresh_cr, double thresh_pd) {
  int    i, n = pchg.length();
  double cp, cm, tmp;

  // t_or, pchg_or, t_pd, pchg_or, t_min, pchg_min
  NumericVector rst(6, 999.0);

  for (i = 0; i < n; i++) {

    cp = pchg[i];
    cm = rst[5];

    //Response
    if (cp < thresh_cr & 999 == rst[0]) {
      rst[0] = i + 1;
      rst[1] = cp;
    }

    //Progression
    tmp = cm + (1 + cm) * thresh_pd;
    if (cp > tmp & 999 == rst[2]) {
      rst[2] = i + 1;
      rst[3] = cp;
    }

    //Minimum
    if (cp < cm) {
      rst[4] = i + 1;
      rst[5] = cp;
    }
  }

  //Return
  for (i = 0; i < rst.length(); i++) {
    if (999 == rst[i]) {
      rst[i] = NA_REAL;
    }
  }

  return(rst);
}
