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
