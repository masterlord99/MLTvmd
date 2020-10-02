// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
ComplexVector fft_arma(const NumericVector & y_) {
  
  const arma::vec & y = as<arma::vec>(wrap(y_));
  arma::cx_vec coef = fft(y);
  
  return(
    as<ComplexVector>(wrap(coef))
  );
}


// [[Rcpp::export]]
NumericVector ifft_arma(const ComplexVector & y_) {
  
  const arma::cx_vec & y = as<arma::cx_vec>(wrap(y_));
  arma::cx_vec coef = ifft(y);
  arma::vec coef_ = arma::real(coef);
  
  return(
    as<NumericVector>(wrap(coef_))
  );
}