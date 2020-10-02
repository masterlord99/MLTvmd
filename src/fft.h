#ifndef fourier_h
#define fouirer_h

#include <RcppArmadillo.h>
using namespace Rcpp;

ComplexVector fft_arma(const NumericVector & y_);
NumericVector ifft_arma(const ComplexVector & y_);

#endif

