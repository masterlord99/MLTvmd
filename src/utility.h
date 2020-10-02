#ifndef utility_vmd_h
#define utility_vmd_h

#include <Rcpp.h>
using namespace Rcpp;

double tmp_dot_prod(NumericVector& x, NumericVector& y);
Rcomplex complex_dot_prod(ComplexVector& x, ComplexVector& y);
Rcomplex complex_rowsum(const ComplexVector x);
ComplexVector complex_rowsum_matrix(const ComplexMatrix x);
NumericVector complex_abs(ComplexVector x);
NumericVector get_real(ComplexVector& x);
NumericVector get_im(ComplexVector& x);
ComplexVector assign_im(ComplexVector& x,const NumericVector& new_im);
ComplexVector assign_re(ComplexVector& x,const NumericVector& new_re);
void add_re(ComplexVector& x,const NumericVector& new_re);
void sub_re(ComplexVector& x,const NumericVector& new_re);
void mul_re(ComplexVector& x,const NumericVector& new_re);
void div_re(ComplexVector& x,const NumericVector& new_re);
void div_vector(ComplexVector& x,const NumericVector& new_);
void div_scalar(ComplexVector& x,double scal);
void mul_scalar(ComplexVector& x,double scal);
void add_im(ComplexVector& x,const NumericVector& new_im);
NumericVector fftshift_num(NumericVector x,bool inverse=false);
ComplexVector fftshift_com(ComplexVector x,bool inverse=false);


#endif