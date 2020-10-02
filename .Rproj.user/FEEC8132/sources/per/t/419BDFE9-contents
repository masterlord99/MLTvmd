#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double tmp_dot_prod(NumericVector& x, NumericVector& y){
  double out;
  for(int i = 0; i < x.size();i++){
    out = out + x[i]*y[i];
  }
  return(out);
}



// [[Rcpp::export]]
Rcomplex complex_dot_prod(ComplexVector& x, ComplexVector& y){
  Rcomplex out;
  ComplexVector new_y = Conj(y);
  for(int i = 0; i < x.size();i++){
    out.r = out.r + x[i].r*new_y[i].r;
    out.i = out.i + x[i].i*new_y[i].i;
  }
  return(out);
}



// [[Rcpp::export]]
Rcomplex complex_rowsum(const ComplexVector x){
  Rcomplex out = x[0];
  for(int i =1;i<x.size();i++){
    out.r = out.r + x[i].r;
    out.i = out.i + x[i].i;
  }
  return(out);
}


// [[Rcpp::export]]
ComplexVector complex_rowsum_matrix(const ComplexMatrix x){
  ComplexVector out(x.rows());
  for(int i=0;i <x.rows();i++){
    out[i] = complex_rowsum(x(i,_));
  }
  return(out);
}


// [[Rcpp::export]]
NumericVector complex_abs(ComplexVector x){
  NumericVector out(x.size());
  for(int i = 0; i < x.size();i++){
    Rcomplex tmp_ = x[i];
    out[i] = (tmp_.i*tmp_.i) + (tmp_.r*tmp_.r);
  }
  return(sqrt(out));
}

// [[Rcpp::export]]
NumericVector get_real(ComplexVector& x){
  NumericVector out(x.size());
  for(int i =0; i < x.size();i++){
    Rcomplex tmp = x[i];
    out[i] = tmp.r;
  }
  return(out);
}

// [[Rcpp::export]]
NumericVector get_im(ComplexVector& x){
  NumericVector out(x.size());
  for(int i =0; i < x.size();i++){
    Rcomplex tmp = x[i];
    out[i] = tmp.i;
  }
  return(out);
}

// [[Rcpp::export]]
ComplexVector assign_im(ComplexVector& x,const NumericVector& new_im){
  for(int i =0; i < x.size();i++){
    Rcomplex& tmp = x[i];
    tmp.i = new_im[i];
  }
  return(x);
}

// [[Rcpp::export]]
ComplexVector assign_re(ComplexVector& x,const NumericVector& new_re){
  for(int i =0; i < x.size();i++){
    Rcomplex& tmp = x[i];
    tmp.r = new_re[i];
  }
  return(x);
}



// [[Rcpp::export]]
void add_re(ComplexVector& x,const NumericVector& new_re){
  for(int i =0; i < x.size();i++){
    Rcomplex& tmp = x[i];
    tmp.r = tmp.r + new_re[i];
  }
}


// [[Rcpp::export]]
void sub_re(ComplexVector& x,const NumericVector& new_re){
  for(int i =0; i < x.size();i++){
    Rcomplex& tmp = x[i];
    tmp.r = tmp.r - new_re[i];
  }
}

// [[Rcpp::export]]
void mul_re(ComplexVector& x,const NumericVector& new_re){
  for(int i =0; i < x.size();i++){
    Rcomplex& tmp = x[i];
    tmp.r = tmp.r * new_re[i];
  }
}

// [[Rcpp::export]]
void div_re(ComplexVector& x,const NumericVector& new_re){
  for(int i =0; i < x.size();i++){
    Rcomplex& tmp = x[i];
    tmp.r = tmp.r / new_re[i];
  }
}


// [[Rcpp::export]]
void div_vector(ComplexVector& x,const NumericVector& new_){
  for(int i =0; i < x.size();i++){
    Rcomplex& tmp = x[i];
    tmp.r = tmp.r / new_[i];
    tmp.i = tmp.i / new_[i];
  }
}


// [[Rcpp::export]]
void div_scalar(ComplexVector& x,double scal){
  for(int i =0; i < x.size();i++){
    Rcomplex& tmp = x[i];
    tmp.r = tmp.r / scal;
    tmp.i = tmp.i / scal;
  }
}

// [[Rcpp::export]]
void mul_scalar(ComplexVector& x,double scal){
  for(int i =0; i < x.size();i++){
    Rcomplex& tmp = x[i];
    tmp.r = tmp.r * scal;
    tmp.i = tmp.i * scal;
  }
}


// [[Rcpp::export]]
void add_im(ComplexVector& x,const NumericVector& new_im){
  for(int i =0; i < x.size();i++){
    Rcomplex& tmp = x[i];
    tmp.i = tmp.i + new_im[i];
  }
}

// [[Rcpp::export]]
NumericVector fftshift_num(NumericVector x,bool inverse=false){
  
  double len = x.size();
  int hw;
  if(inverse){
    hw = ceil(len/2);
  }else{
    hw = floor(len/2);
  }
  NumericVector out = rep(NA_REAL,len);
  NumericVector tmp1 = x[Range(hw,len-1)];
  NumericVector tmp2 = x[Range(0,hw-1)];
  out[Range(0,len -hw -1)] = tmp1;
  // std::cout << tmp1<< "\n";
  out[Range(len -hw ,len-1)] = tmp2;
  return(out);
}


// [[Rcpp::export]]
ComplexVector fftshift_com(ComplexVector x,bool inverse=false){
  double len = x.size();
  int hw;
  if(inverse){
    hw = ceil(len/2);
  }else{
    hw = floor(len/2);
  }
  ComplexVector out =ComplexVector(len);
  ComplexVector tmp1 = x[Range(hw,len-1)];
  ComplexVector tmp2 = x[Range(0,hw-1)];
  out[Range(0,len -hw -1)] = tmp1;
  // std::cout << tmp1<< "\n";
  out[Range(len -hw ,len-1)] = tmp2;
  return(out);
}