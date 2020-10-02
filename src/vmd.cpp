#include "fft.h"

#include <Rcpp.h>
using namespace Rcpp;

#include "utility.h"



// [[Rcpp::export]]
List vmd(
    NumericVector signal,
    double tau = 0,
    int alpha = 2000,
    int K = 3,
    double tol = 0.00001,
    int N = 500,
    int init = 0,
    bool DC = true
) {
  
  int lenOrg = signal.size();
  double fs = 1.0/signal.size();
  // std::cout << fs;
  int hw = lenOrg/2;
  NumericVector lhs = head(signal,hw);
  lhs = rev(lhs);
  NumericVector rhs = tail(signal,lenOrg-  hw);
  rhs = rev(rhs);
  
  NumericVector signalMir = rep(NA_REAL,lhs.size() + signal.size() + rhs.size());
  signalMir[Range(0,lhs.size()-1)] = lhs;
  signalMir[Range(lhs.size(),lhs.size() + signal.size()-1)] = signal;
  signalMir[Range(lhs.size() + signal.size(),signalMir.size())] = rhs;
  
  int lenMir = signalMir.size();
  IntegerVector help = Range(1,lenMir);
  NumericVector t = as<NumericVector>(help);
  t = t/lenMir;
  
  NumericVector freqs = t - 0.5 - (1.0/lenMir);
  IntegerVector Alpha = rep(alpha,K);
  
  ComplexVector f_hat = fft_arma(signalMir);
  f_hat = fftshift_com(f_hat);
  ComplexVector f_hat_plus = clone(f_hat);
  
  Rcomplex helper;
  int i=0;
  while(i <floor(lenMir/2)){
    f_hat_plus[i] = helper;
    i++;
  }
  
  NumericMatrix omega_plus = NumericMatrix(N,K);
  
  if(init==1){
    IntegerVector tmp_k_ = Range(0,K-1);
    NumericVector tmp_k = as<NumericVector>(tmp_k_);
    tmp_k = (0.5/K) * tmp_k;
    omega_plus(0,_)= tmp_k;
  }else if(init== 2){
    //
  }
  
  if(DC)omega_plus(0,0)=0;
  
  ComplexMatrix lambda_hat = ComplexMatrix(N,lenMir);
  ComplexMatrix u_hat_plus = ComplexMatrix(lenMir,K);
  
  IntegerVector ix = Range(floor(lenMir/2),lenMir-1);
  NumericVector freqs_ix = freqs[ix];
  double uDiff = 1e10;
  int n = 0;
  ComplexVector sum_uk(lenMir);
  
  ComplexMatrix u_hat_plus_0 = ComplexMatrix(lenMir,K);
  
  while((uDiff > tol) & (n < (N-1))){
    
    ComplexMatrix u_hat_plus_1 = ComplexMatrix(lenMir,K);
    
    int k=0;
    while(k < K){
      
      //accumulator
      if(k==1){
        sum_uk = u_hat_plus_0(_,K) + sum_uk - u_hat_plus_0(_,k);
      }else{
        sum_uk = u_hat_plus_1(_,k-1) + sum_uk - u_hat_plus_0(_,k);
      }
      
      // accumulator finished
      
      // mode spectrum
      ComplexVector target_col = u_hat_plus_1(_,k);
      
      target_col = f_hat_plus - sum_uk;
      ComplexVector tmp_lambda_hat = lambda_hat(n,_);
      div_scalar(tmp_lambda_hat,2.0);
      target_col = target_col - tmp_lambda_hat;
      
      NumericVector omega_help = freqs - omega_plus(n,k);
      omega_help = pow(omega_help,2);
      omega_help = omega_help * Alpha[k] + 1;
      
      div_vector(target_col,omega_help);
      u_hat_plus_1(_,k) = target_col;

      
      // mode spectrum finished
      
      // center frequencies
      
      if(!DC | (k > 1)){
        
        ComplexVector tmp_u_hat_plus_1 = u_hat_plus_1(_,k);
        tmp_u_hat_plus_1 = tmp_u_hat_plus_1[ix];
        NumericVector tmp_hat_1 = complex_abs(tmp_u_hat_plus_1);
        tmp_hat_1 = pow(tmp_hat_1,2);
        double bot_side = sum(tmp_hat_1);
        double top_side = tmp_dot_prod(freqs_ix,tmp_hat_1);
        omega_plus(n+1,k) = (top_side/bot_side);
        
      }
      
      k++;
    }
    
    // dual ascent
    
    ComplexMatrix::Row target_row = lambda_hat(n+1,_);
    
    ComplexVector row_sums_u_hat_plus_1 = complex_rowsum_matrix(u_hat_plus_1);
    mul_scalar(row_sums_u_hat_plus_1,tau);
    
    target_row = lambda_hat(n,_); - row_sums_u_hat_plus_1 - f_hat_plus;
    
    n++;
    
    // converged yet?
    
    Rcomplex uDiff;
    for(int i=0;i<K;i++){
      ComplexVector a = u_hat_plus_1(_,i) - u_hat_plus_0(_,i);
      Rcomplex tmp_diff = complex_dot_prod(a,a);
      uDiff.r = uDiff.r + tmp_diff.r/lenMir;
      uDiff.i = uDiff.i + tmp_diff.i/lenMir;
    }
    double u_diff = sqrt(uDiff.r*uDiff.r + uDiff.i*uDiff.i);
    
    if(R_IsNA(u_diff)){
      std::cout << "Divergence! Check params!\n";
      return(0);
    }
    
    // reassign holders;
    
    u_hat_plus_0 = u_hat_plus_1;
    
  }
  
  // post-process and cleanup
  
  NumericMatrix omega = omega_plus(Range(0,n-1),_);
  ComplexMatrix u_hat = ComplexMatrix(lenMir,K);
  // here done for today
  
  
  return 0;
}

