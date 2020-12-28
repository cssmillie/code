#include <Rcpp.h>
#include <omp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector dnormpar2(NumericVector x, NumericVector mu, NumericVector sig){

  double c = 1/sqrt(2*PI);  
  int n = x.size();
  int muSize = mu.size();
  int sigSize = sig.size();
  NumericVector ret(n);
  double x0,s0;

  #pragma omp parallel for if(n> 50000) private(x0,s0)
  for(int i=0; i<n; ++i){
    s0 = sig[i % sigSize];
    x0 = x[i]-mu[i % muSize];
    ret[i] = exp(-x0*x0/(2*s0*s0))*c/s0;
  }

  return ret;

}
