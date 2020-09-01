// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

//' @useDynLib qbld
//' @importFrom Rcpp sourceCpp
// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include "qbld.h"
//#include "RcppDist.h"

// [[Rcpp::depends(RcppArmadillo)]]

////////////////////////////////////
///// RALD function!

// [[Rcpp::export]]
arma::vec raldmix(double n,double mu,double sigma,double p)
{
  if(n <= 0 || std::floor(n) !=n) 
    Rcpp::stop("Sample size must be a positive integer value.");
  if(sigma <= 0) 
    Rcpp::stop("sigma (scale parameter) must be a positive number.");
  if(p >= 1 || p <= 0) 
    Rcpp::stop("p must be a real number in (0,1).");
  if(fabs(mu) ==R_PosInf) 
    Rcpp::stop("mu (location parameter) must be a finite real number.");
  
  double theta = (1-2*p)/(p*(1-p));
  double tau = sqrt(2/(p*(1-p)));
  
  arma::vec z = Rcpp::rexp(n,1);
  arma::vec u = arma::randn<arma::vec>(n);
  
  arma::vec r = sigma*(theta*z + tau*((arma::sqrt(z))%u));
  r = r + mu;
  return(r);
}
