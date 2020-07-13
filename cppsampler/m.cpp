// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
//#include "RcppDist.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo,RcppDist)]]

// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//

// [[Rcpp::export]]

arma::mat ok(int tt)
{
  int m = 5;
  int n = 8;
  arma::mat w(m,n,arma::fill::ones);
  arma::mat invD1(m,m,arma::fill::zeros);
  double tau2 = 0.25;
  
  invD1.diag() = 1/(tau2*((w).col(4)));
  
  return(invD1);
  
}