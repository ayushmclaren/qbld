// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include "qbld.h"
//#include "RcppDist.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo,RcppDist)]]


// This function samples the variance of random-effects, varphi2
//--------------------------------------------------------------------------
// Output
// varphi2   : a scalar quantity
// Input
// alpha     : random-effects parameter, size (l x n)
// c1        : IG prior hyperparameter
// d1        : IG prior hyperparameter
//--------------------------------------------------------------------------

double sampleVarphi2(arma::mat alpha, double c1, double d1, int l, int n)
{
  double sum = 0;
  
  sum = arma::accu(arma::square(alpha));
  
  return(1/R::rgamma((n*l+c1)/2,(sum + d1)/2));
}

//--------------------------------------------------------------------------