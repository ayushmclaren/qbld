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


/////////// Varphi2 SAMPLER - start - equation (8) //////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//--------------------------------------------------------------------------

double sampleVarphi2(arma::cube* alpha, double c1, double d1, int l, int n,int sim)
{
  double sum = 0;
  
  sum = arma::accu(arma::square((*alpha).slice(sim)));
  
  return(1/R::rgamma((n*l+c1)/2,2/(sum + d1)));
  //return(1/arma::randg<double>(arma::distr_param((n*l+c1)/2,2/(sum + d1))));
}

/////////// Varphi2 SAMPLER - end - equation (8) //////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

//--------------------------------------------------------------------------