// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include "qbld.h"
//#include "RcppDist.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo,RcppDist)]]


// This function samples the latent weight W
//--------------------------------------------------------------------------
// Output
// w        : matrix of size (m x n) drawn from GIG distribution
// Input
// z         : Latent response variable, matrix of size (m x n)
// x         : covariates including a column of ones, size (k x m x n)
// s         : random-effects covariates, size (k x m x n)
// beta      : fixed-effects parameter, size (m x 1)
// alpha     : random-effects parameter, size (l x n)
// tau2      : paramter from normal-exponential mixture of Laplace dist.
// theta     : paramter from normal-exponential mixture of Laplace dist.
// lambda    : GIG parameter
//-------------------------------------------------------------------------

/////////// W SAMPLER - start - equation (7) //////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
int sampleW(arma::mat*z, arma::cube*X, arma::cube*S, arma::vec beta, arma::cube*alpha, double tau2, double theta, double lambda, int k, int m, int n,arma::mat*w,int sim)
{
  
  double tilde_eta = (pow(theta,2))/(tau2) + 2;
  double tilde_lambda = 0;
  //arma::mat w(m,n,arma::fill::zeros);
  
  for(int i=0; i<n; i++)
  {
    for(int j=0; j<m; j++)
    {
      tilde_lambda = (pow(((*z)(j,i) - ((((*X).slice(i)).col(j)).t())*beta - ((((*S).slice(i)).col(j)).t())*(((*alpha).slice(sim)).col(i))),2)/tau2).eval()(0,0);
      if(tilde_lambda < 0.00000001) tilde_lambda = 0.00000001;
      
      (*w)(j,i) = rgig(1,lambda,tilde_eta,tilde_lambda)(0);
    }
  }
  return(0);
}

/////////// W SAMPLER - end - equation (7) //////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

