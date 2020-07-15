// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include "qbld.h"
//#include "RcppDist.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo,RcppDist)]]


// This function draws the fixed effects parameter beta marginally (unblocked)
//--------------------------------------------------------------------------
// Output
// beta      : Gibbs draw of beta, a column vector (k x 1)
// Input
// z         : Latent response variable, matrix of size (m x n)
// x         : covariates including a column of ones, size (k x m x n)
// s         : random-effects covariates, size (l x m x n)
// alpha     : random-effects parameter, matrix of size (l x n)
// w         : matrix, size (m x n)
// varphi2   : scalar quantity
// tau2      : parameter from normal-exponential mixture of Laplace dist.
// theta     : parameter from normal-exponential mixture of Laplace dist.
// invB0     : inverse of the prior covariance matrix
// invB0b0   : prior precision times prior mean
//--------------------------------------------------------------------------

arma::vec sampleBeta_2(arma::mat*z, arma::cube*X, arma::cube*S, arma::mat*w, arma::mat alpha ,double varphi2, double tau2, double theta, arma::mat*invB0, arma::mat*invB0b0, int k, int m, int n)
{
  arma::mat sumvar(k,k,arma::fill::zeros);
  arma::vec summean(k,arma::fill::zeros);
  arma::mat vari(k,k,arma::fill::zeros);
  arma::vec meani(k,arma::fill::zeros);
  arma::mat inv_phi(m,m,arma::fill::zeros);
  
  for(i in 1:n)
  {
    inv_phi.diag() = 1/(tau2*((*w).col(i)));
    
    vari = ((*X).slice(i))*inv_phi*(((*X).slice(i)).t());
    meani = ((*X).slice(i))*inv_phi*((*z).col(i) - theta*((*w).col(i)) - (((*S).slice(i)).t())*alpha.col(i));
    
    sumvar  = sumvar + vari;
    summean = summean + meani;
  }
  
  vari  = arma::inv_sympd((*invB0) + sumvar);
  meani = vari*((*invB0b0) + summean);
    
  return(arma::mvnrnd(meani,vari,1));
}
