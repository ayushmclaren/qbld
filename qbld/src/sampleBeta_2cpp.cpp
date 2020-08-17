// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include "qbld.h"
//#include "RcppDist.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]


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

/////////// BETA SAMPLER - start - equation (9) //////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
int sampleBeta_2(arma::mat*z, arma::cube*X, arma::cube*S, arma::mat*w, arma::cube *alpha ,double varphi2, double tau2, double theta, arma::mat*invB0, arma::mat*invB0b0, int k, int m, int n,arma::mat* beta,int sim)
{
  arma::mat sumvar(k,k,arma::fill::zeros);
  arma::vec summean(k,arma::fill::zeros);
  //arma::mat vari(k,k,arma::fill::zeros);
  //arma::vec meani(k,arma::fill::zeros);
  arma::mat inv_phi(k,m,arma::fill::zeros);
  
  for(int i=0; i<n ; i++)
  {
    inv_phi = ((*X).slice(i))*(arma::diagmat(tau2*((*w).col(i))).i());
    
    //vari 
    sumvar  += inv_phi*(((*X).slice(i)).t());
    //meani 
    summean  += inv_phi*((*z).col(i) - theta*((*w).col(i)) - (((*S).slice(i)).t())*((*alpha).slice(sim-1)).col(i));
    
    // sumvar  = sumvar + vari;
    //  summean = summean + meani;
  }
  
  sumvar  = arma::inv_sympd((*invB0) + sumvar);
  summean = sumvar*((*invB0b0) + summean);
  
  (*beta).col(sim) = arma::mvnrnd(summean,sumvar,1);
  return(0);
}
/////////// BETA SAMPLER - end - equation (9) //////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
