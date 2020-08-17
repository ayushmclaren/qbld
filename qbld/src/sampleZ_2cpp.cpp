// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include "qbld.h"
#include "RcppDist.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo,RcppDist)]]

// This function samples the latent variable Z marginally
//--------------------------------------------------------------------------
// Output
// z        : matrix of size (m x n) latent variables
// Input
// y         : observed response variable, matrix of size (m x n)
// x         : covariates including a column of ones, size (k x m x n)
// s         : random-effects covariates, size (l x m x n)
// beta      : fixed-effects parameter, size (k x 1)
// varphi2   : scalar, variance of random effects
// w         : matrix, size (m x n)
// tau2      : parameter from normal-exponential mixture of Laplace dist.
// theta     : parameter from normal-exponential mixture of Laplace dist.
//--------------------------------------------------------------------------

/////////// Z SAMPLER - start - equation (10) //////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
int sampleZ_2(arma::mat*y, arma::cube*X, arma::vec beta, arma::cube*S, arma::cube* alpha ,double theta, arma::mat*w, double varphi2, double tau2, int m, int n,arma::mat*z,int sim)
{
  //arma::mat z(m,n,arma::fill::zeros);
  double mean_comp=0,sd_comp=0;
  
  for(int i=0; i<n; i++)
  {
    for(int j=0; j<m; j++)
    {
      mean_comp = (((((*X).slice(i)).col(j)).t())*beta + ((((*S).slice(i)).col(j)).t())*(((*alpha).slice(sim)).col(i)) + theta*(*w)(j,i)).eval()(0,0);
      sd_comp  = sqrt(tau2*(*w)(j,i));
      
      if((*y)(j,i) == 0)
        (*z)(j,i) = r_truncnorm(mean_comp, sd_comp, R_NegInf,0);
      
      if((*y)(j,i) == 1)
        (*z)(j,i) = r_truncnorm(mean_comp, sd_comp,0,R_PosInf);
    }
  }
  return(0);
}
/////////// Z SAMPLER - end - equation (10) //////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
