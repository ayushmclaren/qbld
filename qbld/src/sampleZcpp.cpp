// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include "qbld.h"
//#include "RcppDist.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]


// This function samples the latent variable Z marginally
//--------------------------------------------------------------------------
// Output
// z        : matrix of size (m x n) latent variables
// Input
// zPrev     : Previous draw in the MCMC
// y         : observed response variable, matrix of size (m x n)
// x         : covariates including a column of ones, size (k x m x n)
// s         : random-effects covariates, size (l x m x n)
// beta      : fixed-effects parameter, size (k x 1)
// varphi2   : scalar, variance of random effects
// w         : matrix, size (m x n)
// tau2      : parameter from normal-exponential mixture of Laplace dist.
// theta     : parameter from normal-exponential mixture of Laplace dist.
//--------------------------------------------------------------------------

/*
 double erf_my(double x, bool inverse) 
 {
 if (inverse) {
 double ans = R::qnorm((x+1)/2,0,1,1,0) / sqrt(2);
 
 if(x < -1) ans = NA_REAL;
 if(x > 1) ans = NA_REAL;
 if(x == -1) ans = R_NegInf;
 if(x == 1) ans = R_PosInf;
 
 return(ans);
 } 
 else {
 return(2*(R::pnorm(x*sqrt(2),0,1,1,0)) - 1); 
 }
 }
 */

/////////// Z SAMPLER - start - equation (5) //////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

double erfc_my(double x, bool inverse) 
{
  if (inverse) {
    double ans = R::qnorm(x/2,0,1,false,false)/sqrt(2);
    
    if(x <  0) ans = NA_REAL;
    if(x > 2) ans = NA_REAL;
    if(x == 2) ans = R_NegInf;
    if(x == 0) ans = R_PosInf;
    
    return(ans);
  } 
  else {
    return(2*(R::pnorm(x*sqrt(2),0,1,false,false)));
  }
}


double normcdf(double x)
{
  return(0.5*erfc_my(-x/sqrt(2),FALSE));
}

double norminv(double x)
{
  return(-sqrt(2)*erfc_my(2*x,TRUE));
}

double limits(double y,int which)
{
  if(which == 0)  //lowerlimits
  {
    if(y == 0) return(R_NegInf);
    else return(0);
  }
  
  else if(which==1) //upperlimits
  {
    if(y == 0) return(0);
    else return(R_PosInf);
  }
  return(0);
}

arma::vec rtruncnorm_gwk(arma::vec z0,arma::vec*mu,arma::mat*sigma,arma::mat*y,int m,int idx)
{
  double newLL=0, newUL=0, condMu=0, sq_condVar=0, mu_i=0, sigma11=0;
  arma::mat inv_sigma22(m-1,m-1,arma::fill::zeros);
  arma::vec u = arma::randu(m); 
  
  for (int i=0; i<m; i++)
  {
    mu_i = (*mu)(i);
    
    arma::vec mu_less_i = *mu;
    mu_less_i.shed_row(i);        
    
    sigma11 = (*sigma)(i,i);
    arma::rowvec sigma12 = (*sigma).row(i);
    sigma12.shed_col(i);
    
    arma::mat sigma22 = *sigma;
    sigma22.shed_row(i);
    sigma22.shed_col(i);
    inv_sigma22 = sigma22.i();
    
    sq_condVar = sqrt((sigma11 - sigma12*inv_sigma22*(sigma12.t())).eval()(0,0)); 
    
    arma::vec z_less_i = z0;
    z_less_i.shed_row(i);        
    
    condMu = (mu_i + sigma12*inv_sigma22*(z_less_i - mu_less_i)).eval()(0,0);
    //(mu_i + sigma12*inv_sigma22*(z_less_i - mu_less_i)).eval()(0,0);
    
    newLL  = (limits((*y)(i,idx),0) - condMu)/sq_condVar;
    newUL= (limits((*y)(i,idx),1) - condMu)/sq_condVar;
    
    z0(i) =  condMu + (sq_condVar)*(norminv(u(i)*(normcdf(newUL) - normcdf(newLL)) + normcdf(newLL)));
    //z0[i] = z[i]
  }
  
  return(z0);
}

int sampleZ(arma::mat*zprev, arma::mat*y, arma::cube*X, arma::vec beta, arma::cube*S, double theta, arma::mat*w, double varphi2, double tau2, int m, int n, arma::mat*z)
{
  
  //arma::mat z(m,n,arma::fill::zeros);
  arma::vec meani(m,arma::fill::zeros);
  arma::mat VarZi(m,m,arma::fill::zeros);
  //arma::mat D1(m,m,arma::fill::zeros);
  
  for(int i=0; i<n; i++)
  {
    //D1.diag() = tau2*((*w).col(i));
    meani = (((*X).slice(i)).t())*beta + theta*((*w).col(i));
    VarZi = (varphi2*((((*S).slice(i)).t())*((*S).slice(i))) + arma::diagmat(tau2*((*w).col(i))));
    
    (*z).col(i) = rtruncnorm_gwk((*zprev).col(i),&meani,&VarZi,y,m,i);
  }
  return(0);
  
}

/////////// Z SAMPLER - end - equation (5) //////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////