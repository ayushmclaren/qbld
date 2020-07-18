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

arma::mat sampleAlpha(arma::mat*z, arma::cube*X, arma::cube*S, arma::vec beta, arma::mat*w, double tau2, double theta, double varphi2, int l, int m, int n)
{
  
  arma::mat alpha_out(l,n,arma::fill::zeros);
  arma::mat invDvarphi2 = (arma::eye(l,l))/varphi2;
  arma::mat invD1(m,m,arma::fill::zeros);
  arma::vec atilde(l,arma::fill::zeros);
  arma::mat Atilde(l,l,arma::fill::zeros);
  
  for(int i=0; i<n; i++)
  {
    invD1.diag() = 1/(tau2*((*w).col(i)));
    
    Atilde = inv_sympd((((*S).slice(i))*invD1*(((*S).slice(i)).t())) + invDvarphi2);
    atilde =  Atilde*(((*S).slice(i))*invD1*((*z).col(i) - (((*X).slice(i)).t())*beta - theta*(*w).col(i)));      
    
    alpha_out.col(i) = arma::mvnrnd(atilde,Atilde,1);
  }
  return(alpha_out);
}



// [[Rcpp::export]]

arma::vec ok2(int tt)
{
  int l  = 5;
  double varphi = 2;
  
  return(Rcpp::rnorm(l,0,varphi));
}
