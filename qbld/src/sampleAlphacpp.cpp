// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include "qbld.h"
//#include "RcppDist.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]


// This function samples the random effects parameter alpha marginally
//--------------------------------------------------------------------------
// Output
// alpha     : Gibbs draw of alpha, matrix of size (l x n)
// Input
// z         : Latent response variable, matrix of size (m x n)
// x         : covariates including a column of ones, size (k x m x n)
// s         : random-effects covariates, size (l x m x n)
// beta      : fixed-effects parameter, size (k x 1)
// w         : matrix, size (m x n)
// tau2      : parameter from normal-exponential mixture of Laplace dist.
// theta     : parameter from normal-exponential mixture of Laplace dist.
// varphi2   : scalar quantity
//--------------------------------------------------------------------------


/////////// Alpha SAMPLER - start - equation (6) //////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

int sampleAlphafast(arma::mat*z, arma::cube*X, arma::cube*S, arma::vec beta, arma::mat*w, double tau2, double theta, double varphi2, int l, int m, int n,arma::cube*alpha,int sim)
{
  arma::mat invD1(m,m,arma::fill::zeros);
  arma::mat Im(m,m,arma::fill::eye);
  arma::mat phii(m,m);
  arma::mat phii_t(m,m); //saving transpose saves time apparently in my tests, hence!
  arma::vec u(l,arma::fill::randn); //randn is fast and does quick N(0,1) samples
  u = sqrt(varphi2)*u; 
  arma::vec delta(m,arma::fill::randn);
  //arma::vec v(m);
  arma::vec alph(m);
  arma::vec w_i(m);
  
  //u = Rcpp::rnorm(l,0,varphi);
  //delta = Rcpp::rnorm(m,0,1); 
  //double varphi2 = pow(varphi,2);
  
  for(int i=0; i<n; i++)
  {
    invD1.diag() = 1/(sqrt(tau2*(*w).col(i)));
    phii = invD1*(((*S).slice(i)).t());
    //phii_t = phii.t();
    //v = phii*u + delta;
    alph = invD1*((*z).col(i) - (((*X).slice(i)).t())*beta - theta*(*w).col(i)) - phii*u - delta; //-v done in this step
    w_i = solve(varphi2*(phii*phii.t()) + Im, alph,arma::solve_opts::fast); //solve_optss_fast disables some  unnecessary checks and hence is faster
    ((*alpha).slice(sim)).col(i) = u + varphi2*(phii.t())*w_i; 
  }
  return(0);
}


int sampleAlpha(arma::mat*z, arma::cube*X, arma::cube*S, arma::vec beta, arma::mat*w, double tau2, double theta, double varphi2, int l, int m, int n,arma::cube*alpha,int sim)
{
  
  //arma::mat alpha_out(l,n,arma::fill::zeros);
  arma::mat invDvarphi2 = (arma::eye(l,l))/varphi2;
  arma::mat invD1(l,m,arma::fill::zeros);
  arma::vec atilde(l,arma::fill::zeros);
  arma::mat Atilde(l,l,arma::fill::zeros);
  
  for(int i=0; i<n; i++)
  {
    invD1 = ((*S).slice(i))*(arma::diagmat(tau2*((*w).col(i))).i());
    //diagmat(1/(tau2*((*w).col(i))));
    
    //  Rcpp::Rcout << "nsim_alpha:" << i << "\n";
    //  Rcpp::Rcout << "varphi2" << varphi2 << "\n";
    //  Rcpp::Rcout << "invD1:" << invD1 << "\n";
    //  Rcpp::Rcout << "invDvarphi2:" << invDvarphi2 << "\n";
    //  Rcpp::Rcout << "*S:" << (*S).slice(i) << "\n";
    
    Atilde = ((invD1*(((*S).slice(i)).t())) + invDvarphi2).i();
    //arma::inv_sympd((((*S).slice(i))*invD1*(((*S).slice(i)).t())) + invDvarphi2);
    // Rcpp::Rcout << "Atilde:" << Atilde << "\n";
    //Atilde = arma::inv_sympd(Atilde);
    atilde =  Atilde*(invD1*((*z).col(i) - (((*X).slice(i)).t())*beta - theta*(*w).col(i)));      
    
    ((*alpha).slice(sim)).col(i) = arma::mvnrnd(atilde,Atilde,1);
  }
  return(0);
}

/////////// Alpha SAMPLER - end - equation (6) //////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////


