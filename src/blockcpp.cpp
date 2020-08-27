//' @useDynLib qbld
//' @importFrom Rcpp sourceCpp
// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
#include "qbld.h"
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

// picks out every jth column from a matrix from start to end
/////////// Combined blocked SAMPLER - start  //////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

arma::mat subset_mat(arma::mat* X, int start, int j) 
{
  arma::uvec IDX = arma::regspace<arma::uvec>(start,  j,  (*X).n_cols-1);
  return (*X).cols(IDX);
}

// y is the output variable, x is fixed, s is random 
// [[Rcpp::export]]
Rcpp::List qbldf(int nsim, double p, arma::mat y, arma::mat datax, arma::mat datas, arma::vec b0, arma::mat B0, double c1, double d1, int m, int n, int k, int l, bool verbose)
{
  
  // int m = y.n_rows;
  //int n = y.n_cols;
  //int k =0,l=0;
  
  //if(x_intercept == FALSE)
    //k = floor(datax.n_cols/n); 
  //else
    //k = floor(datax.n_cols/n) + 1; //add an intercept column in model matrix
  
  //if(s_intercept == FALSE)
    //  l = floor(datas.n_cols/n);
  //else
    //l = floor(datas.n_cols/n) + 1; //add an intercept column in model matrix
  
  arma::mat z = y - 0.5; //start for z
  arma::mat zPrev = z;
  arma::cube X(k,m,n);  //model matrix fixed
  arma::cube S(l,m,n);  //model matrix random
  
  for(int i=0;i<n;i++)
  {
    X.slice(i) = subset_mat(&datax,i,n).t(); 
    S.slice(i) = subset_mat(&datas,i,n).t(); 
  }
  
  
  /// MCMC and burnins
  //  nsim  = 12000;
  //int burn  = 0;
  //if(burnin == TRUE)
  //{
    //burn = 0.25*nsim;
  //}
  
  int MCMC  = nsim;    ///Total number of simulation
  int count = MCMC/10;
  
  //// Prior Distributions and Initializations!
  /// -------------------------------------------------------------------------
  
  // Beta
  ////--------------------
  // arma::vec b0(k,arma::fill::zeros);
  // arma::mat B0(k,k,arma::fill::eye); B0 = 10*B0;
  arma::mat invB0 = arma::inv_sympd(B0);
  arma::mat invB0b0 = invB0*b0;
  arma::mat beta_out(k,MCMC); // k x 1 matrices
  beta_out.col(0) = b0;
  
  /// varphi2
  ////--------------------
  //c1 = 10;
  //d1 = 9;
  arma::vec varphi2(MCMC);
  varphi2(0) = 1.0;
  //1/R::rgamma(c1/2,d1/2);
  
  // Alpha
  ////-------------------
  arma::vec a0(l,arma::fill::zeros);
  arma::mat A0(l,l,arma::fill::eye);
  A0 = A0*varphi2(0);
  arma::cube alpha(l,n,MCMC);
  alpha.slice(0) = arma::mvnrnd(a0,A0,n);
  
  // w
  ////-------------------
  arma::mat w(Rcpp::rexp(m*n,1.0));
  w.reshape(m,n);
  
  
  /// Parameters for Normal-Exponential ALD mixture
  //double varp  = (1- 2*p + 2*pow(p,2))/((pow(p,2))*pow(1-p,2));
  double theta = (1-2*p)/(p*(1-p));
  double tau2  = 2/(p*(1-p));
  // Index parameter for GIG distribution
  double lambdaIndex=0.5;  
  
  
  // Setting the limits for sampling from Truncated Multivariate Normal
  // will replace with function!!
  
  //  arma::mat lowerLimits(m,n,arma::fill::zeros);
  //  arma::mat upperLimits(m,n,arma::fill::zeros);
  
  //    for (int j=0;j<m;j++)
  //    {
  //      for (int i=0;i<n;i++)
  //      {
  //        if (y(j,i) == 0)
  //        {
  //          lowerLimits(j,i) = R_NegInf;
  //          upperLimits(j,i) = 0;
  //        }
  //        
  //        else if(y(j,i) == 1)
  //        {
  //          lowerLimits(j,i) = 0;
  //          upperLimits(j,i) = R_PosInf;
  //        }
  //      }
  //    }
  
  
  Rcpp::Rcout << "Please wait while we're processing your request.\n";
  Rcpp::Rcout << "I recommend listening to Vienna by Billy Joel while you wait.\n";
  Rcpp::Rcout << "https://music.apple.com/in/album/vienna/158617952?i=158618071\n";
  //Rcpp::Rcout << "OR Everywhere by Fleetwood Mac.\n";
  //Rcpp::Rcout << "https://music.apple.com/in/album/everywhere/202271826?i=202272247\n";
  
  
  ///-------------------------------------------------------------------------
  //// GIBBS SAMPLING
  //--------------------------------------------------------------------------------
  // The sequence of sampling is important and follows Algorithm 5 in Chib and Carlin.
  //---------------------------------------------------------------------------------
  
  for(int sim=1;sim<MCMC;sim++)
  {
    //int sim = 1;
    
    if(verbose && sim+1 == count)
    {
      Rcpp::Rcout << "No. of sim: " << sim+1 << "\n";
      count += MCMC/10;
    }
    
    ////--------- Sample beta,z marginally of alpha in a block --------------
    //beta_out.col(sim) = 
    sampleBeta(&z,&X,&S,&w,varphi2(sim-1),tau2,theta,&invB0,&invB0b0,k,m,n,&beta_out,sim);
    
    //////--------- Sample z, marginally of alpha -----------------------------
    ///// Draws random numbers from trucnated MVN (Geweke, 1991).
    zPrev = z;
    //z = 
    sampleZ(&zPrev,&y,&X,beta_out.col(sim),&S,theta,&w,varphi2(sim-1),tau2,m,n,&z);
    
    ////---------- Sample alpha conditionally on beta,z -------------------
    //alpha.slice(sim) = 
    sampleAlpha(&z,&X,&S,beta_out.col(sim),&w,tau2,theta,varphi2(sim-1),l,m,n,&alpha,sim);
    //sampleAlphafast(&z,&X,&S,beta_out.col(sim),&w,tau2,theta,varphi2(sim-1),l,m,n,&alpha,sim);
    
    /////--------- Sample w ---------------------------
    //w = 
    sampleW(&z,&X,&S,beta_out.col(sim),&alpha,tau2,theta,lambdaIndex,k,m,n,&w,sim);
    
    ////--------- Sample varphi2 ---------------------
    varphi2(sim) = sampleVarphi2(&alpha,c1,d1,l,n,sim);
    
  }
  
  return (Rcpp::List::create(Rcpp::Named("Beta", beta_out.t()),Rcpp::Named("Alpha", alpha),Rcpp::Named("Varphi2", varphi2))); 
  // return (Rcpp::List::create(Rcpp::Named("w", w))); 
}


/////////// Combined blocked SAMPLER - end  //////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
