// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//' @useDynLib Gqbldcpp
//' @importFrom Rcpp sourceCpp

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include "qbld.h"
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

////////////////////////////////////////////////////
// Data Generation
// [[Rcpp::export]]
Rcpp::List datagen(int n,int m,double p)
{
  arma::vec alph_mean(2,arma::fill::zeros);
  arma::mat alph_cov(2,2,arma::fill::eye);
  arma::mat dataalpha = arma::mvnrnd(alph_mean,alph_cov,n);
  
  arma::mat x1(m,n,arma::fill::ones);
  arma::mat x2(m,n,arma::fill::randu);
  arma::mat x3(m,n,arma::fill::randu);
  arma::mat datax = arma::join_rows(x1,x2,x3);
  arma::mat datas(m,n,arma::fill::randu);
  
  arma::mat beta_true(3,1);
  beta_true(0,0) = -5.0;
  beta_true(1,0) = 6.0;
  beta_true(2,0) = 4.0;
  
  arma::mat z(m,n);
  arma::mat y(m,n);
  
  //double p = 0.25;
  arma::mat epsilon(raldmix(m*n,0,1,p));
  epsilon.reshape(m,n);
  
  
  for(int i=0; i<n; i++)
  {
    z.col(i) =  (arma::join_rows(x1.col(i),x2.col(i),x3.col(i)))*beta_true + dataalpha(0,i) + (datas.col(i))*dataalpha(1,i) + epsilon.col(i);
  }
  
  for(int i=0;i<n;i++)
  {
    for(int j=0;j<m;j++)
    {
      if(z(j,i) > 0)
        y(j,i) = 1;
      
      else
        y(j,i) = 0;
    }
  }
  
  return (Rcpp::List::create(Rcpp::Named("y", y),Rcpp::Named("x", datax),Rcpp::Named("s", datas)));
}

