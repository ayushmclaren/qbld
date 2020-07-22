// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
//#include "RcppEigen.h"
//#include "RcppDist.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo,RcppEigen)]]

// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//

 
 /*
int yo(arma::mat*A, arma::mat*D,int m,int n)
  {
    Eigen::MatrixXd B = Eigen::Map<Eigen::MatrixXd>(*A.memptr(),m,n);
    Eigen::MatrixXd C = Eigen::MatrixXd(A.n_rows, A.n_rows).setZero().selfadjointView<Eigen::Lower>().rankUpdate(B);
  //  Rcpp::Rcout << std::endl << C << std::endl;
    arma::mat c  = arma::mat(C.data(), C.rows(), C.cols(),false, false);
    //Rcpp::Rcout << std::endl << c << std::endl;
    *D = c;
    //Rcpp::Rcout << std::endl << *D << std::endl;
    return(0);
  }
*/
 
 
 // [[Rcpp::export]]

 arma::mat ok1(int tt)
 {
   int m = 10;
   int n = 5000;
   
  // arma::mat out(m,n,arma::fill::zeros);
   
   return(arma::mvnrnd(arma::vec(10,arma::fill::ones),arma::eye(10,10),5000));
  // return(out);
 }
 
 
 // [[Rcpp::export]]
 
Rcpp::RObject ok2(int n, arma::vec mu,arma::mat sigma,int ncores,bool ischol)
 {
   //int m = 10;
   //int n = 5000;
   
   //arma::mat out(m,n,arma::fill::zeros);
   Rcpp::Environment pkg = Rcpp::Environment::namespace_env("mvnfast");
   Rcpp::Function f = pkg["rmvnCpp"];
   //bool fal = FALSE;
   
   return(f(n,mu,sigma,ncores,ischol,NULL));
 }
 
 
 
 
 