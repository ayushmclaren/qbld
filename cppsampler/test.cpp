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
 
 arma::mat subset_mat(arma::mat* X, int start, int j, bool intercept) 
 {
    arma::uvec IDX = arma::regspace<arma::uvec>(start,  j,  (*X).n_cols-1);
    
    if(intercept == TRUE)
    {
       int rw = (*X).n_rows;
       return(arma::join_rows(arma::vec(rw,arma::fill::ones),(*X).cols(IDX)));
    }
    
    return (*X).cols(IDX);
 }
 
    arma::mat ok(int tt)
    {
       int m = 3;
       int n = 4;
       arma::mat w(Rcpp::rexp(m*n,1.0));
       w.reshape(m,n);
       return(w);
    }
 
 
 
 