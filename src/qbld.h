#ifndef QBLD_H
#define QBLD_H

#include <RcppArmadillo.h>

//GIG distbn functions
double mode(double lambda,double omega);
int rgig_noshift(arma::vec*out, int n, double lambda, int check, double omega, double alpha);
int rgig_shift(arma::vec*out, int n, double lambda, int check, double omega, double alpha);
int rgig_conc(arma::vec*out, int n, double lambda, int check, double omega, double alpha);
arma::vec rgig(double n,double lambda,double a,double b);

//AL distbn function
arma::vec raldmix(double n,double mu,double sigma,double p);

// Samplers

//sampleBeta - Block
int sampleBeta(arma::mat*z, arma::cube*X, arma::cube*S, arma::mat*w, double varphi2, double tau2, double theta, arma::mat*invB0, arma::mat*invB0b0, int k, int m, int n, arma::mat* beta,int sim);

//sampleZ - Block
int sampleZ(arma::mat*zprev, arma::mat*y, arma::cube*X, arma::vec beta, arma::cube*S, double theta, arma::mat*w, double varphi2, double tau2, int m, int n, arma::mat*z);
//double erf_my(double x, bool inverse);
double erfc_my(double x, bool inverse);
double normcdf(double x);
double norminv(double x);
double limits(double y,int which);
arma::vec rtruncnorm_gwk(arma::vec z0,arma::vec*mu,arma::mat*sigma,arma::mat*y,int m,int idx); 

//sampleAlpha
int sampleAlphafast(arma::mat*z, arma::cube*X, arma::cube*S, arma::vec beta, arma::mat*w, double tau2, double theta, double varphi2, int l, int m, int n,arma::cube*alpha,int sim);
int sampleAlpha(arma::mat*z, arma::cube*X, arma::cube*S, arma::vec beta, arma::mat*w, double tau2, double theta, double varphi2, int l, int m, int n,arma::cube*alpha,int sim);

//sampleW
int sampleW(arma::mat*z, arma::cube*X, arma::cube*S, arma::vec beta, arma::cube*alpha, double tau2, double theta, double lambda, int k, int m, int n,arma::mat*w,int sim);
  
//sampleVarphi2
double sampleVarphi2(arma::cube* alpha, double c1, double d1, int l, int n,int sim);

//sampleZ - Unblock
int sampleZ_2(arma::mat*y, arma::cube*X, arma::vec beta, arma::cube*S, arma::cube* alpha ,double theta, arma::mat*w, double varphi2, double tau2, int m, int n,arma::mat*z,int sim);

//sampleBeta - Unblock
int sampleBeta_2(arma::mat*z, arma::cube*X, arma::cube*S, arma::mat*w, arma::cube *alpha ,double varphi2, double tau2, double theta, arma::mat*invB0, arma::mat*invB0b0, int k, int m, int n,arma::mat* beta,int sim);

//FinalBlock
arma::mat subset_mat(arma::mat* X, int start, int j);
Rcpp::List qbldf(int nsim, double p, arma::mat y, arma::mat datax, arma::mat datas, arma::vec b0, arma::mat B0, double c1, double d1, int m, int n, int k, int l,bool verbose);
Rcpp::List qbldunblock(int nsim, double p, arma::mat y, arma::mat datax, arma::mat datas, arma::vec b0, arma::mat B0, double c1, double d1,int m, int n, int k, int l,bool verbose);

#endif