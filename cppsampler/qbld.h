#ifndef QBLD_H
#define QBLD_H

#include <Rcpp.h>

//GIG distbn functions
double minn(double z, double zz);
double psi(std::vector<double> k);
double dpsi(std::vector<double> m);
double g(double x,double sd,double td,double f1,double f2);
std::vector<double> rgig(double P, double a, double b, int n);

//AL distbn function
arma::vec rald_mix(double n,double mu,double sigma,double p);

// Generic data generator for the model 
Rcpp::List datagen(int n,int m,double p);

// Samplers

//sampleBeta - Block
arma::vec sampleBeta(arma::mat*z, arma::cube*X, arma::cube*S, arma::mat*w, double varphi2, double tau2, double theta, arma::mat*invB0, arma::mat*invB0b0, int k, int m, int n);

//sampleZ - Block
arma::mat sampleZ(arma::mat*zprev, arma::mat*y, arma::cube*X, arma::vec beta, arma::cube*S, double theta, arma::mat*w, double varphi2, double tau2, int m, int n);
double erf_my(double x, bool inverse);
double erfc_my(double x, bool inverse);
double normcdf(double x);
double norminv(double x);
double limits(arma::mat*y,int which,int idx1,int idx2);
arma::vec rtruncnorm_gwk(arma::vec z0,arma::vec*mu,arma::mat*sigma,arma::mat*y,int m,int idx);  

//sampleAlpha
arma::mat sampleAlpha(arma::mat*z, arma::cube*X, arma::cube*S, arma::vec beta, arma::mat*w, double tau2, double theta, double varphi2, int l, int m, int n);

//sampleW
arma::mat sampleW(arma::mat*z, arma::cube*X, arma::cube*S, arma::vec beta, arma::mat alpha, double tau2, double theta, double lambda, int k, int m, int n);

//sampleVarphi2
double sampleVarphi2(arma::mat alpha, double c1, double d1, int l, int n);

#endif