// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//' @useDynLib qbld
//' @importFrom Rcpp sourceCpp

////////////////////////////////////
///// DGIG function!

#include <Rcpp.h>
using namespace Rcpp;

double bessel_k (double x, double nu, bool islog, bool expon_scaled){
  
#define M_LNPI 1.14472988584940017414342735135      /* ln(pi) */
  
  double z;                   /* rescaled argument for K_nu() */
  double sz, t, t2, eta;      /* auxiliary variables */
  double d, u1t,u2t,u3t,u4t;  /* (auxiliary) results for Debye polynomials */
  double res;                 /* value of log(K_nu(x)) [= result] */
  
  /* rescale: we comute K_nu(z * nu) */
  z = x / nu;
  
  /* auxiliary variables */
  sz = hypot(1,z);   /* = sqrt(1+z^2) */
  t = 1. / sz;
  t2 = t*t;
  
  eta = (expon_scaled) ? (1./(z + sz)) : sz;
  eta += log(z) - log1p(sz);                  /* = log(z/(1+sz)) */
  
  /* evaluate Debye polynomials u_j(t) */
  u1t = (t * (3. - 5.*t2))/24.;
  u2t = t2 * (81. + t2*(-462. + t2 * 385.))/1152.;
  u3t = t*t2 * (30375. + t2 * (-369603. + t2 * (765765. - t2 * 425425.)))/414720.;
  u4t = t2*t2 * (4465125.
                   + t2 * (-94121676.
                   + t2 * (349922430.
                   + t2 * (-446185740.
                   + t2 * 185910725.)))) / 39813120.;
                   d = (-u1t + (u2t + (-u3t + u4t/nu)/nu)/nu)/nu;
                   
                   /* log(K_nu(x)) */
                   res = log(1.+d) - nu*eta - 0.5*(log(2.*nu*sz) - M_LNPI);
                   
                   return (islog ? res : exp(res));
}

// [[Rcpp::export]]
std::vector<double> dgig(std::vector<double> x, double a, double b, double p,bool log_density)
{
  int nx = x.size();
  p = fabs(p);
  if(nx==0)
    stop("Either x is NULL or NaNs produced.");
  if(a<0 || b<0 || p == R_PosInf || p == R_NegInf || (a==0 && p<=0) || (b==0 &&p>=0))
    stop("invalid parameters for GIG distribution.");
  
  std::vector<double> ans(nx);
  std::vector<double> log_ans(nx);
  double cons_log = 0;
  double cons = 0.5*pow(a/b,p/2)/bessel_k(sqrt(a*b),p,FALSE,FALSE);
  
  if(log_density==TRUE)
    cons_log = 0.5*pow(a/b,p/2)/bessel_k(sqrt(a*b),p,TRUE,FALSE);
  
  for(int i=0;i<nx;i++)
  {
    if(x[i]<=0||std::isnan(x[i]))
      stop("X can't be non positive.");
    
    ans[i] = cons*pow(x[i],p)*exp(-0.5*(a*x[i]+b/x[i]));
    
    if(log_density==TRUE)
      log_ans[i] = cons_log + (p-1)*x[i] - 0.5*(a*x[i]+b/x[i]);
  }
  
  if(log_density==TRUE)
    return(log_ans);
  
  return(ans);
}

