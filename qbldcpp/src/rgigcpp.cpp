// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//' @useDynLib qbldcpp
//' @importFrom Rcpp sourceCpp
#include <RcppArmadillo.h>
#include"qbld.h"
// [[Rcpp::depends(RcppArmadillo)]]
////////////////////////////////////
///// GIG function!

double mode(double lambda,double omega)
{
  if(lambda>=1)
    return((sqrt((lambda-1)*(lambda-1)+omega*omega)+(lambda-1))/omega);
  
  return(omega/(sqrt((1-lambda)*(1-lambda)+omega*omega)+(1-lambda)));
}

int rgig_noshift(arma::vec*out, int n, double lambda, int check, double omega, double alpha)
{
  double xm,nc,ym,um,s,t,U,V,X;
  
  t = 0.5*(lambda-1);
  s = 0.25*omega;
  
  xm = mode(lambda,omega);
  nc = t*log(xm) - s*(xm + 1/xm);
  ym = ((lambda+1) + sqrt((lambda+1)*(lambda+1) + omega*omega))/omega;
  um = exp(0.5*(lambda+1)*log(ym) - s*(ym + 1/ym) - nc);
  
  for(int i=0; i<n; i++)
  {
    do{
      U = um*arma::randu<double>();
      V = arma::randu<double>();
      X = U/V;
    } while ((log(V)) > (t*log(X) - s*(X+1/X)- nc));
    (*out)(i) = (check==1) ? (alpha/X) : (alpha*X);
  }
  return(0);
}

int rgig_shift(arma::vec*out, int n, double lambda, int check, double omega, double alpha)
{
  double xm,nc,s,t,U,V,X,a,b,c,p,q,fi,fak,y1,y2,uplus,uminus;
  
  t = 0.5*(lambda-1);
  s = 0.25*omega;
  
  xm = mode(lambda,omega);
  nc = t*log(xm) - s*(xm + 1/xm); 
  
  a = -(2*(lambda+1)/omega +xm);
  b = (2*(lambda-1)*xm/omega -1);
  c = xm;
  
  p = b - a*a/3;
  q = (2*a*a*a)/27 - (a*b)/3 + c;
  
  fi = acos(-q/(2*sqrt(-p*p*p/27)));
  fak = 2*sqrt(-p/3);
  y1 = fak*cos(fi/3) - a/3;
  y2 = fak*cos(fi/3 + (4./3.)*M_PI) - a/3;
  
  uplus = (y1-xm)*exp(t*log(y1) - s*(y1 + 1/y1) -nc);
  uminus =  (y2-xm)*exp(t*log(y2) - s*(y2 + 1/y2) -nc);
  
  for(int i=0; i<n; i++)
  {
    do{
      U = uminus + arma::randu<double>() * (uplus-uminus);
      V = arma::randu<double>();
      X = U/V + xm;
    } while ((X<=0) || ((log(V)) > (t*log(X) - s*(X+1/X)- nc)));
    (*out)(i) = (check==1) ? (alpha/X) : (alpha*X);
  }
  return(0);
}


int rgig_conc(arma::vec*out, int n, double lambda, int check, double omega, double alpha)
{
  arma::vec A(3);
  double Atot,k0,k1,k2,xm,x0,a,U,V,X,hx;
  
  if(lambda >=1 || omega > 1)
    Rcpp::stop("Invalid parameters.");
  
  xm = mode(lambda,omega);
  x0 = omega/(1-lambda);
  
  k0 = exp((lambda-1)*log(xm) - 0.5*omega*(xm + 1/xm));
  A(0) = k0*x0;
  
  if(x0 >= 2/omega)
  {
    k1 = 0;
    A(1) = 0;
    k2 = pow(x0,lambda-1);
    A(2) = k2*2*exp(-omega*x0/2)/omega;
  }
  
  else
  {
    k1 = exp(-omega);
    A(1) = (lambda==0) ? (k1*log(2/(omega*omega))) : ((k1/lambda)*(pow(2/omega,lambda) - pow(x0,lambda)));
    k2 = pow(2/omega,lambda-1);
    A(2) = k2*2*exp(-1)/omega;
  }
  
  Atot = A(0) + A(1) + A(2);
  
  for(int i=0; i<n; i++)
  {
    do{
      V = Atot*(arma::randu<double>());
      
      do{
        
        if(V <= A(0))
        {
          X = x0*V/A(0);
          hx = k0;
          break;
        }
        
        V -= A(0);
        if (V <= A(1)) {
          if (lambda == 0) {
            X = omega * exp(exp(omega)*V);
            hx = k1 / X;
          }
          else {
            X = pow(pow(x0, lambda) + (lambda / k1 * V), 1/lambda);
            hx = k1 * pow(X, lambda-1);
          }
          break;
        }
        
        V -= A(1);
        a = (x0 > 2/omega) ? x0 : 2/omega;
        X = -2/omega * log(exp(-omega/2 * a) - omega/(2*k2) * V);
        hx = k2 * exp(-omega/2 * X);
        break;
        
      } while(0);
      
      U = hx*(arma::randu<double>());
      
      if(log(U) <= (lambda-1)*log(X) - omega/2 * (X+1/X))
      {
        (*out)(i) = (check==1) ? (alpha/X) : (alpha*X);
        break;
      }
    } while(1);
  }
  return(0);
}


// [[Rcpp::export]]
arma::vec rgig(double n,double lambda,double a,double b)
{
  arma::vec out(n);
  int check = 0;
  
  if(n<=0||std::floor(n) !=n)
    Rcpp::stop("sample size 'n' must be a positive integer");
  
  if ( !(R_FINITE(lambda) && R_FINITE(b) && R_FINITE(a)) ||
       (b <  0. || a < 0)      || 
       (b == 0. && lambda <= 0.) ||
       (a == 0. && lambda >= 0.) ) 
    Rcpp::stop("Invalid Parameters");
  
  if(b==0)
  {
    if(lambda>0)
      return(Rcpp::rgamma(n,lambda,2/a));
    else
      return(1/Rcpp::rgamma(n,-lambda,2/a));
  }
  
  else if(a==0)
  {
    if(lambda>0)
      return(Rcpp::rgamma(n,lambda,2/b));
    else
      return(1/Rcpp::rgamma(n,-lambda,2/b));
  }
  
  else
  {
    
    if(lambda<0)
    {
      lambda = -lambda;
      check = 1;
    }
    double alpha = sqrt(b/a);
    double omega = sqrt(a*b);
    
    if(lambda > 2 || omega > 3)
    {
      //RoU shift
      rgig_shift(&out,n,lambda,check,omega,alpha);
      return(out);
    }
    
    if(lambda >= 1 - 2.25*a*b || omega > 0.2)
    {
      //RoU no shift
      rgig_noshift(&out,n,lambda,check,omega,alpha);
      return(out);
    }
    
    if(lambda>=0 && omega>0)
    {
      //log-concave
      rgig_conc(&out,n,lambda,check,omega,alpha);
      return(out);
    }
    
    Rcpp::stop("Invalid parameters.");
  }
  return(out);
}