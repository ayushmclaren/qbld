//' @useDynLib GIGpack
//' @importFrom Rcpp sourceCpp

#include <Rcpp.h>
using namespace Rcpp;

double minn(double z, double zz)
{
  if(z<=zz)
    return(z);
  return(zz);
}

double psi(std::vector<double> k)
{
  return(-k[1]*(cosh(k[0])-1) -k[2]*(exp(k[0])-k[0]-1));
}

double dpsi(std::vector<double> m)
{
  return(-m[1]*(sinh(m[0])) -m[2]*(exp(m[0])-1));
}

double g(double x,double sd,double td,double f1,double f2)
{
  double a = 0;
  double b = 0;
  double c = 0;

  if((x >= -sd) && (x <= td))
    a = 1;

  else if(x > td)
    b = f1;

  else if(x < -sd)
    c = f2;

  return(a+b+c);
}

// [[Rcpp::export]]
std::vector<double> rgig(double P, double a, double b, int n) {

  std::vector<double> X(n);
  std::vector<double> y(3);
  std::vector<double> UVW(3);

  // we sample from the two parameter version of the GIG(alpha,omega)
  double lambda = P;
  double omega = sqrt(a*b);
  int check = 0;
  double t=0,s=0,f1=0,f2=0;

  if(lambda<0)
  {
    lambda = lambda*-1;
    check = 1;
  }

  if (n<=0)
    stop("sample size n must be positive integer.");

  //if(a<0 || b<0 || P == R_PosInf || P == R_NegInf || (a==0 && P<=0) || (b==0 &&P>=0))
   // stop("invalid parameters for GIG distribution.");     check this line dude some of these are not needed!

  double alpha = sqrt(pow(omega,2) + pow(lambda,2)) - lambda;
  y[1] = alpha;
  y[2] = lambda;

  // Find t
  y[0] = 1;
  double x = -psi(y);

  if((x >= 0.5) && (x <= 2))
    t = 1;
  else if(x > 2)
    t = sqrt(2 / (alpha + lambda));
  else if(x < 0.5)
    t = log(4/(alpha + 2*lambda));

  // Find s
  y[0] = -1;
  x = -psi(y);
  if((x >= 0.5) && (x <= 2))
    s = 1;
  else if(x > 2)
    s = sqrt(4/(alpha*cosh(1) + lambda));
  else if(x < 0.5)
    s = minn(1/lambda, log(1 + 1/alpha + sqrt(1/pow(alpha,2)+2/alpha)));

  y[0] = t;
  double eta = -psi(y);
  double zeta = -dpsi(y);
  y[0] = s;
  double theta = -psi(y);
  double xi = dpsi(y);

  double p = 1/xi;
  double r = 1/zeta;
  double td = t - r*eta;
  double sd = s - p*theta;
  double q = td + sd;

  for(int i=0;i<n;i++)
  {
    while(1)
    {
      UVW[0] = (float)std::rand()/RAND_MAX ;
      UVW[1] = (float)std::rand()/RAND_MAX ;
      UVW[2] = (float)std::rand()/RAND_MAX ;

      if(UVW[0] < (q / (p + q + r)))
        X[i] = -sd + q*UVW[1];

      else if(UVW[0] < ((q + r) / (p + q + r)))
        X[i] = td - r*log(UVW[1]);

      else
        X[i] = -sd + p*log(UVW[1]);

      f1 = exp(-eta - zeta*(X[i]-t));
      f2 = exp(-theta + xi*(X[i]+s));

      y[0] = X[i];
      if( UVW[2]*g(X[i], sd, td, f1, f2) <= exp(psi(y)) )
        break;
    }
  }

  for(int j=0;j<n;j++)
  {
    X[j] = exp(X[j]) * (lambda / omega + sqrt(1 + pow((lambda/omega),2)));

    if(check==1)
      X[j] = 1/X[j];

    X[j] = X[j] * sqrt(b/a);
  }
  return(X);
}
