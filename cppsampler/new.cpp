// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

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

//double o(arma::cube *a)
//{
//  return(((*a).slice(1))(1,1));
//  }

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

std::vector<double> rgig(double lambda, double a, double b, int n) {
  
  std::vector<double> X(n);
  std::vector<double> y(3);
  std::vector<double> UVW(3);
  
  // we sample from the two parameter version of the GIG(alpha,omega)
 // double lambda = P;
  double omega = sqrt(a*b);
  int check = 0;
  double t=0,s=0,f1=0,f2=0;
  
  if(lambda<0)
  {
    lambda = lambda*-1;
    check = 1;
  }
  
  //if (n<=0)
    //Rcpp::stop("sample size n must be positive integer.");
  
  //if(a<0 || b<0 || P == R_PosInf || P == R_NegInf || (a==0 && P<=0) || (b==0 &&P>=0))
  // stop("invalid parameters for GIG distribution.");     check this line dude some of these are not needed!
  
  y[1] = sqrt(pow(omega,2) + pow(lambda,2)) - lambda;
  //y[1] = alpha;
  y[2] = lambda;
  
  // Find t
  y[0] = 1;
  double x = -psi(y);
  
  if((x >= 0.5) && (x <= 2))
    t = 1;
  else if(x > 2)
    t = sqrt(2 / (y[1] + lambda));
  else if(x < 0.5)
    t = log(4/(y[1] + 2*lambda));
  
  // Find s
  y[0] = -1;
  x = -psi(y);
  if((x >= 0.5) && (x <= 2))
    s = 1;
  else if(x > 2)
    s = sqrt(4/(y[1]*cosh(1) + lambda));
  else if(x < 0.5)
    s = minn(1/lambda, log(1 + 1/y[1] + sqrt(1/pow(y[1],2)+2/y[1])));
  
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
      UVW[0] = R::runif(0,1);
      UVW[1] = R::runif(0,1);
      UVW[2] = R::runif(0,1);
      
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


int sampleBeta(arma::mat*z, arma::cube*X, arma::cube*S, arma::mat*w, double varphi2, double tau2, double theta, arma::mat*invB0, arma::mat*invB0b0, int k, int m, int n, arma::mat* beta,int sim)
{
  arma::mat sumvar(k,k,arma::fill::zeros);
  arma::vec summean(k,arma::fill::zeros);
  // arma::mat vari(k,k,arma::fill::zeros);
  //  arma::vec meani(k,arma::fill::zeros);
  //arma::mat D1(m,m,arma::fill::zeros);
  arma::mat inv_omega(k,m,arma::fill::zeros);
  
  for(int i=0;i<n;i++)
  {
    //D1.diag() = tau2*((*w).col(i));
    inv_omega = ((*X).slice(i))*((varphi2*((((*S).slice(i)).t())*((*S).slice(i))) + arma::diagmat(tau2*((*w).col(i)))).i());
    
    //vari 
    sumvar  += inv_omega*(((*X).slice(i)).t());
    //meani 
    summean += inv_omega*((*z).col(i) - theta*((*w).col(i)));
    
    // sumvar  += vari;
    //  summean += meani;
  }
  
  sumvar  = ((*invB0) + sumvar).i();
    //arma::inv_sympd((*invB0) + sumvar);
  summean = sumvar*((*invB0b0) + summean);
  
  (*beta).col(sim) = arma::mvnrnd(summean,sumvar,1);
  return(0);
}

/*
double erf_my(double x, bool inverse) 
{
  if (inverse) {
    double ans = R::qnorm((x+1)/2,0,1,1,0) / sqrt(2);
    
    if(x < -1) ans = NA_REAL;
    if(x > 1) ans = NA_REAL;
    if(x == -1) ans = R_NegInf;
    if(x == 1) ans = R_PosInf;
    
    return(ans);
  } 
  else {
    return(2*(R::pnorm(x*sqrt(2),0,1,1,0)) - 1); 
  }
}
*/

double erfc_my(double x, bool inverse) 
{
  if (inverse) {
    double ans = R::qnorm(x/2,0,1,0,0)/sqrt(2);
    
    if(x <  0) ans = NA_REAL;
    if(x > 2) ans = NA_REAL;
    if(x == 2) ans = R_NegInf;
    if(x == 0) ans = R_PosInf;
    
    return(ans);
  } 
  else {
    return(2*(R::pnorm(x*sqrt(2),0,1,0,0)));
  }
}


double normcdf(double x)
{
  return(0.5*erfc_my(-x/sqrt(2),FALSE));
}

double norminv(double x)
{
  return(-sqrt(2)*erfc_my(2*x,TRUE));
}

double limits(double y,int which)
{
  if(which == 0)  //lowerlimits
  {
    if(y == 0) return(R_NegInf);
    else return(0);
  }
  
  else if(which==1) //upperlimits
  {
    if(y == 0) return(0);
    else return(R_PosInf);
  }
  return(0);
}

arma::vec rtruncnorm_gwk(arma::vec z0,arma::vec*mu,arma::mat*sigma,arma::mat*y,int m,int idx)
{
  double newLL=0, newUL=0, condMu=0, sq_condVar=0, mu_i=0, sigma11=0;
  arma::mat inv_sigma22(m-1,m-1,arma::fill::zeros);
  arma::vec u = arma::randu(m); 
  
  for (int i=0; i<m; i++)
  {
    mu_i = (*mu)(i);
    
    arma::vec mu_less_i = *mu;
    mu_less_i.shed_row(i);        
    
    sigma11 = (*sigma)(i,i);
    arma::rowvec sigma12 = (*sigma).row(i);
    sigma12.shed_col(i);
    
    arma::mat sigma22 = *sigma;
    sigma22.shed_row(i);
    sigma22.shed_col(i);
    inv_sigma22 = sigma22.i();
    
    sq_condVar = sqrt((sigma11 - sigma12*inv_sigma22*(sigma12.t())).eval()(0,0)); 
    
    arma::vec z_less_i = z0;
    z_less_i.shed_row(i);        
    
    condMu = (mu_i + sigma12*inv_sigma22*(z_less_i - mu_less_i)).eval()(0,0);
    
    newLL  = (limits((*y)(i,idx),0) - condMu)/sq_condVar;
    newUL= (limits((*y)(i,idx),1) - condMu)/sq_condVar;
    
    z0(i) =  condMu + (sq_condVar)*(norminv(u(i)*(normcdf(newUL) - normcdf(newLL)) + normcdf(newLL)));
    //z0[i] = z[i]
  }
  
  return(z0);
}

int sampleZ(arma::mat*zprev, arma::mat*y, arma::cube*X, arma::vec beta, arma::cube*S, double theta, arma::mat*w, double varphi2, double tau2, int m, int n, arma::mat*z)
{
  
  //arma::mat z(m,n,arma::fill::zeros);
  arma::vec meani(m,arma::fill::zeros);
  arma::mat VarZi(m,m,arma::fill::zeros);
  //arma::mat D1(m,m,arma::fill::zeros);
  
  for(int i=0; i<n; i++)
  {
    //D1.diag() = tau2*((*w).col(i));
    meani = (((*X).slice(i)).t())*beta + theta*((*w).col(i));
    VarZi = (varphi2*((((*S).slice(i)).t())*((*S).slice(i))) + arma::diagmat(tau2*((*w).col(i))));
    
    (*z).col(i) = rtruncnorm_gwk((*zprev).col(i),&meani,&VarZi,y,m,i);
  }
  return(0);
  
}


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
    invD1 = ((*S).slice(i))*diagmat(1/(tau2*((*w).col(i))));
    
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


int sampleW(arma::mat*z, arma::cube*X, arma::cube*S, arma::vec beta, arma::cube*alpha, double tau2, double theta, double lambda, int k, int m, int n,arma::mat*w,int sim)
{
  
  double tilde_eta = (pow(theta,2))/(tau2) + 2;
  double tilde_lambda = 0;
  //arma::mat w(m,n,arma::fill::zeros);
  
  for(int i=0; i<n; i++)
  {
    for(int j=0; j<m; j++)
    {
      tilde_lambda = (pow(((*z)(j,i) - ((((*X).slice(i)).col(j)).t())*beta - ((((*S).slice(i)).col(j)).t())*(((*alpha).slice(sim)).col(i))),2)/tau2).eval()(0,0);
      if(tilde_lambda < 0.00000001) tilde_lambda = 0.00000001;
      
      (*w)(j,i) = rgig(lambda,tilde_eta,tilde_lambda,1)[0];
    }
  }
  return(0);
}

//--------------------------------------------------------------------------

double sampleVarphi2(arma::cube* alpha, double c1, double d1, int l, int n,int sim)
{
  double sum = 0;
  
  sum = arma::accu(arma::square((*alpha).slice(sim)));
  
  //return(1/R::rgamma((n*l+c1)/2,(sum + d1)/2));
  return(1/arma::randg<double>(arma::distr_param((n*l+c1)/2,2/(sum + d1))));
}

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


// [[Rcpp::export]]

// y is the output variable, x is fixed, s is random 

Rcpp::List qbldcpp_f(int nsim, double p, arma::mat y, arma::mat datax, arma::mat datas, bool x_intercept, bool s_intercept, arma::vec b0, arma::mat B0, double c1, double d1,bool burnin)
{
  
  int m = y.n_rows;
  int n = y.n_cols;
  int k =0,l=0;
  
  if(x_intercept == FALSE)
    k = floor(datax.n_cols/n); 
  else
    k = floor(datax.n_cols/n) + 1; //add an intercept column in model matrix
  
  if(s_intercept == FALSE)
    l = floor(datas.n_cols/n);
  else
    l = floor(datas.n_cols/n) + 1; //add an intercept column in model matrix
  
  arma::mat z = y - 0.5; //start for z
  arma::mat zPrev = z;
  arma::cube X(k,m,n);  //model matrix fixed
  arma::cube S(l,m,n);  //model matrix random
  
  for(int i=0;i<n;i++)
  {
    X.slice(i) = subset_mat(&datax,i,n,x_intercept).t(); 
    S.slice(i) = subset_mat(&datas,i,n,s_intercept).t(); 
  }
  
  
  /// MCMC and burnins
  //  nsim  = 12000;
  int burn  = 0;
  if(burnin == TRUE)
  {
    burn = 0.25*nsim;
  }
  
  int MCMC  = burn + nsim;    ///Total number of simulation
  
  
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
  arma::mat w(Rcpp::rexp(m*n,1));
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
  
  
  Rcpp::Rcout << "The sampler generates about 1250-1350 samples a minute, please wait while we're processing your request.\n";
  Rcpp::Rcout << "I recommend listening to Vienna by Billy Joel while you wait for every 5000 samples.\n";
  Rcpp::Rcout << "https://music.apple.com/in/album/vienna/158617952?i=158618071\n";
  Rcpp::Rcout << "OR Everywhere by Fleetwood Mac, Make it 10,000 XD.\n";
  Rcpp::Rcout << "https://music.apple.com/in/album/everywhere/202271826?i=202272247\n";
  
  
  ///-------------------------------------------------------------------------
  //// GIBBS SAMPLING
  //--------------------------------------------------------------------------------
  // The sequence of sampling is important and follows Algorithm 5 in Chib and Carlin.
  //---------------------------------------------------------------------------------
  
  for(int sim=1;sim<MCMC;sim++)
  {
    //int sim = 1;
    
    // if(sim/10.0 == sim/10)
    //  Rcpp::Rcout << "No. of sim: " << sim << "\n";
    
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
  
  return (Rcpp::List::create(Rcpp::Named("Beta", beta_out),Rcpp::Named("Alpha", alpha),Rcpp::Named("Varphi2", varphi2))); 
  // return (Rcpp::List::create(Rcpp::Named("w", w))); 
}

