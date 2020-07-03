##################################################################
### This is a skeleton, please DON'T RUN for now!
### DO NOT BUILD THIS PACKAGE, for internal purpose only!
##################################################################
### This essentially brings all the sampling steps together!
### The whole thing, with all the samplers will be migrated
### to Rcpp/Armadillo at a later stage for faster and efficient sampling!
##################################################################
### Convergence Diagnostics using mcmcse!
### as well as other outputs will be dealt with next!
##################################################################

# Loading data
datax = read.csv("./datax.csv")
datas = read.csv("./datas.csv")
dataalpha = read.csv("./dataalpha.csv")
y = read.csv("./y25.csv")

l     = 2                       #num of attributes-random effects s1 and 1(intercept)
n     = ncol(y)                  #number of people at a time
m     = nrow(y)                 # number of time periods
x1    = datax[ ,1:n]        #attribute 1 - fixed effect variable 1
x2    = datax[ ,n+1:2*n]     #attribute 1 - fixed effect variable 1
x3    = datax[ ,2*n+1 : 3*n]    #attribute 1 - fixed effect variable 1
s1    = datas                #variable of random effects 1
alpha = dataalpha              #true params for random effects
#####
beta     = matrix(c(-5,6,4),nrow=3,ncol=1)    #true params for fixed effects

####################################################

## SAMPLING starts below ###

####################################################

k  = ncol(datax)/n  ### for number of attributes in fixed effects x1, x2, x3
x  = array(0, dim = c(k,m,n))
s  = array(0, dim = c(l,m,n))

for (i in 1:n)
{
  x[,,i] = t(cbind(x1[,i],x2[,i],x3[,i]))
  s[,,i] = t(cbind(matrix(1,nrow=m,ncol=1), s1[,i]))
}

#### Prior Distributions
#### -------------------------------------------------------------------------
## Beta
##--------------------

b0 = matrix(0,nrow=k,ncol=1)
B0 = 10*diag(k)
invB0   = solve(B0)
invB0b0 = invB0%*%b0

## Alpha
###-------------------

a0 = matrix(0,nrow=l,ncol=1)
A0 = diag(l)

### varphi2
#--------------------

c1 = 10
d1 = 9

### MCMC and burnins
mc      = 12000
burn    = 0.25*mc
MCMC    = burn + mc         ####Total number of simulation

#Storage matrices
#z       = zeros(m,n,N);    # only if i want store the latent values
betap     = matrix(0,nrow=k,ncol=MCMC) # k x 1 matrices
alpha     = array(0, dim=c(l,n,MCMC)) # l x n matrices
varphi2   = matrix(0,nrow=MCMC,ncol=1)    #scalars


### Initial values of z, alpha, w, varphi2
z               = y - 0.5
#z               = z25
alpha[,,1]      = t(rmvnorm(n,rep(0,l),diag(l)))
w               = abs(matrix(rnorm(m*n,mean=2,sd=1), nrow=m,ncol=n))
varphi2[1,1]    = 1

#### Parameters for Normal-Exponential mixture
varp  = (1-2*p + 2*p^2)/((p^2)*(1-p)^2)
theta = (1-2*p)/(p*(1-p))
tau2  = 2/(p*(1-p))

lambdaIndex=0.5  # Index parameter for GIG distribution

## Setting the limits for sampling from Truncated Multivariate Normal

lowerLimits = matrix(0,nrow=m,ncol=n)
upperLimits = matrix(0,nrow=m,ncol=n)

for (j in 1:m)
{
  for (i in 1:n)
  {
    if (y[j,i] == 0)
    {
      lowerLimits[j,i] = -Inf
      upperLimits[j,i] = 0
    }

    else if (y[j,i] ==1)
    {
      lowerLimits[j,i] = 0
      upperLimits[j,i] = Inf
    }
  }
}


###-------------------------------------------------------------------------
#### GIBBS SAMPLING
##--------------------------------------------------------------------------------
## The sequence of sampling is important and follows Algorithm 5 in Chib and Carlin.
##---------------------------------------------------------------------------------

for(nsim in 2:MCMC)
{

  #---------- Sample beta, alpha in a block ------------------------
  betap[,nsim] = sampleBeta_2(z,x,s,w,alpha[,,nsim-1],varphi2[nsim-1,1],tau2,theta,invB0, invB0b0)


  #---------- Sample alpha conditionally on beta,z -------------------
  alpha[,,nsim] = sampleAlpha(z,x,s,betap[,nsim],w,tau2,theta,varphi2[nsim-1,1])


  #---------- Sample w ---------------------------
  w = sampleW(z,x,s,betap[,nsim],alpha[,,nsim],tau2,theta,lambdaIndex)


  #---------- Sample varphi2 ---------------------
  varphi2[nsim,1] = sampleVarphi2(alpha[,,nsim],c1,d1)


  #---------- Sample z, marginally of alpha ---------------------------
  z = sampleZ_2(y,x,s,betap[,nsim],alpha[,,nsim-1],varphi2[nsim,1],w,tau2,theta,lowerLimits,upperLimits)

  }
