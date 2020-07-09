#' @export

finalblock <- function(mc,p)
{


### just in case, don't know if it will build automatically so!
#source("./R/rald_mix.R")
#source("./R/sampleBeta.R")
#source("./R/sampleZ.R")
#source("./R/sampleAlpha.R")
#source("./R/sampleW.R")
#source("./R/sampleVarphi2.R")

####################
## Data generation 
n = 500
m = 10

alpha    = t(mvtnorm::rmvnorm(n,rep(0,2),diag(2))) #true params for random effects
x1       = matrix(1,nrow=m,ncol=n)        #attribute 1 - fixed effect variable 1 (intercept)
x2       = matrix(runif(m*n),nrow=m,ncol=n)  #attribute 2 - fixed effect variable 2
x3       = matrix(runif(m*n),nrow=m,ncol=n)  #attribute 3 - fixed effect variable 3
datas    = matrix(runif(m*n),nrow=m,ncol=n)    #variable of random effects 1            
datax    = cbind(x1,x2,x3) #collated data
beta     = matrix(c(-2,3,4),nrow=3,ncol=1)    #true params for fixed effects
z        = matrix(0,nrow=m,ncol=n) #latent variable

### 25th Quantile
### ------------------

#p = 0.25
epsilon  = matrix(rald_mix(m*n,0,1,p),nrow=m,ncol=n) ### ALD distbn on error terms

#alpha[1,i]*1 is for the intercept in random effects
for(i in 1:n)
  z[,i] = cbind(x1[,i],x2[,i],x3[,i])%*%beta + alpha[1,i] + datas[,i]*alpha[2,i] + epsilon[,i]

## Creating binary indicator outcome variables
y = matrix(0,nrow=m,ncol=n)

for(i in 1:n)
{
  for(j in 1:m)
  {
    if(z[j,i] > 0)
      y[j,i] = 1
  }
}

#######################
## Load Data - no need to load few things again, already in this environment (only for testing)

l     = 2                       #num of attributes-random effects s1 and 1(intercept)
n     = ncol(y)                 #number of people at a time
m     = nrow(y)                 # number of time periods
z     = y - 0.5 
######################################################################################

####################################################

## SAMPLING starts below ###

####################################################
k  = ncol(datax)/n  ### for number of attributes in fixed effects x1, x2, x3
x  = array(0, dim = c(k,m,n))
s  = array(0, dim = c(l,m,n))

for (i in 1:n)
{
  x[,,i] = t(cbind(x1[,i],x2[,i],x3[,i]))   ##model matrix fixed + intercept
  s[,,i] = t(cbind(matrix(1,nrow=m,ncol=1), datas[,i])) ##model matrix random + intercept
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
#mc      = 12000
burn    = 0.25*mc
MCMC    = burn + mc         ####Total number of simulation

#Storage matrices
#z       = zeros(m,n,N);    # only if i want store the latent values
betap     = matrix(0,nrow=k,ncol=MCMC) # k x 1 matrices
alpha     = array(0, dim=c(l,n,MCMC)) # l x n matrices
varphi2   = matrix(0,nrow=MCMC,ncol=1)    #scalars


### Initial values of z, alpha, w, varphi2
alpha[,,1]      = t(mvtnorm::rmvnorm(n,rep(0,l),diag(l)))
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
  
  print(nsim)
  #--------- Sample beta,z marginally of alpha in a block --------------
  betap[,nsim] = sampleBeta(z,x,s,w,varphi2[nsim-1,1],tau2,theta,invB0,invB0b0)
  
  
  #--------- Sample z, marginally of alpha -----------------------------
  # Draws random numbers from trucnated MVN (Geweke, 1991).
  zPrev = z
  z = sampleZ(zPrev,y,x,betap[,nsim],s,theta,w,varphi2[nsim-1,1],tau2,lowerLimits,upperLimits)
  
  
  #---------- Sample alpha conditionally on beta,z -------------------
  alpha[,,nsim] = sampleAlpha(z,x,s,betap[,nsim],w,tau2,theta,varphi2[nsim-1,1])
  
  
  #---------- Sample w ---------------------------
  w = sampleW(z,x,s,betap[,nsim],alpha[,,nsim],tau2,theta,lambdaIndex)
  
  
  #---------- Sample varphi2 ---------------------
  varphi2[nsim,1] = sampleVarphi2(alpha[,,nsim],c1,d1)
}

return(list(betap,varphi2))

}



