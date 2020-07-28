library(RcppDist)
library(Rcpp)
library(RcppArmadillo)
#library(mvtnorm)
library(microbenchmark)
#library(FinalBlock)

#rald_mix <- function(n,mu=0,sigma=1,p=0.5,plot=FALSE)
#{
 # if(length(n) != 1) stop("Provide integer sample size n.")
#  if(n <= 0 || n%%1 !=0) stop("Sample size must be a positive integer value.")
 # if(sigma <= 0) stop("sigma (scale parameter) must be a positive number.")
  #if(p >= 1 || p <= 0) stop("p must be a real number in (0,1).")
  #if(abs(mu) ==Inf) stop("mu (location parameter) must be a finite real number.")
  
#  theta <- (1-2*p)/(p*(1-p))
 # tau <- sqrt(2/(p*(1-p)))
  
#  z <- rexp(n,rate=1)
 # u <- rnorm(n,mean=0,sd=1)
  
  #r <- theta*z + tau* sqrt(z) * u
  #r <- r + rep(mu,n)
  
#  if(plot == TRUE) #creates a sample plot for the ALD
 #   plot(density(rald_mix(1e5,mu,sigma,p)),main=paste0("ALD(", mu,",",sigma,",",p,")"),col="blue")
  
#  return(r)
#}


#datagenR <- function(n,m,p)
#{
  ####################
  ## Data generation 
  #n = 500
  #m = 10
  #p = 0.25
 # alpha    = t(mvtnorm::rmvnorm(n,rep(0,2),diag(2))) #true params for random effects
#  x1       = matrix(1,nrow=m,ncol=n)        #attribute 1 - fixed effect variable 1 (intercept)
 # x2       = matrix(runif(m*n),nrow=m,ncol=n)  #attribute 2 - fixed effect variable 2
#  x3       = matrix(runif(m*n),nrow=m,ncol=n)  #attribute 3 - fixed effect variable 3
 # datas    = matrix(runif(m*n),nrow=m,ncol=n)    #variable of random effects 1            
#  datax    = cbind(x1,x2,x3) #collated data
#  beta     = matrix(c(-5,6,4),nrow=3,ncol=1)    #true params for fixed effects
#  z        = matrix(0,nrow=m,ncol=n) #latent variable
  
  ### 100*p th Quantile
  ### ------------------
  
  #p = 0.25
 # epsilon  = matrix(rald_mix(m*n,0,1,p),nrow=m,ncol=n) ### ALD distbn on error terms
  
  #alpha[1,i]*1 is for the intercept in random effects
  #for(i in 1:n)
   # z[,i] = cbind(x1[,i],x2[,i],x3[,i])%*%beta + alpha[1,i] + datas[,i]*alpha[2,i] + epsilon[,i]
  
  ## Creating binary indicator outcome variables
  #y = matrix(0,nrow=m,ncol=n)
  
  #for(i in 1:n)
  #{
   # for(j in 1:m)
    #{
     # if(z[j,i] > 0)
      #  y[j,i] = 1
    #}
  #}
  #return(list(y,datax,datas))
#}


sourceCpp("./m.cpp") ## data generation is here
sourceCpp("./old.cpp") ##sampler is here
sourceCpp("./new.cpp"); #fast sampler is here

#dataR <- datagenR(500,10,0.25)
datacpp <- datagen(500,10,0.25) 

#x_r <- dataR[[2]]
x_cpp <- datacpp[[2]]

#y_r <- dataR[[1]]
y_cpp <- datacpp[[1]]
  
#s_r <- dataR[[3]]
s_cpp <- datacpp[[1]]
  

nsim <- 400 
x_intercept <- FALSE
s_intercept <- TRUE
b0 <- rep(0,3)
B0 <- 10*diag(3)
c1 <- 10
d1 <- 9
p <- 0.25

#runR <- qbldcpp(nsim, p, y_r, x_r, s_r, x_intercept, s_intercept, b0, B0, c1, d1)
#runcpp <- qbldcpp(nsim, p, y_cpp, x_cpp, s_cpp, x_intercept, s_intercept, b0, B0, c1, d1)
mbm <- microbenchmark(old = qbldcpp(nsim, p, y_cpp, x_cpp, s_cpp, x_intercept, s_intercept, b0, B0, c1, d1),new = qbldcpp_f(nsim, p, y_cpp, x_cpp, s_cpp, x_intercept, s_intercept, b0, B0, c1, d1),times=1)
print(mbm)
### change nsim! 500 samples just to test.
### compare with finalblock

##mbm <- microbenchmark(r = finalblock(400,0.25), cpp=qbldcpp(nsim, p, y_r, x_r, s_r, x_intercept, s_intercept, b0, B0, c1, d1),times=1)
##print(mbm)