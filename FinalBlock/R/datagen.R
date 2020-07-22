#' @export


datagenR <- function(n,m,p)
{
####################
## Data generation 
#n = 500
#m = 10

alpha    = t(mvtnorm::rmvnorm(n,rep(0,2),diag(2))) #true params for random effects
x1       = matrix(1,nrow=m,ncol=n)        #attribute 1 - fixed effect variable 1 (intercept)
x2       = matrix(runif(m*n),nrow=m,ncol=n)  #attribute 2 - fixed effect variable 2
x3       = matrix(runif(m*n),nrow=m,ncol=n)  #attribute 3 - fixed effect variable 3
datas    = matrix(runif(m*n),nrow=m,ncol=n)    #variable of random effects 1            
datax    = cbind(x1,x2,x3) #collated data
beta     = matrix(c(-5,6,4),nrow=3,ncol=1)    #true params for fixed effects
z        = matrix(0,nrow=m,ncol=n) #latent variable

### 100*p th Quantile
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
return (list(y,datax,datas))
}