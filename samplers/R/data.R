## Data Generation: QBLD
## This file creates data sets for  Bayesian quantile regression for longitudinal
##  data with binary responses.

## Simulation Study
##
##Data Generation

n = 500
m = 10

alpha = t(rmvnorm(n,rep(0,2),diag(2)))
x1       = matrix(1,nrow=m,ncol=n)
x2       = matrix(runif(m*n),nrow=m,ncol=n)
x3       = matrix(runif(m*n),nrow=m,ncol=n)
s        = matrix(runif(m*n),nrow=m,ncol=n)
x        = cbind(x1,x2,x3)


## Storing the data
datax     = x
datas     = s
dataalpha = alpha

write.csv(datax,'datax.csv')
write.csv(datas,'datas.csv')
write.csv(dataalpha,'dataalpha.csv')

##################################################
#load('datax')
#load('datas')
#load('dataalpha')

# Loading data

n     = 500                  #number of people at a time
m     = 10                  # number of time periods
x1    = datax[ ,1:n]        #attribute 1 - fixed effect variable 1
x2    = datax[ ,n+1 : 2*n]     #attribute 1 - fixed effect variable 1
x3    = datax[ ,2*n+1 : 3*n]    #attribute 1 - fixed effect variable 1
s1    = datas                #variable of random effects
alpha = dataalpha              #true params for random effects
#####
beta     = matrix(c(-1,2,3),nrow=3,ncol=1)       #true params for fixed effects
z        = matrix(0,nrow=m,ncol=n) #latent variable

### 25th Quantile
### ------------------

p = 0.25
epsilon  = matrix(rald_mix(m*n,0,1,p),nrow=m,ncol=n) ### ALD distbn on error terms

#alpha[1,i]*1 is for the intercept in random effects
for(i in 1:n)
z[,i] = cbind(x1[,i],x2[,i],x3[,i])%*%beta + alpha[1,i] + s1[,i]*alpha[2,i] + epsilon[,i]

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

### Saving the data
y25 = y
write.csv(y25,'y25.csv')
z25 = z
write.csv(z25,'z25.csv')


### 50th Quantile
### ------------------

p = 0.50
epsilon  = matrix(rald_mix(m*n,0,1,p),nrow=m,ncol=n) ### ALD distbn on error terms

#alpha[1,i]*1 is for the intercept in random effects
for(i in 1:n)
  z[,i] = cbind(x1[,i],x2[,i],x3[,i])%*%beta + alpha[1,i] + s1[,i]*alpha[2,i] + epsilon[,i]

## Creating binary indicator outcome variables
y = matrix(0,nrow=m,ncol=n)

for(i in 1:n)
{
  for(j in 1:m)
  {
    if (z[j,i] > 0)
    y[j,i] = 1
  }
}

### Saving the data
y50 = y;
write.csv(y50,'y50.csv')
z50 = z;
write.csv(z50,'z50.csv')



### 75th Quantile
### ------------------

p = 0.75
epsilon  = matrix(rald_mix(m*n,0,1,p),nrow=m,ncol=n) ### ALD distbn on error terms

#alpha[1,i]*1 is for the intercept in random effects
for(i in 1:n)
  z[,i] = cbind(x1[,i],x2[,i],x3[,i])%*%beta + alpha[1,i] + s1[,i]*alpha[2,i] + epsilon[,i]

## Creating binary indicator outcome variables
y = matrix(0,nrow=m,ncol=n)

for(i in 1:n)
{
  for(j in 1:m)
  {
    if (z[j,i] > 0)
    y[j,i] = 1
  }
}

### Saving the data
y75 = y;
write.csv(y75,'y75.csv')
z75 = z;
write.csv(z75,'z75.csv')

###--------------------------------------------------------------------------
