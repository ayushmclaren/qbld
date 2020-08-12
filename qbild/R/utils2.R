model.qbld2 <- function(nsim, p=0.25, y, fixed, random, b0, B0, c1, d1,
                       fixed_intercept=TRUE, random_intercept=TRUE, burn=0,
                       method="block", names_fixed, summarize=FALSE, verbose=FALSE)
{
  
  if(missing(nsim)||nsim==0||nsim%%1!=0)
    stop("n must be specified correctly.")
  
  if(missing(y)||missing(fixed)||missing(random)
     ||(sum(y[is.na(y)==TRUE]) > 0)||(sum(fixed[is.na(fixed)==TRUE]) > 0)
     ||(sum(random[is.na(random)==TRUE]) > 0))
    stop("Either the data isn't provided or contains NAs.")
  
  if (!(all(sapply(fixed, is.numeric))) || !(all(sapply(random, is.numeric)))) {
    stop ("Data frame contains non-numeric values. Consider conversion of Factor variables.")
  }
  
  
  if(burn<0 || burn>1)
    stop("burn should be a number in [0,1].")
  
  m = nrow(y)
  n = ncol(y)
  
  if(fixed_intercept)
    k = floor(ncol(fixed)/n) + 1 ##add an intercept column in model matrix
  else
    k = floor(ncol(fixed)/n) 
  
  if(random_intercept)
    l = floor(ncol(random)/n) + 1 ##add an intercept column in model matrix
  else
    l = floor(ncol(random)/n)
  
  if(missing(b0))
    b0 = rep(0,k) #start from 0
  if(missing(B0))
    B0 = diag(k) 
  
  if(missing(c1))
    c1 = 9
  if(missing(d1))
    d1 = 10
  
  if(missing(names_fixed)||is.null(names_fixed))
  {
    if(fixed_intercept)
      varnames.fixed = c(paste("beta", 1:(k-1), sep = ""),"Varphi2") #Naming Beta1 to Betak
    
    else
      varnames.fixed = c(paste("beta", 1:k, sep = ""),"Varphi2") #Naming Beta1 to Betak
  }
  else
    varnames.fixed = c(names_fixed,"Varphi2") #Take names form the data frame
  
  if(fixed_intercept)
    varnames.fixed = c("Intercept",varnames.fixed) #If intercept
  
  if(regexec("block",method,ignore.case=TRUE)[[1]][1]==1) #blocked
    out = (qbldf(nsim, p, y, fixed, random, b0, B0, c1, d1, m, n, k, l, verbose))
  
  if(regexec("unblock",method,ignore.case=TRUE)[[1]][1]==1) #unblocked
    out = (qbldunblock(nsim, p, y, fixed, random, b0, B0, c1, d1, m, n, k, l, verbose))
  
  if(is.null(out))
    stop("Method is either 'block' or 'unblock'.")
  
  out = make.qbld(out, p, nsim, burn, varnames.fixed, method) 
  
  if(summarize==TRUE)
  {
    beta  = matrix(colMeans(out[[1]]),ncol=1)
    alpha = rowMeans(out[[2]],dims=2)
    varphi2 = mean(out[[3]])
    
    model_fit = mofit(y, fixed, random, beta, alpha, varphi2, p, 
                      fixed_intercept, random_intercept, k, l, m, n)
    
    print(summary(out))
    
    cat("\n3. Model Selection Criterion\n")
    cat("Log likelihood =", model_fit[1], "\n")
    cat("AIC =", model_fit[2], "\n")
    cat("BIC =", model_fit[3], "\n")
  }
  
  return(out)
}