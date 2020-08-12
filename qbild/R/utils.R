#' @rdname qbild
#' @export

model.qbld <- function(data, id="id", fixed_formula, random_formula, nsim, p=0.25, b0, B0, 
                       c1, d1, method="block", burn=0, summarize=FALSE, verbose=FALSE)
{
  
  if(missing(nsim)||nsim==0||nsim%%1!=0)
    stop("n must be specified correctly.")
  
  if(missing(data)||(sum(is.na(data)) > 0)|| !is.data.frame(data))
    stop("Either the data isn't provided as a data frame or contains NAs.")
  
  if(is.null(names(data)))
    stop("objects in data.frame must have a name")
  
  if(burn<0 || burn>1)
    stop("burn should be a number in [0,1].")
  
  if(missing(fixed_formula) || missing(random_formula))
    stop("Please provide valid formulas for both fixed and random effects.")
  
  must_convert<-sapply(data,is.factor)       # logical vector telling if a variable needs to be displayed as numeric
  if(sum(must_convert))
  {
    df2 <-sapply(data[,must_convert,drop=FALSE],unclass) # data.frame of all categorical variables now displayed as numeric
    #colnames(df2) <-names(data)[must_convert]
    data<-cbind(data[,!must_convert],df2)  #converted back appropriately 
  }
  
  if(!missing(id)) {indvs<-as.vector(data[[id]])}
  if (missing(id)){ if (all(is.na(match(names(data), "id")))) stop ("id must be defined")
    else indvs<-as.vector(data$id)}
  
  
  m = sum(indvs==1) #no. of rows for my data by counting no. of entries for each individual 
  
  fixed_formula = formula(fixed_formula)
  expr1 = terms(fixed_formula, data=data)
  expr <- attr(expr1, "variables")
  var.names <- attr(expr1, "term.labels")
  
  if(any(is.na(match(all.vars(expr), names(data)))))
    stop("Variables in formula not contained in the data.frame")
  
  response <- as.character(expr[[2]])
  if(attr(expr1,"response")==0 || is.null(response))
    stop("Please provide response variable in fixed_formula correctly.")
  y = matrix(as.matrix(data[response]),nrow=m)
  n = ncol(y) # second dimension
  
  fixed_intercept = attr(expr1,"intercept")
  if(fixed_intercept == 1)
    var.names = c("Intercept",var.names)
  
  fixed = model.matrix(fixed_formula,data=data)
  fixed = matrix(fixed,nrow=m)
  k = floor(ncol(fixed)/n)
  
  random_formula = formula(random_formula)
  expr2 = terms(random_formula, data=data)
  random_intercept = attr(expr2,"intercept")
  random = model.matrix(random_formula,data=data)
  random = matrix(random,nrow=m)
  l = floor(ncol(random)/n)
  
  
  if(missing(b0))
    b0 = rep(0,k) #start from 0
  if(missing(B0))
    B0 = diag(k) 
  
  if(missing(c1))
    c1 = 9
  if(missing(d1))
    d1 = 10
  
  
  var.names = c(var.names,"Varphi2") #Take names form the data frame and add varphi2
  
  if(regexec("block",method,ignore.case=TRUE)[[1]][1]==1) #blocked
    out = (qbldf(nsim, p, y, fixed, random, b0, B0, c1, d1, m, n, k, l, verbose))
  
  if(regexec("unblock",method,ignore.case=TRUE)[[1]][1]==1) #unblocked
    out = (qbldunblock(nsim, p, y, fixed, random, b0, B0, c1, d1, m, n, k, l, verbose))
  
  if(is.null(out))
    stop("Method is either 'block' or 'unblock'.")
  
  out = make.qbld(out, p, nsim, burn, var.names, method) 
  
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



pald <- function(q,mu=0,sigma=1,p=0.5,lower.tail=TRUE)
{
  if(length(q) == 0) stop("Provide quantile q.")
  if(sigma <= 0) stop("sigma (scale parameter) must be a positive number.")
  if(p >= 1 || p <= 0) stop("p must be a real number in (0,1).")
  if(abs(mu) ==Inf) stop("mu (location parameter) must be a finite real number.")
  if(lower.tail != TRUE && lower.tail != FALSE) stop("lower.tail must be specified.")
  
  pald = ifelse(test= (q < mu), yes = p*exp((1-p)*(q-mu)/sigma),
                no = 1-(1-p)*exp(-(p)*(q-mu)/sigma))
  
  ifelse(test=lower.tail == TRUE,yes=return(pald),no=return(1-pald))
}

subset_mat <- function(X,start,j,m)  
{
  idx = seq(start,ncol(X),j)
  return(matrix(X[,idx],nrow=m))
}

#' @rdname qbild
#' @export
#### loglike + rho + aic + bic
mofit <- function(y, fixed, random, beta, alpha, varphi2, p, 
                      fixed_intercept, random_intercept, k, l, m, n)
{
  
  x  = array(0, dim = c(k,m,n))
  s  = array(0, dim = c(l,m,n))
  loglike = 0
  
  for (i in 1:n)
  {
    x[,,i] = t(subset_mat(fixed,i,n,m))   ##model matrix fixed + intercept
    s[,,i] = t(subset_mat(random,i,n,m)) ##model matrix random + intercept
    meani = subset_mat(fixed,i,n,m)%*%beta + subset_mat(random,i,n,m)%*%alpha[,i]
    
    for(j in 1:m)
    {
      if (y[j,i] == 1)
        loglike = loglike + log(1 - pald(0, meani[j,1], 1, p))
    else
        loglike = loglike + log(pald(0, meani[j,1], 1, p))
    }
  }
  
  x = t(matrix(x,k,m*n))
  s = t(matrix(s,l,m*n))
  
  xtx  = crossprod(x)
  xts  = crossprod(x,s)
  stx  = t(xts)
  sts  = crossprod(s)
  
  variance      = (1-2*p + 2*p^2)/((p*(1-p))^2);  
  invDstar      = (variance/varphi2)*diag(l);
  
  Hatcomp1  = rbind(cbind(xtx,xts),cbind(stx,sts+invDstar))
  Hatcomp2  = rbind(cbind(xtx,xts),cbind(stx,sts))
  
  rho = sum(diag(solve(Hatcomp1,Hatcomp2)))
  
  AIC = -2*loglike + 2*(rho+1);                          # AIC values are adjusted for mixed effect models
  BIC = -2*loglike + (l+1)*log(n) + (k)*log(m*n);        # BIC values are adjusted for mixed effect models
  
  return(c(loglike,AIC,BIC))
}
  
#' @rdname qbild
#' @export
"make.qbld" <- function(data,p=0.25,nsim=0,burn=0,varnames.fixed="missing",which="Blocked")
{
  if(missing(data))
    stop("Data must be provided.")
  
  attr(data,"burn") <- FALSE
  
  if(burn!=0)
  {
    burn = floor(burn*nsim)
    data[[1]] = data[[1]][burn:nsim,]
    data[[2]] = data[[2]][,,burn:nsim]
    data[[3]] = data[[3]][burn:nsim]
    attr(data,"burn") <- TRUE
  }
  attr(data,"which") <- which
  attr(data,"nsim") <- nsim
  attr(data,"varnames") <- varnames.fixed
  attr(data,"class") <- "qbld"
  attr(data,"quantile") <- p
  return(data)
}


"is.qbld" <- function (x) 
{
  if (inherits(x, "qbld")) 
    return(TRUE)
  return(FALSE)
}

"as.qbld" <- function (x, ...) 
  UseMethod("as.qbld")

"as.qbld.default" <- function (x, ...) 
  if (is.qbld(x)) x else make.qbld(x)

#' @rdname qbild
#' @export
"summary.qbld" <-
  function (data, quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975), ...) 
  {
    x <- as.qbld(data)
    statnames <- c("Mean", "SD", "MCSE", "ESS", "GR Diagnostic")
    
    varstats <- matrix(0,nrow = length(attr(data,"varnames")), ncol = length(statnames),
                       dimnames = list(attr(data,"varnames"),statnames)) # dimnames = list(varnames(x), statnames)
    ## RGA replaced with safespec0
    #sp0 <- function(x) spectrum0(x)$spec
    if (is.list(data)) {
      xmean <- round(c(colMeans(data[[1]]),mean(data[[3]])),2)
      xsd <- c(round(apply(data[[1]], 2, sd),2),round(sd(data[[3]]),2))
      xmcse <- c(round(mcse.mat(data[[1]]),3)[,2],round(mcse.mat(data[[3]]),3)[2])
      xess <- round(c(ess(data[[1]]), ess(data[[3]])),2)
      xsgr <- c(stable.GR((data[[1]]))$psrf,stable.GR((data[[3]]))$psrf)
      varquant <- round(rbind(t(apply(data[[1]], 2, quantile, quantiles)), t(apply(data[[3]], 2, quantile, quantiles))),3)
      rownames(varquant) <- attr(data,"varnames")
    }
    else {
      stop("Error in summary, Not a List!")
    }
    varstats[, 1] <- xmean
    varstats[, 2] <- xsd
    varstats[, 3] <- xmcse
    varstats[, 4] <- xess
    varstats[, 5] <- xsgr
    varstats <- drop(varstats)
    nsim <- attr(data,"nsim")
    burnin <- attr(data,"burn")
    wich <- attr(data,"which")
    quantile = attr(data,"quantile")
    multiess = c(multiESS(data[[1]]),multiESS(data[[3]]))

    out <- list(statistics = varstats, quantums = varquant, run=nsim, burn=burnin, 
                block=wich, quant=quantile, muless=multiess)
    
    class(out) <- "summary.qbld"
    return(out)
  }

#' @rdname qbild
#' @export
"print.summary.qbld" <-function (x, ...) 
  {
    cat("\n", "Quantile used = ",x$quant,"\n", sep = "")
    cat("\n", "No. of Iterations = ",x$run," samples\n", sep = "")
    cat("Type of Sampler =", x$block, "\n")
    cat("Burn-in Used? =", x$burn, "\n")
    cat("\n1. Statistics for each variable,\n")
    print(x$statistics, ...)
    cat("\n")
    cat("\n2. Quantiles for each variable,\n")
    print(x$quantums)
    cat("\n")
    cat("MultiESS value =", x$muless, "\n")
    invisible(x)
  }

#' @rdname qbild
#' @export
"plot.qbld" <- function (data, trace = TRUE, density = TRUE, 
                        auto.layout = TRUE, ask = dev.interactive(), ...) 
{
  data <- as.qbld(data)
  pars <- NULL
  on.exit(par(pars))
  vars <- attr(data,"varnames")
  nvars = length(vars)
  
  if (auto.layout) {
    mfrow <- set.mfrow(Nparms = nvars, nplots = trace + density)
    pars <- par(mfrow = mfrow)
  }
  
  for (i in 1:(nvars-1))
    {
      nam <- vars[i]
      beta = ts(data[[1]][,i])
      myplot(beta, l = nam, trace, density)
      
      if (i==1)
        pars <- c(pars, par(ask=ask))
  }
  
    nam = vars[nvars]
    myplot((data[[3]]), l = nam, trace, density)
  
}



"myplot" <- function (dat, l="covariate", trace, density) 
  {
        if(trace)
        {
          plot.ts(x=dat, main=paste("Trace of", l))
        #  abline(h=mn,col="red")
        }
  
        if(density)
        {
          plot(density(dat), main=paste("Density of", l))
         # abline(v=mn,col="blue")
          rug(dat,col="red")
        }
  }


"set.mfrow" <-function (Nparms = 1, nplots = 1) 
  {
    ## Set up dimensions of graphics window: 
    ## If only density plots OR trace plots are requested, dimensions are: 
    ##	1 x 1	if Nparms = 1 
    ##	1 X 2 	if Nparms = 2 
    ##	2 X 2 	if Nparms = 3 or 4 
    ##	3 X 2 	if Nparms = 5 or 6 or 10 - 12 
    ##	3 X 3 	if Nparms = 7 - 9 or >= 13 
    ## If both density plots AND trace plots are requested, dimensions are: 
    ##	1 x 2	if Nparms = 1 
    ##	2 X 2 	if Nparms = 2 
    ##	3 X 2 	if Nparms = 3, 4, 5, 6, 9 
    ##	4 x 2	if Nparms otherwise 
      if (nplots==1) {
        ## One plot per variable
        mfrow <- switch(min(Nparms,13),
                        c(1,1),
                        c(1,2),
                        c(2,2),
                        c(2,2),
                        c(3,2),
                        c(3,2),
                        c(3,3),
                        c(3,3),
                        c(3,3),
                        c(3,2),
                        c(3,2),
                        c(3,2),
                        c(3,3))
      }
      else {
        ## Two plot per variable
        ##
        mfrow <- switch(min(Nparms, 13),
                        c(1,2),
                        c(2,2),
                        c(3,2),
                        c(3,2),
                        c(3,2),
                        c(3,2),
                        c(4,2),
                        c(4,2),
                        c(3,2),
                        c(4,2),
                        c(4,2),
                        c(4,2),
                        c(4,2))
      }
    return(mfrow)
  }



