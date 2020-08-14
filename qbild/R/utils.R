#' @title QBLD sampler!
#'
#' @description Runs the QBLD sampler and outputs a qbld class object which consists of a list of markov chains
#' for Beta(the fixed effects estimate), Alpha(the random effects estimate), and Varphi2 (as per the model),
#' of which Beta and Varphi2 are of interest.  
#'
#' @details For a detailed information on the sampler, please check README.pdf 
#' data are contained in a data.frame. Each element of the data argument must be identifiable by
#' a name. The simplest situation occurs when all subjects are observed at the same time points. The
#' response variable represent the individual profiles of each subject, it is expected a variable in the
#' data.frame that identifies the correspondence of each component of the response variable to the
#' subject that it belongs, by default is named id variable. Hence NA values are not valid.
#' 
#' `qbld` object contains markov chains and sampler run information as attributes , and is compatible 
#' with S3 methods like summary,plot. make.qbld function can be used to convert a similar
#' type-object to `qbld` class. Check qbild/qbild2.pdf and the code to check for different additional options
#' in summary and plots.
#'
#' @name model.qbld
#' @aliases qbld
#' @param fixed_formula : a description of the model to be fitted of the form 
#' response~fixed effects predictors i.e Xi in the model.
#' @param data : data frame, NAs not allowed and should throw errors, factor variables are auto-converted, 
#' find airpollution.rda and locust.rda built into the package.
#' @param id : variable name in the dataset, that specifies individual profile. By default, id="id" and
#' data is expected to conatin an id variable.  
#' @param random_formula : a description of the model to be fitted of the form 
#' response~random effects predictors i.e Si in the model. This defaults to Si being only an intercept.
#' @param p : quantile for the AL distribution on the error term, p=0.25 by default.
#' @param nsim : number of simultions to run the sampler.
#' @param b0,B0 : Prior model parameters for Beta. See README for details.
#' @param c1,d1 : Prior model parameters for Varphi2. See README for details.
#' @param method : Choose between the "block" vs "Unblock" sampler, regex is used to counter 
#' the effects of small/large cases.
#' @param burn : number between (0,1). Burn-in values are discarded and not used for summary calculations
#' @param summarize : Outputs a summary table (same as summary(output)), in addition also prints Model fit
#' AIC/BIC/Log-likelihood values (same as calling the 'mofit' function with appropriate inputs). False by default.
#' @param verbose : False by default. Spits out progress reports while the sampler is running.
#'
#' @return
#' \itemize{
#' \item {\code{model.qbld}} {Returns `qbld` class object.}
#' }
#' 
#' @examples
#' data = load("airpollution")
#' data2 = load("locust")
#' output <- (fixed_formula = wheeze~smoking+I(age^2)-1, data = data, id="id", 
#'             random_formula = ~-., p=0.25, 
#'            nsim=1000, method="block", burn=0, 
#'            summarize=FALSE, verbose=FALSE)
#' summary(output)
#' plot(output) 
#' output2 <- (fixed_formula = move~sex+I(time^2), data = data2, id="id", 
#'             random_formula = ~-., p=0.50, 
#'            nsim=1000, method="block", burn=0, 
#'            summarize=TRUE, verbose=FALSE)
#'
#' @references
#'
#' Rahman and Vossmyer 2019, https://home.iitk.ac.in/~marshad/RahmanVossmeyer2019.pdf?
#'
#'


#' @rdname model.qbld
#' @export

model.qbld <- function(fixed_formula, data, id="id", random_formula = ~-., p=0.25, 
                       nsim, b0, B0, c1, d1, method="block", burn=0, 
                       summarize=FALSE, verbose=FALSE)
{
  
  if(missing(nsim)||nsim==0||nsim%%1!=0)
    stop("n must be specified correctly.")
  
  if(missing(data)||(sum(is.na(data)) > 0)|| !is.data.frame(data))
    stop("Either the data isn't provided as a data frame or contains NAs.")
  
  if(is.null(names(data)))
    stop("objects in data.frame must have a name")
  
  if(burn<0 || burn>1)
    stop("burn should be a number in [0,1].")
  
  if(missing(fixed_formula))
    stop("Please provide valid formulas for fixed effects.")
  
  must_convert<-sapply(data,is.factor)       # logical vector telling if a variable needs to be displayed as numeric
  if(sum(must_convert))
  {
    df2 <-sapply(data[,must_convert,drop=FALSE],unclass) # data.frame of all categorical variables now displayed as numeric
    #colnames(df2) <-names(data)[must_convert]
    data<-cbind(data[,!must_convert],df2)  #converted and joined back appropriately 
  }
   
  if(!missing(id)) {indvs<-as.vector(data[[id]])} #to creeate individual profile.
  if (missing(id)){ if (all(is.na(match(names(data), "id")))) stop ("id must be defined")
    else indvs<-as.vector(data$id)}
  
  
  m = sum(indvs==1) #no. of rows for my data by counting no. of entries for each individual 
  
  fixed_formula = formula(fixed_formula) # processing fixed formula
  expr1 = terms(fixed_formula, data=data)
  expr <- attr(expr1, "variables")
  var.names <- attr(expr1, "term.labels")
  
  if(any(is.na(match(all.vars(expr), names(data)))))
    stop("Variables in formula not contained in the data.frame")
  
  response <- as.character(expr[[2]])
  if(attr(expr1,"response")==0 || is.null(response)) #response for fixed is needed
    stop("Please provide response variable in fixed_formula correctly.")
  y = matrix(as.matrix(data[response]),nrow=m) 
  n = ncol(y) # second dimension
  
  fixed_intercept = attr(expr1,"intercept") 
  if(fixed_intercept == 1)
    var.names = c("Intercept",var.names)
  
  fixed = model.matrix(fixed_formula,data=data) #model matrix
  if(ncol(fixed)==0)
    stop("Invalid/Empty fixed_formula resulting in empty model matrix.")
  fixed = matrix(fixed,nrow=m)
  k = floor(ncol(fixed)/n) #no. of covariates
  
  random_formula = formula(random_formula) #processing random formula 
  expr2 = terms(random_formula, data=data)
  random_intercept = attr(expr2,"intercept")
  random = model.matrix(random_formula,data=data) #random model matrix
  if(ncol(random)==0)
    stop("Invalid/Empty random_formula resulting in empty model matrix.")
  random = matrix(random,nrow=m)
  l = floor(ncol(random)/n) #no. of covariates
  
  
  if(missing(b0))
    b0 = rep(0,k) #start from 0 if not specified
  if(length(b0)!=k)
    stop("b0 dimensions are not compatible with fixed_formula.")
  if(missing(B0))
    B0 = diag(k) #identity matrix
  if(nrow(B0)!=k ||ncol(B0)!=k)
    stop("B0 dimensions are not compatible with fixed_formula.")
  
  if(missing(c1))
    c1 = 9
  if(missing(d1))
    d1 = 10
  
  
  var.names = c(var.names,"Varphi2") #Take names form the data frame and add varphi2
  
  if(regexec("block",method,ignore.case=TRUE)[[1]][1]==1) #blocked
    out = (qbldf(nsim, p, y, fixed, random, b0, B0, c1, d1, m, n, k, l, verbose)) #check block.cpp
  
  if(regexec("unblock",method,ignore.case=TRUE)[[1]][1]==1) #unblocked
    out = (qbldunblock(nsim, p, y, fixed, random, b0, B0, c1, d1, m, n, k, l, verbose)) #check unblock.cpp
  
  
  if(is.null(out))
    stop("Method is either 'block' or 'unblock'.")
  
  out = make.qbld(out, p, nsim, burn, var.names, method)  #make qbld object
  
  if(summarize==TRUE)
  {
    beta  = matrix(colMeans(out[[1]]),ncol=1) #for model fit
    alpha = rowMeans(out[[2]],dims=2)
    varphi2 = mean(out[[3]])
    
    model_fit = mofit(y, fixed, random, beta, alpha, varphi2, p, 
                      fixed_intercept, random_intercept, k, l, m, n) #AIC/BIC
    
    print(summary(out))
    
    cat("\n3. Model Selection Criterion\n")
    cat("Log likelihood =", model_fit[1], "\n")
    cat("AIC =", model_fit[2], "\n")
    cat("BIC =", model_fit[3], "\n")
  }
  
  return(out)
}


#### AL distribution CDF
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
"make.qbld" <- function(data,p=0.25,nsim=0,burn=0,varnames.fixed="missing",which="Blocked") #make qbld object
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
  function (object, quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975), ...) 
  {
    object <- as.qbld(object)
    nsim <- attr(object,"nsim")
    statnames <- c("Mean", "SD", "MCSE", "ESS", "Gelman-Rubin","Signif.") #create matrix to output
    
    varstats <- matrix(0,nrow = length(attr(object,"varnames")), ncol = length(statnames),
                       dimnames = list(attr(object,"varnames"),statnames)) # dimnames = list(varnames(x), statnames)
  
    
    if (is.list(object)) {
      xmean <- round(c(colMeans(object[[1]]),mean(object[[3]])),2) #Mean
      xsd <- c(round(apply(object[[1]], 2, sd),2),round(sd(object[[3]]),2)) #SD
      xmcse <- c(round(mcse.mat(object[[1]]),3)[,2],round(mcse.mat(object[[3]]),3)[2]) #MCSE
      xess <- floor(c(ess(object[[1]]), ess(object[[3]]))) #ESS
      xsgr <- sqrt((nsim-1/nsim) +  1/xess) #Rhat via ESS
      varquant <- round(rbind(t(apply(object[[1]], 2, quantile, quantiles)), t(apply(object[[3]], 2, quantile, quantiles))),3)
      rownames(varquant) <- attr(object,"varnames")
    }
    else {
      stop("Error in summary, Not a List!")
    }
    targetpsrfuni = target.psrf(m = 1, p = 1, epsilon = 0.10)$psrf #targte psrf
    varstats[, 1] <- xmean 
    varstats[, 2] <- xsd
    varstats[, 3] <- xmcse
    varstats[, 4] <- xess
    varstats[, 5] <- xsgr
    varstats[, 6] <- ifelse( test = xsgr < targetpsrfuni, yes = 1, no = 0) #to add the signif stars
    varstats <- drop(varstats)
    burnin <- attr(object,"burn")
    wich <- attr(object,"which")
    quantile = attr(object,"quantile")
    multiess = multiESS(cbind(object[[1]],object[[3]]))
    multisgr = sqrt((nsim-1)/nsim +  1/multiess)
    targetpsrfmulti = target.psrf(m = 1, p = 1, epsilon = 0.10)$psrf #to add the signif stars
    stars = multisgr < targetpsrfmulti
    out <- list(statistics = varstats, quantums = varquant, run=nsim, burn=burnin, 
                block=wich, quant=quantile, muless=multiess, musgr = multisgr, star = stars)
    
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
    cat("MultiESS value =", x$muless, "\n")
    cat("Multi Gelman-Rubin Diagnostic =", x$musgr)
    if(x$star == TRUE) cat(" *")
    else cat(" #")
    cat("\nNote : 1(0) indicates samples enough (not) for the covariate,*(#) indicates samples enough (not) for the whole sampler.\n")
    cat("\n2. Quantiles for each variable,\n")
    print(x$quantums)
    cat("\n")
    invisible(x)
  }

#' @rdname qbild
#' @export
"plot.qbld" <- function (x, trace = TRUE, density = TRUE, 
                        auto.layout = TRUE, ask = dev.interactive(), ...) #trace and density can be chosen whether or not to output
{
  x <- as.qbld(x)
  pars <- NULL
  on.exit(par(pars))
  vars <- attr(x,"varnames")
  nvars = length(vars)
  
  if (auto.layout) {
    mfrow <- set.mfrow(Nparms = nvars, nplots = trace + density)
    pars <- par(mfrow = mfrow)
  }
  
  for (i in 1:(nvars-1))
    {
      nam <- vars[i]
      beta = ts(x[[1]][,i])
      myplot(beta, l = nam, trace, density)
      
      if (i==1)
        pars <- c(pars, par(ask=ask))
  }
  
    nam = vars[nvars]
    myplot((x[[3]]), l = nam, trace, density)
  
}



"myplot" <- function (dat, l="covariate", trace, density) 
  {
        if(trace)
        {
          plot.ts(x=dat, main=paste("Trace of", l),ylab="Trace")
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



