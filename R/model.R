#' @title QBLD Sampler
#'
#' @description Runs the QBLD sampler as in Rahman and Vossmeyer(2019) and outputs a `qbld' class object which 
#'  consists of Markov chains for Beta(the fixed effects estimate), Alpha(the random effects estimate),
#'  and Varphi2 (as per the model), of which Beta and Varphi2 are of interest.  
#'
#' @details For a detailed information on the sampler, please check the vignette.
#' Data are contained in a data.frame. Each element of the data argument must be identifiable by
#' a name. The simplest situation occurs when all subjects are observed at the same time points. The
#' id variable represent the individual profiles of each subject, it is expected a variable in the
#' data.frame that identifies the correspondence of each component of the response variable to the
#' subject that it belongs, by default is named id variable. Hence NA values are not valid.
#' For very low \eqn{(<=0.025)} or very high \eqn{(>=0.970)} values of \eqn{p}, sampler forces to unblock version to avoid errors.
#' Block version in this case may lead to machine tolerance issues.
#' 
#' `qbld' object contains markov chains and sampler run information as attributes , and is compatible 
#' with S3 methods like summary,plot. make.qbld function can be used to convert a similar
#' type-object to `qbld' class.
#'
#' @name model.qbld
#' @usage model.qbld(fixed_formula, data, id = "id", random_formula = ~1, p = 0.25, 
#'                   b0 = 0, B0 = 1, c1 = 9, d1 = 10, method = c("block","unblock"), 
#'                   nsim, burn = 0, summarize = FALSE, verbose = FALSE)
#' @param fixed_formula : a description of the model to be fitted of the form 
#' response~fixed effects predictors i.e \eqn{Xi} in the model. See vignette for more information.
#' @param data : data frame, NAs not allowed and should throw errors, factor variables are auto-converted, 
#' find airpollution.rda and locust.rda built into the package.
#' @param id : variable name in the dataset that specifies individual profile. By default, \code{id="id"} and
#' data is expected to contain an id variable. This is omitted while modelling.
#' @param random_formula : a description of the model to be fitted of the form 
#' response~random effects predictors i.e \eqn{Si} in the model. This defaults to \eqn{Si} being only an intercept.
#' See vignette for more information.
#' @param p : quantile for the AL distribution on the error term, \eqn{p=0.25} by default. For very low \eqn{(<=0.025)} or
#' very high \eqn{(>=0.970)} values of p, sampler forces to unblock version to avoid errors.
#' @param b0,B0 : Prior model parameters for Beta. These are defaulted to 0 vector, and Identity matrix.
#' @param c1,d1 : Prior model parameters for Varphi2. These are defaulted to 9,10 (arbitrary) respectively.
#' @param method : Choose between the "Block" vs "Unblock" sampler, Block is slower but produces lower correlation.
#' @param nsim : number of simultions to run the sampler.
#' @param burn : Burn in percentage, number between (0,1). Burn-in values are discarded and not used for summary calculations.
#' @param summarize : Outputs a summary table (same as \code{summary(output)}), in addition also prints Model fit
#' AIC/BIC/Log-likelihood values. False by default.
#' @param verbose : False by default. Spits out progress reports while the sampler is running.
#'
#' @return Returns `qbld' class object. `qbld' class contains the following :
#' \itemize{
#' \item {\code{Beta:}} { Matrix of MCMC samples of fixed-effects parameters.}
#' \item {\code{Alpha:}} { 3-dimensional matrix of MCMC samples of random-effects parameters.}
#' \item {\code{Varphi2:}} { Matrix of MCMC samples for varphi2.}
#' \item {\code{nsim:}} { Attribute; No. of simulations of chain run.}
#' \item {\code{burn:}} { Attribute; Whether or not burn-in used.}
#' \item {\code{which:}} {Attribute; "block" or "unblock" sampler used}
#' }
#' 
#' @examples
#' 
#' data(airpollution)
#' 
#' output <- model.qbld(fixed_formula = wheeze~smoking+I(age^2), data = airpollution, id="id", 
#'                      random_formula = ~1, p=0.25, nsim=1000, method="block", burn=0, 
#'                      summarize=TRUE, verbose=FALSE)
#'            
#' plot(output)
#'
#'
#'
#' @references
#' Rahman, Mohammad Arshad and Angela Vossmeyer, “Estimation and Applications of 
#' Quantile Regression for Binary Longitudinal Data,” 
#' Advances in Econometrics, 40B, 157-191, 2019. 
#'
#' 
#' 
#' @seealso A qbld object may be summarized by the summary function and visualized with the plot function.  
#' 
#' \code{\link{summary.qbld}}, \code{\link{plot.qbld}}  
#' 
#' Datasets : \code{\link{airpollution}}, \code{\link{locust}}


#' @rdname model.qbld
#' @export

model.qbld <- function(fixed_formula, data, id = "id", random_formula = ~1, p = 0.25, 
                       b0 = 0, B0 = 1, c1 = 9, d1 = 10, method = c("block","unblock"), 
                       nsim, burn = 0, summarize = FALSE, verbose = FALSE)
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
  
  if(p>=1 || p<=0)
    stop("p shpuld be between (0,1)")
  
  
  if(!missing(id)) {indvs<-as.vector(data[[id]])} #to creeate individual profile.
  if (missing(id)){ if (all(is.na(match(names(data), "id")))) stop ("id must be defined")
    else indvs<-as.vector(data$id)}
  
  
  data = data[ , ! colnames(data) %in% id ]  # remove the "id" variable 
  must_convert<-sapply(data,is.factor)       # logical vector telling if a variable needs to be displayed as numeric
  if(sum(must_convert))
  {
    df2 <-sapply(data[,must_convert,drop=FALSE],unclass) # data.frame of all categorical variables now displayed as numeric
    #colnames(df2) <-names(data)[must_convert]
    data<-cbind(data[,!must_convert],df2)  #converted and joined back appropriately 
  }
   
  
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
  var.names = attr(fixed,"dimnames")[[2]]
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
  
  var.names = c(var.names,"Varphi2") #Take names form the data frame and add varphi2
  
  flag = 0
  if(p <= 0.025 || p >= 0.970)
      flag = 1   #for these values call only unblock
  
  run = FALSE #whether or not this has run
  if(flag || regexec("unblock",method,ignore.case=TRUE)[[1]][1]==1) #unblocked
    {
        out = (qbldunblock(nsim, p, y, fixed, random, b0, B0, c1, d1, m, n, k, l, verbose)) #check unblock.cpp
        run = TRUE
        method = "unblock"
    }
  
  if(!run && regexec("block",method,ignore.case=TRUE)[[1]][1]==1) #blocked
  {
    out = (qbldf(nsim, p, y, fixed, random, b0, B0, c1, d1, m, n, k, l, verbose)) #check block.cpp
    method="block"
    run=TRUE
  }
  
  if(run==FALSE)
    stop("Sampler did not run, please check inputs and try again.")
  
  if(is.null(out))
    stop("Method is either 'block' or 'unblock'.")
  
  out = make.qbld(out, p, nsim, burn, var.names, method)  #make qbld object
  
  if(summarize==TRUE)
  {
    beta  = matrix(colMeans(out[[1]]),ncol=1) #for model fit
    alpha = rowMeans(out[[2]],dims=2)
    varphi2 = mean(out[[3]])
    
    print(summary(out))
    
    tryCatch({
      model_fit = mofit(y, fixed, random, beta, alpha, varphi2, p, 
                        fixed_intercept, random_intercept, k, l, m, n) #AIC/BIC
      
      cat("\n3. Model Selection Criterion\n")
      cat("Log likelihood =", model_fit[1], "\n")
      cat("AIC =", model_fit[2], "\n")
      cat("BIC =", model_fit[3], "\n")
    }, error=function(e){cat("ERROR in model-fit :",conditionMessage(e), "\n")})
  }
  
  return(out)
}


subset_mat <- function(X,start,j,m)  
{
  idx = seq(start,ncol(X),j)
  return(matrix(X[,idx],nrow=m))
}


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
        loglike = loglike + log(1 - paldmix(0, meani[j,1], 1, p))
    else
        loglike = loglike + log(paldmix(0, meani[j,1], 1, p))
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


## usethis namespace: start
#' @useDynLib qbld, .registration = TRUE
## usethis namespace: end
NULL

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL
qbldf <- function(nsim, p, y, datax, datas, b0, B0, c1, d1, m, n, k, l, verbose) {
  .Call(`_qbld_qbldf`, nsim, p, y, datax, datas, b0, B0, c1, d1, m, n, k, l, verbose)
}


## usethis namespace: start
#' @useDynLib qbld, .registration = TRUE
## usethis namespace: end
NULL

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL
qbldunblock <- function(nsim, p, y, datax, datas, b0, B0, c1, d1, m, n, k, l, verbose) {
  .Call(`_qbld_qbldunblock`, nsim, p, y, datax, datas, b0, B0, c1, d1, m, n, k, l, verbose)
}




