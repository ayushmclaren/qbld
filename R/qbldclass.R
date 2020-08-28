"make.qbld" <- function(data,p=0.25,nsim=0,burn=0,varnames.fixed="missing",which="Blocked") #make qbld object
{
  if(missing(data))
    stop("Data must be provided.")
  
  attr(data,"burn") <- FALSE
  attr(data,"nsim") <- nsim*(1-burn)
  if(burn!=0)
  {
    burn = floor(burn*nsim)
    data[[1]] = as.matrix(data[[1]][burn:nsim,])
    data[[2]] = data[[2]][,,burn:nsim,drop=FALSE]
    data[[3]] = as.matrix(data[[3]][burn:nsim])
    attr(data,"burn") <- TRUE
  }
  attr(data,"which") <- which
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


#' @title QBLD Summary Class
#' @name summary.qbld
#' @description Outputs a `summary.qbld' class object, and prints as described.
#' 
#' @details `qbld.summary' class summarizes the outputs of the model.qbld function. 
#' Markov Std Error (MCSE), Effective sample size (ESS) are calculated using mcmcse package. 
#' Gelman-Rubin diagnostic (R hat), and significance stars are indicated using Vats and Knudson et. al.
#' 
#' @usage \method{summary}{qbld}(object,quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975),epsilon = 0.05,...)
#' 
#' @param object : `qbld' class object
#' @param x : (for print.summary.qbld) `qbld.summary' class object
#' @param quantiles : Vector of quantiles for summary of the covariates, 
#' defaulted to \code{c(0.025, 0.25, 0.5, 0.75, 0.975)}
#' @param epsilon : epsilon value for calculating significance stars (see details), 0.05 by default.
#' @param ... : Other summary arguments
#' 
#' @return  summary.qbld produces following sets of summary statistics for each variable:
#' \itemize{
#' \item {\code{statistics:}} { Contains the mean, sd, markov std error, ess and Gelman-Rubin diagnostic}
#' \item {\code{quantiles:}}  { Contains quantile estimates for each variable}
#' \item {\code{nsim:}}       { No. of simulations run}
#' \item {\code{burn:}}       { Burn-in used or not}
#' \item {\code{which:}}      { Block, or Unblock version of sampler}
#' \item {\code{p:}}          { quantile for the AL distribution on the error term}
#' \item {\code{multiess:}}   { multiess value for the sample}
#' \item {\code{multigelman:}} { multivariate version of Gelman-Rubin}
#' }
#' 
#' @seealso \code{\link{plot.qbld}}, \code{\link{model.qbld}}  
#' 
#' Additional functions : \code{\link[mcmcse]{mcse.mat}}, \code{\link[mcmcse]{ess}}, \code{\link[mcmcse]{multiESS}}, 
#' \code{\link[stableGR]{stable.GR}}, \code{\link[stableGR]{target.psrf}}
#' 
#' @references 
#' Vats, Dootika and Christina Knudson. “Revisiting the Gelman-Rubin Diagnostic.” arXiv 
#' 
#' James M. Flegal, John Hughes, Dootika Vats, and Ning Dai. (2020). mcmcse: Monte Carlo Standard Errors for MCMC. R
#' package version 1.4-1. Riverside, CA, Denver, CO, Coventry, UK, and Minneapolis, MN.
#' 
#' Christina Knudson and Dootika Vats (2020). stableGR: A Stable Gelman-Rubin Diagnostic for Markov Chain Monte Carlo. R
#' package version 1.0.
#' 


#' @rdname summary.qbld
#' @export
"summary.qbld" <- function (object,quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975),epsilon = 0.05,...) 
  {
    object <- as.qbld(object)
    nsim <- attr(object,"nsim")
    
    if (is.list(object)) {
      xmean <- round(c(colMeans(object[[1]]),mean(object[[3]])),3) #Mean
      xsd <- c(round(apply(object[[1]], 2, sd),2),round(sd(object[[3]]),3)) #SD
      xmcse <- c(round(mcmcse::mcse.mat(object[[1]]),3)[,2],round(mcmcse::mcse.mat(object[[3]]),3)[2]) #MCSE
      xess <- floor(c(mcmcse::ess(object[[1]]), mcmcse::ess(object[[3]]))) #ESS
      xsgr <- sqrt((nsim-1)/nsim +  1/xess) #Rhat via ESS
      varquant <- round(rbind(t(apply(object[[1]], 2, quantile, quantiles)), t(apply(object[[3]], 2, quantile, quantiles))),3)
      rownames(varquant) <- attr(object,"varnames")
    }
    else {
      stop("Error in summary, Not a List!")
    }
    targetpsrfuni = stableGR::target.psrf(m = 1, p = 1, epsilon = epsilon)$psrf #targte psrf
    varstats <- data.frame(Mean=xmean, SD=xsd, MCSE=xmcse, ESS=xess, GelmanRubin=xsgr, 
                           Signif.= ifelse( test = xsgr < targetpsrfuni, yes = '*', no = ''),
                           row.names=attr(object,"varnames")) #to add the signif stars
    colnames(varstats) <- c("Mean", "SD", "MCSE", "ESS", "Gelman-Rubin"," ") #create matrix to output
    
    varstats <- drop(varstats)
    burnin <- attr(object,"burn")
    which <- attr(object,"which")
    quantile = attr(object,"quantile")
    multiess = mcmcse::multiESS(cbind(object[[1]],object[[3]]))
    multisgr = sqrt((nsim-1)/nsim +  1/multiess)
    targetpsrfmulti = stableGR::target.psrf(m = 1, p = nrow(varstats), epsilon = epsilon)$psrf #to add the signif stars
    stars = multisgr < targetpsrfmulti
    out <- list(statistics = varstats, quantiles = varquant, nsim=nsim, burn=burnin, 
                which=which, p=quantile, multiess=multiess, multigelman = multisgr, foo = stars)
    
    class(out) <- "summary.qbld" 
    return(out)
  }


#' @rdname summary.qbld
#'@export
"print.summary.qbld" <-function (x, ...) 
{
  cat("\n", "Quantile used = ",x$p,"\n", sep = "")
  cat("\n", "No. of Iterations = ",x$nsim," samples\n", sep = "")
  cat("Type of Sampler =", x$which, "\n")
  cat("Burn-in Used? =", x$burn, "\n")
  cat("\n1. Statistics for each variable,\n")
  print(x$statistics, ...)
  cat("\n")
  cat("MultiESS value =", x$multiess, "\n")
  cat("Multi Gelman-Rubin =", x$multigelman)
  if(x$foo == TRUE) cat(" ***")
  cat("\nNote : * indicates enough samples for the covariate\n")
  cat("       *** indicates enough samples for the whole sampler.\n")
  cat("\n2. Quantiles for each variable,\n")
  print(x$quantiles)
  cat("\n")
  invisible(x)
}


#' @title Plot QBLD
#' @name plot.qbld
#' @description Plots `qbld' class object.
#' 
#' @usage \method{plot}{qbld}(x,trace = TRUE,density = TRUE,auto.layout = TRUE,ask = dev.interactive(),...) 
#' 
#' @param x : `qbld' class object to plot.
#' @param trace : Whether or not to plot trace plots for covariates, TRUE by default
#' @param density : Whether or not to plot density for covariates, TRUE by default.
#' @param auto.layout : Auto set layout or not, TRUE as default. Plots according to the local settings if false.
#' @param ask : For Interactive plots
#' @param ... : Other plot arguments
#' 
#' @return Plots as specified.
#' 
#' @seealso \code{\link{summary.qbld}}, \code{\link{model.qbld}} 


#' @export
"plot.qbld" <- function (x,trace = TRUE,density = TRUE,auto.layout = TRUE,ask = dev.interactive(),...) 
  #trace and density can be chosen whether or not to output
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

