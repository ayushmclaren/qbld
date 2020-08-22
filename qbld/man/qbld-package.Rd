\name{qbld-package}
\alias{qbld-package}
\alias{qbld}
\docType{package}
\title{
\packageTitle{qbld}
}
\description{
\packageDescription{qbld}
}

\details{
\tabular{ll}{
Package: \tab qbld\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2020-08-17\cr
License: \tab GPL (>= 3)\cr
}

The package contains the following functions: 

\itemize{
\item \code{\link{model.qbld}} {: Runs the QBLD sampler as in Rahman and Vossmeyer(2019) and outputs a `qbld` class object.}
\item \code{\link{summary.qbld}} {: S3 method that summarizes the outputs of the model.qbld function.}
\item \code{\link{plot.qbld}} {: S3 method that plots `qbld` class object. }
\item \code{\link{aldmix}} {: Cumulative density, Probability distribution function, Quantile function and Random generation for the Asymmetric Laplace distribution.}
\item \code{\link{gig}} {: Probability distribution function, Quantile Random generation for the Generalised Inverse Gaussian.}
\item \code{\link{airpollution}}, \code{\link{locust}} {: In-built datasets} 
}

}

\author{
\packageAuthor{qbld}

Maintainer: \packageMaintainer{qbld}
}
\references{
Rahman, Mohammad & Vossmeyer, Angela. (2018). 
Estimation and Applications of Quantile Regression for Binary Longitudinal Data. 
Advances in Econometrics. 40.   
URL: \url{https://home.iitk.ac.in/~marshad/RahmanVossmeyer2019.pdf?}

Vats, Dootika and Christina Knudson. “Revisiting the Gelman-Rubin Diagnostic.” arXiv: Computation (2018): n. pag. 

Keming Yu & Jin Zhang (2005) A Three-Parameter Asymmetric Laplace Distribution
and Its Extension, Communications in Statistics - Theory and Methods.

Kobayashi, Genya. (2011). Gibbs Sampling Methods for Bayesian Quantile Regression.
J Stat Comput Simul.

Devroye, L. Random variate generation for the generalized inverse Gaussian distribution. 
Stat Comput 24, 239–246 (2014).

Wolfgang Hörmann and Josef Leydold (2013). 
Generating generalized inverse Gaussian random variates, Statistics and Computing (to appear).

J. S. Dagpunar (1989). An easily implemented generalised inverse Gaussian generator, 
Comm. Statist. B – Simulation Comput. 18, 703–710.
}

\examples{

# Dataset
data(airpollution)

# output will be a qbld class object, summarize is TRUE
output <- model.qbld(fixed_formula = wheeze~smoking+I(age^2)-1, data = airpollution, id="id", 
            random_formula = ~1, p=0.25, 
           nsim=1000, method="block", burn=0, 
           summarize=TRUE, verbose=FALSE)
           
# plots           
plot(output)

# GIG sampler
rgig(n = 1, lambda = 0.5, a = 1, b = 2)

# ALD sampler
raldmix(n = 10, mu = 5, sigma = 10, p = 0.5)
}

\keyword{package}
\keyword{qbld}

