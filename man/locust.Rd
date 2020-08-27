\name{locust}
\alias{locust}
\docType{data}
\title{Dataset}
\description{This data set was presented by MacDonald and Raubenheimer (1995) 
and analyze the effect of hunger on locomotory behaviour of 24 locust (\emph{Locusta migratoria}) 
observed at 161 time points. The subjects were divided in two treatment groups ("fed" and "not fed"), 
and within each of the two groups, the subjects were alternatively "male" and "female". 
For the purpose of this analysis the categories of the response variable 
were "moving" and "not moving". During the observation period, the behavior of each of the subjects was 
registered every thirty seconds.}
\usage{data(locust)}
\format{
  A data frame with 3864 observations on the following 7 variables.
  \describe{
    \item{\code{id}}{a numeric vector that identifies de number of the individual profile.}
    \item{\code{move}}{a numeric vector representing the response variable.}
    \item{\code{sex}}{a factor with levels \code{1} for "male" and \code{0} for "female". }
    \item{\code{time}}{a numeric vector that identifies de number of the time points observed. 
    The \code{time} vector considered was obtained dividing (1:161) by 120 (number of observed periods in 1 hour).}
    \item{\code{feed}}{a factor with levels \code{0} "no" and \code{1} "yes".}
  }
}
\details{The response variable, \code{move} is the binary type coded as \code{1} for "moving" and \code{0} for "not moving". 
The \code{sex} covariate was coded as \code{1} for "male" and \code{0} for "female". The \code{feed} covariate indicating the treatment group,  
was coded as \code{1} for "fed" and \code{0} for "not fed". Azzalini and Chiogna (1997) also have analyze this 
data set using their \code{S-plus} package \code{rm.tools}. }
\source{MacDonald, I. and Raubenheimer, D. (1995). Hidden Markov models and animal behaviour. 
\emph{Biometrical Journal}, 37, 701-712}
\references{Azzalini, A. and Chiogna, M. (1997). S-Plus Tools for the Analysis of Repeated Measures Data. 
\emph{Computational Statistics}, 12, 53-66}