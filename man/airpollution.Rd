\name{airpollution}
\alias{airpollution}
\docType{data}
\title{Dataset}
\description{This example is a subset of data from Six Cities study, a longitudinal study of the health effects of 
air pollution (Ware, J. H. et al., 1984).}
\usage{data(airpollution)}
\format{
  A data frame with 128 observations on the following 5 variables.
  \describe{
    \item{\code{id}}{identifies de number of the individual profile. This vector contains observations of 537 individual profiles.}
    \item{\code{wheeze}}{a numeric vector that identify the wheezing status (1="yes", 0="no") of a child at each occasion.}
    \item{\code{age}}{a numeric vector corresponding to the \code{age} in years.}
    \item{\code{smoking}}{a factor that identify if the mother smoke (1="smoke", 0="no smoke").}
    \item{\code{counts}}{a numeric vector corresponding to the replications of each individual profile.}
  }
}
\details{The data set presented by Fitzmaurice and Laird (1993) contains complete records
on 537 children from Steubnville, Ohio, each woman was examined annually at ages 7 through 10. 
The repeated binary response is the wheezing status (1="yes", 0="no") of a child at each occasion. 
Although mother's smoking status could vary with time, it was determined in the first interview and 
was treated as a time-independent covariate. Maternal smoking was categorized as 1 if the mother
smoked regularly and 0 otherwise.}
\source{ Fitzmaurice, G. M. and Laird, N. M. (1993). 
A Likelihood-Based Method for analyzing Longitudinal Binary Response. 
\emph{Biometrika}, 80, 141-51.}

\references{  Ware, J. H., Dockery, D. W., Spiro,  A. III, Speizer, F. E. and Ferris, B. G., Jr. (1984). 
Passive smoking, gas cooking and respiratory health in children 
living in six cities. \emph{Am. Rev. Respir. dis.}, 129, 366-74.
}
