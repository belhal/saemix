% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SaemixData.R
\name{transformContCov}
\alias{transformContCov}
\alias{transform.SaemixData}
\title{Transform covariates}
\usage{
transformContCov(object, covariate, transformation = function(x) x,
  centering = "median", verbose = FALSE)
}
\arguments{
\item{object}{saemixData object}

\item{covariate}{name of the covariate}

\item{transformation}{transformation function. Defaults to no transformation}

\item{centering}{string, giving the value used to center the covariate; can be "mean" or "median", in which case this value will be computed from the data, 'none' or 0 for no centering, or a value given by the user. Defaults to the median value over the dataset.}

\item{verbose}{a boolean, prints messages during the execution of the function if TRUE. Defaults to FALSE.}
}
\value{
an object of class \code{"\linkS4class{SaemixData}"}
}
\description{
Transform and/or center continuous covariates
}
\examples{
# TODO
}
\keyword{data}
