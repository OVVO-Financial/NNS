% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Central_tendencies.R
\name{NNS.mode}
\alias{NNS.mode}
\title{NNS mode}
\usage{
NNS.mode(x, discrete = FALSE, multi = TRUE)
}
\arguments{
\item{x}{vector of data.}

\item{discrete}{logical; \code{FALSE} (default) for discrete distributions.}

\item{multi}{logical; \code{TRUE} (default) returns multiple mode values.}
}
\value{
Returns a numeric value representing the mode of the distribution.
}
\description{
Mode of a distribution, either continuous or discrete.
}
\examples{
\dontrun{
set.seed(123)
x <- rnorm(100)
NNS.mode(x)
}
}
\author{
Fred Viole, OVVO Financial Systems
}
