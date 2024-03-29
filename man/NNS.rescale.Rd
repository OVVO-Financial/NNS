% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Central_tendencies.R
\name{NNS.rescale}
\alias{NNS.rescale}
\title{NNS rescale}
\usage{
NNS.rescale(x, a, b)
}
\arguments{
\item{x}{vector of data.}

\item{a}{numeric; lower limit.}

\item{b}{numeric; upper limit.}
}
\value{
Returns a rescaled distribution within provided limits.
}
\description{
Rescale min-max scaling output between two numbers.
}
\examples{
\dontrun{
set.seed(123)
x <- rnorm(100)
NNS.rescale(x, 5, 10)
}
}
\author{
Fred Viole, OVVO Financial Systems
}
