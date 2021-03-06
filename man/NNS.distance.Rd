% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NNS_Distance.R
\name{NNS.distance}
\alias{NNS.distance}
\title{NNS Distance}
\usage{
NNS.distance(rpm, rpm_class, dist.estimate, type, k, n)
}
\arguments{
\item{rpm}{REGRESSION.POINT.MATRIX from \link{NNS.reg}}

\item{rpm_class}{integer \code{rpm}.}

\item{dist.estimate}{Vector to generate distances from.}

\item{type}{"L1", "L2", "DTW" or "FACTOR"}

\item{k}{\code{n.best} from \link{NNS.reg}}

\item{n}{number of observations.}
}
\value{
Returns sum of weighted distances.
}
\description{
Internal kernel function for NNS multivariate regression \link{NNS.reg} parallel instances.
}
