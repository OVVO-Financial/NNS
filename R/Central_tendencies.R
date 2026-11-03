#' NNS mode
#'
#' Mode of a distribution, either continuous or discrete.
#'
#' @param x vector of data.
#' @param discrete logical; \code{FALSE} (default) for continuous distributions, \code{TRUE} for discrete distributions.
#' @param multi logical; \code{TRUE} (default) returns multiple mode values.  \code{FALSE} returns the highest peak; among equally high peaks, the one nearest \link{NNS.gravity} of the data (then the lower value), so the result is always one of the modes.
#' @return Returns the mode of the distribution: a numeric vector of all modes when \code{multi = TRUE} (default), otherwise a single numeric value.
#' @details The continuous mode finds the peaks of a smoothed histogram with bin width IQR / sqrt(n), peaks in the first and last bins included, and centres each peak on the mean of the observations within three bins of it, so a peak on a single value returns that value.  When every gap between consecutive distinct values is wider than a bin, the data are treated as discrete-valued: the modes are the values more frequent than their neighbouring values, a run of equally frequent values counting together.  \code{discrete = TRUE} instead rounds to integers and returns the most frequent integers.
#' @author Fred Viole, OVVO Financial Systems
#' @examples
#' \dontrun{
#' set.seed(123)
#' x <- rnorm(100)
#' NNS.mode(x)
#' }
#' @export


NNS.mode <- function(x, discrete = FALSE, multi = TRUE) {
    .Call(`_NNS_NNS_mode_cpp`, as.numeric(x), as.logical(discrete), as.logical(multi))
}


mode <- function(x) NNS.mode(x, discrete = FALSE, multi = FALSE)

mode_class <- function(x) NNS.mode(x, discrete = TRUE, multi = FALSE)


#' NNS gravity
#'
#' Central tendency averaging the median, the mean, the midpoint of the quartiles and a histogram mode; less sensitive to outliers than the mean.
#'
#' @param x vector of data.
#' @param discrete logical; \code{FALSE} (default) for continuous distributions, \code{TRUE} for discrete distributions.
#' @return Returns a numeric value representing the central tendency of the distribution.
#' @details \deqn{(Q_2 + m + \bar{x} + (Q_1 + Q_3)/2) / 4,} with the quartiles \eqn{Q_1, Q_2, Q_3} of \code{quantile(x, type = 7)} and, as the mode \eqn{m}, the mean of the observations in the tallest bin of a histogram with bin width IQR / sqrt(n) and its two neighbouring bins (all observations when several bins are equally tall).  Up to three values give the median.  The mean term remains unbounded under one-point contamination: a single 1e6 among 99 zeros gives 2500, against 10000 for the mean. The other components can also change, so the one-quarter factor is an example, not a universal sensitivity bound. Its asymptotic replacement breakdown point is zero; it is not a high-breakdown estimator.
#' @author Fred Viole, OVVO Financial Systems
#' @examples
#' \dontrun{
#' set.seed(123)
#' x <- rnorm(100)
#' NNS.gravity(x)
#' }
#' @export

NNS.gravity <- function(x, discrete = FALSE) {
  .Call(`_NNS_NNS_gravity_cpp`, as.numeric(x), as.logical(discrete))
}

gravity <- function(x) NNS.gravity(x, discrete = FALSE)

gravity_class <- function(x) NNS.gravity(x, discrete = TRUE)


#' NNS rescale
#'
#' Rescale a vector using either min-max scaling or risk-neutral adjustment.
#'
#' @param x numeric vector; data to rescale (e.g., terminal prices for risk-neutral method).
#' @param a numeric; defines the scaling target:
#'   - For \code{method = "minmax"}: the lower limit of the output range (e.g., 5 to scale to [5, b]).
#'   - For \code{method = "riskneutral"}: the initial price \( S_0 \) (must be positive, e.g., 100), used to set the target mean.
#' @param b numeric; defines the scaling range or rate:
#'   - For \code{method = "minmax"}: the upper limit of the output range (e.g., 10 to scale to [a, 10]).
#'   - For \code{method = "riskneutral"}: the risk-free rate \( r \) (e.g., 0.05), used with \( T \) to adjust the mean.
#' @param method character; scaling method: \code{"minmax"} (default) for min-max scaling, or \code{"riskneutral"} for risk-neutral adjustment.
#' @param T numeric; time to maturity in years (required for \code{method = "riskneutral"}, ignored otherwise; e.g., 1). Default is NULL.
#' @param type character; for \code{method = "riskneutral"}: \code{"Terminal"} (default) or \code{"Discounted"} (mean = \( S_0 \)).
#' @return Returns a rescaled distribution:
#'   - For \code{"minmax"}: values scaled linearly to the range \code{[a, b]}.
#'   - For \code{"riskneutral"}: values scaled multiplicatively to a risk-neutral mean (\( S_0 e^(rT) \) if \code{type = "Terminal"}, or \( S_0 \) if \code{type = "Discounted"}).
#' @author Fred Viole, OVVO Financial Systems
#' @examples
#' \dontrun{
#' set.seed(123)
#' # Min-max scaling: a = lower limit, b = upper limit
#' x <- rnorm(100)
#' NNS.rescale(x, a = 5, b = 10, method = "minmax")  # Scales to [5, 10]
#' 
#' # Risk-neutral scaling (Terminal): a = S_0, b = r  # Mean approx 105.13
#' prices <- 100 * exp(cumsum(rnorm(100, 0.001, 0.02)))
#' NNS.rescale(prices, a = 100, b = 0.05, method = "riskneutral", T = 1, type = "Terminal")
#' 
#' # Risk-neutral scaling (Discounted): a = S_0, b = r  # Mean approx 100
#' NNS.rescale(prices, a = 100, b = 0.05, method = "riskneutral", T = 1, type = "Discounted")
#' }
#' @export

NNS.rescale <- function(x, a, b, method = "minmax", T = NULL, type = "Terminal") {
  .Call(`_NNS_NNS_rescale_cpp`, as.numeric(x), as.numeric(a), as.numeric(b),
        as.character(method), if (is.null(T)) NULL else as.numeric(T), as.character(type))
}