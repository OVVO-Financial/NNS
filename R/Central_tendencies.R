#' NNS mode
#'
#' Mode of a distribution, either continuous or discrete.
#'
#' @param x vector of data.
#' @param discrete logical; \code{FALSE} (default) for discrete distributions.
#' @param multi logical; \code{TRUE} (default) returns multiple mode values.
#' @return Returns a numeric value representing the mode of the distribution.
#' @author Fred Viole, OVVO Financial Systems
#' @examples
#' \dontrun{
#' set.seed(123)
#' x <- rnorm(100)
#' NNS.mode(x)
#' }
#' @export


NNS.mode <- function (x, discrete = FALSE, multi = TRUE)
{
  x <- as.numeric(x)
  l <- length(x)
  if (l <= 3)
    return(median(x))
  if (length(unique(x)) == 1)
    return(x[1])
  x_s <- x[order(x)]
  range <- abs(x_s[l] - x_s[1])
  if (range == 0)
    return(x[1])
  z <- NNS_bin(x_s, range/128, origin = x_s[1], missinglast = FALSE)
  lz <- length(z$counts)
  max_z <- z$counts == max(z$counts)
  z_names <- seq(x_s[1], x_s[l], z$width)
  if (sum(max_z) > 1) {
    z_ind <- 1:lz
    if (multi)
      return(z_names[max_z])
  }
  else {
    z_c <- which.max(z$counts)
    z_ind <- max(1, (z_c - 1)):min(lz, (z_c + 1))
  }
  final <- sum(z_names[z_ind] * z$counts[z_ind])/sum(z$counts[z_ind])
  if (discrete) {
    final <- ifelse(final%%1 < 0.5, floor(final), ceiling(final))
    return(final)
  }
  else {
    if (multi) {
      return(final)
    }
    else {
      return(mean(final))
    }
  }
}



#' NNS gravity
#'
#' Alternative central tendency measure more robust to outliers.
#'
#' @param x vector of data.
#' @param discrete logical; \code{FALSE} (default) for discrete distributions.
#' @return Returns a numeric value representing the central tendency of the distribution.
#' @author Fred Viole, OVVO Financial Systems
#' @examples
#' \dontrun{
#' set.seed(123)
#' x <- rnorm(100)
#' NNS.gravity(x)
#' }
#' @export

NNS.gravity <- function (x, discrete = FALSE)
{
  l <- length(x)
  if (l <= 3) return(median(x))
  if (length(unique(x)) == 1) return(x[1])

  x_s <- x[order(x)]
  range <- abs(x_s[l] - x_s[1])

  if (range == 0)  return(x[1])

  l_25 = l*.25
  l_50 = l*.5
  l_75 = l*.75

  if(l%%2==0){
    q1 <- x_s[l_25]
    q2 <- x_s[l_50]
    q3 <- x_s[l_75]
  } else {
    f_l_25 = floor(l_25)
    f_l_75 = floor(l_75)

    q1 <- sum(x_s[f_l_25]+(l_25%%1 * (x_s[ceiling(l_25)] - x_s[f_l_25])))
    q2 <- (x_s[floor(l_50)]+x_s[ceiling(l_50)])/2
    q3 <- sum(x_s[f_l_75]+((l_75)%%1 * (x_s[ceiling(l_75)] - x_s[f_l_75])))
  }

  width = (q3 - q1) * l^(-1/2)
  
  z <- tryCatch(
    {
      NNS_bin(x_s, width, origin = x_s[1], missinglast = FALSE)
    },
    error = function(e) {
        return(NNS_bin(x_s, range / 128, origin = x_s[1], missinglast = FALSE))
      } 
    
  )
  
  lz <- length(z$counts)
  max_z <- z$counts == max(z$counts)
  if (sum(max_z) > 1) {
    z_ind <- 1:lz
  }
  else {
    z_c <- which.max(z$counts)
    z_ind <- max(1, (z_c - 1)):min(lz, (z_c + 1))
  }
  z_names <- seq(x_s[1], x_s[l], z$width)
  m <- sum(z_names[z_ind] * z$counts[z_ind])/sum(z$counts[z_ind])
  mu <- sum(x)/l
  res <- (q2 + m + mu + mean(c(q1, q2, q3)))/4
  if (is.na(res))
    final <- q2
  else final <- res
  if (discrete)
    return(ifelse(final%%1 < 0.5, floor(final), ceiling(final)))
  else return(final)
}



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
  x <- as.numeric(x)
  method <- tolower(method)
  type <- tolower(type)
  
  if (method == "minmax") {
    # Original min-max scaling
    if (max(x) == min(x)) stop("Cannot rescale: max(x) equals min(x)")
    output <- a + (b - a) * (x - min(x)) / (max(x) - min(x))
  } else if (method == "riskneutral") {
    # Risk-neutral scaling
    if (is.null(T)) stop("T (time to maturity) must be provided for riskneutral method")
    if (a <= 0) stop("S_0 (a) must be positive for riskneutral method")
    S_0 <- a  # Initial price
    r <- b    # Risk-free rate
    
    if (type == "terminal") {
      # Scale to S_0 * e^(r * T)
      theta <- log(S_0 * exp(r * T) / mean(x))
      output <- x * exp(theta)  # Mean = S_0 * e^(r * T)
    } else if (type == "discounted") {
      # Scale to S_0
      theta <- log(S_0 / mean(x))
      output <- x * exp(theta)  # Mean = S_0
    } else {
      stop("Invalid type: use 'Terminal' or 'Discounted' for riskneutral method")
    }
  } else {
    stop("Invalid method: use 'minmax' or 'riskneutral'")
  }
  
  return(output)
}