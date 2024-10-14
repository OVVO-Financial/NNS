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

  z <- NNS_bin(x_s, range/128, origin = x_s[1], missinglast = FALSE)
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
#' Rescale min-max scaling output between two numbers.
#'
#' @param x vector of data.
#' @param a numeric; lower limit.
#' @param b numeric; upper limit.
#' @return Returns a rescaled distribution within provided limits.
#' @author Fred Viole, OVVO Financial Systems
#' @examples
#' \dontrun{
#' set.seed(123)
#' x <- rnorm(100)
#' NNS.rescale(x, 5, 10)
#' }
#' @export


NNS.rescale <- function (x, a, b) {
  x <- as.numeric(x)
  output <- a + (b - a) * (x - min(x))/(max(x) - min(x))
  return(output)
}


