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
#' set.seed(123)
#' x <- rnorm(100)
#' NNS.mode(x)
#' @export


NNS.mode <- function(x, discrete = FALSE, multi = TRUE){
  x <- as.numeric(x)
  l <- length(x)
  if(l <= 3) return(median(x))
  if(length(unique(x))==1) return(x[1])
  x_s <- x[order(x)]
  range <- abs(x_s[l]-x_s[1])
  if(range==0) return(x[1])
  
  z <- MESS::bin(x_s, range/128, origin = x_s[1], missinglast = FALSE)
  lz <- length(z$counts)
  max_z <- z$counts==max(z$counts)
  z_names <- seq(x_s[1], x_s[l], z$width)
  
  if(sum(max_z)>1){
    z_ind <- 1:lz
    if(multi) return(z_names[max_z])
  } else {
    z_c <- which.max(z$counts)
    z_ind <- max(1, (z_c - 1)):min(lz,(z_c + 1))
  }
  
  
  
  final <- sum(z_names[z_ind] * z$counts[z_ind] )/sum(z$counts[z_ind])
  if(discrete) return(ifelse(final%%1 < .5, floor(final), ceiling(final))) else return(final)
}