#' NNS Causation
#'
#' Returns the causality from observational data between two variables.
#'
#' @param x a numeric vector, matrix or data frame.
#' @param y \code{NULL} (default) or a numeric vector with compatible dimsensions to \code{x}.
#' @param factor.2.dummy logical; \code{TRUE} (default) Automatically augments variable matrix with numerical dummy variables based on the levels of factors.  Includes dependent variable \code{y}.
#' @param tau options: ("cs", "ts", integer); Number of lagged observations to consider (for time series data).  Otherwise, set \code{(tau = "cs")} for cross-sectional data.  \code{(tau = "ts")} automatically selects the lag of the time series data, while \code{(tau = [integer])} specifies a time series lag.
#' @param time.series logical; \code{FALSE} (default) If analyzing time series data with \code{tau = [integer]}, select \code{(time.series = TRUE)}.  Not required when \code{(tau = "ts")}.
#' @param plot logical; \code{FALSE} (default) Plots the raw variables, tau normalized, and cross-normalized variables.
#' @return Returns the directional causation (x ---> y) or (y ---> x) and net quantity of association.  For causal matrix, directional causation is returned as ([column variable] ---> [row variable]).  Negative numbers represent causal direction attributed to [row variable].
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' ## x clearly causes y...
#' \dontrun{
#' set.seed(123)
#' x <- rnorm(100) ; y <- x ^ 2
#' NNS.caus(x, y, tau = "cs")
#'
#' x <- 1:100 ; y <- x^2
#' NNS.caus(x, y, tau = "ts", time.series = TRUE)}
#'
#' ## Causal matrix
#' \dontrun{
#' NNS.caus(data.matrix(iris), tau = 0)
#' }
#' @export

NNS.caus <- function(x, y,
                     factor.2.dummy = TRUE,
                     tau,
                     time.series=FALSE,
                     plot=FALSE){

  orig.tau <- tau
  orig.time.series <- time.series
  orig.plot <- plot

  if(factor.2.dummy){
      if(!is.null(dim(x))){
          if(!is.numeric(x)){
              x <- sapply(x,factor_2_dummy_FR)
          } else {
              x <- apply(x,2,as.double)
          }
          if(is.list(x)){
              x <- do.call(cbind,x)
              x <- apply(x,2,as.double)
          }

      } else {
          x <- factor_2_dummy(x)
          if(is.null(dim(x))){
              x <- as.double(x)
          } else {
              x <- apply(x,2,as.double)
          }
      }
  }

  if(!missing(y)){
      if(is.numeric(tau)){
          Causation.x.given.y <- Uni.caus(x,y,tau=tau,plot = FALSE,time.series=time.series)
          Causation.y.given.x <- Uni.caus(y,x,tau=tau,plot = FALSE,time.series=time.series)

      if(Causation.x.given.y == Causation.y.given.x |
         Causation.x.given.y == 0 | Causation.y.given.x == 0){
            Causation.x.given.y <- Uni.caus(x, y, tau = tau, plot = FALSE, scale = TRUE, time.series = time.series)
            Causation.y.given.x <- Uni.caus(y, x, tau = tau, plot = FALSE, scale = TRUE, time.series = time.series)
      }
    }

    if(tau == "cs"){
        Causation.x.given.y <- Uni.caus(x, y, tau = 0, plot = FALSE)
        Causation.y.given.x <- Uni.caus(y, x, tau = 0, plot = FALSE)

        if(Causation.x.given.y == Causation.y.given.x |
            Causation.x.given.y == 0 | Causation.y.given.x == 0){
                Causation.x.given.y <- Uni.caus(x, y, tau = 0, plot = FALSE, scale = TRUE)
                Causation.y.given.x <- Uni.caus(y, x, tau = 0, plot = FALSE, scale = TRUE)
        }
    }

    if(tau == "ts"){
        Causes.xy <- numeric()
        Causes.yx <- numeric()

        for(i in 0 : 4){
        # Populate scaling taus and calculate uni causation
            Causes.xy[i] <- Uni.caus(x, y, tau = ceiling((2 ^ i) / 100 * length(x)), plot = FALSE, time.series = TRUE) / (2 ^ i)
            Causes.yx[i] <- Uni.caus(y, x, tau = ceiling((2 ^ i) / 100 * length(x)) , plot = FALSE, time.series = TRUE) / (2 ^ i)
        }

        Causation.x.given.y <- mean(Causes.xy)

        Causation.y.given.x <- mean(Causes.yx)
    }


    if(abs(Causation.x.given.y) <= abs(Causation.y.given.x)){
        if(plot){
            # For plotting only
            if(tau == "cs"){
                tau <- 0
            }
            if(tau == "ts"){
                tau <- ceiling(0.03 * length(x))
            }
            Uni.caus(y, x, tau = tau, plot = plot)
        }
        return(c(Causation.x.given.y = Causation.x.given.y,
               Causation.y.given.x = Causation.y.given.x,
               "C(y--->x)" =  Causation.y.given.x - Causation.x.given.y))
    } else {
        if(plot){
            # For plotting only
            if(tau == "cs"){
                tau <- 0
            }
            if(tau == "ts"){
                tau <- ceiling(0.03 * length(x))
            }
            Uni.caus(x, y, tau = tau, plot = plot)
        }
    return(c(Causation.x.given.y = Causation.x.given.y,
               Causation.y.given.x = Causation.y.given.x,
               "C(x--->y)" = Causation.x.given.y - Causation.y.given.x))
    }
  } else {

    NNS.caus.matrix(x, tau = orig.tau, time.series = orig.time.series, plot = orig.plot)
  }


}
