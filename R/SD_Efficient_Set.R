#' NNS SD Efficient Set
#'
#' Determines the set of stochastic dominant variables for various degrees.
#'
#' @param x a numeric matrix or data frame.
#' @param degree numeric options: (1, 2, 3); Degree of stochastic dominance test from (1, 2 or 3).
#' @param type options: ("discrete", "continuous"); \code{"discrete"} (default) selects the type of CDF.
#' @param status logical; \code{TRUE} (default) Prints status update message in console.
#' @return Returns set of stochastic dominant variable names.
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2016) "LPM Density Functions for the Computation of the SD Efficient Set." Journal of Mathematical Finance, 6, 105-126.  DOI: \doi{10.4236/jmf.2016.61012}.
#'
#' Viole, F. (2017) "A Note on Stochastic Dominance." \url{https://www.ssrn.com/abstract=3002675}.
#' @examples
#' set.seed(123)
#' x <- rnorm(100) ; y<-rnorm(100) ; z<-rnorm(100)
#' A <- cbind(x, y, z)
#' NNS.SD.efficient.set(A, 1)
#' @export



NNS.SD.efficient.set <- function(x, degree, type = "discrete", status = TRUE) {
  type <- tolower(type)
  if(!any(type %in% c("discrete", "continuous")))
    warning("type needs to be either 'discrete' or 'continuous'")
  if(!any(degree %in% c(1,2,3)))
    warning("degree needs to be 1, 2, or 3")
  if(any(class(x) %in% c("tbl","data.table"))) 
    x <- as.data.frame(x)

  n <- ncol(x)
  max_target <- max(x)
  LPM_order <- numeric()
  Dominated_set <- numeric()
  current_base <- numeric()

  if(is.null(colnames(x))) 
    colnames(x) <- paste0("X_",1:ncol(x))

  LPM_order <- sapply(1 : n, function(i) LPM(1, max_target, x[ , i]))
  final_ranked <- x[ , order(LPM_order)]
  all_variables <- colnames(final_ranked)
  current_base <- 1
  for (i in 1:(n-1)) {
    if(status) 
      message("Checking ", i, " of ", (n-1), "\r", appendLF=FALSE)
    if(i == (n-1) & status) 
      message("                                        ", appendLF=TRUE)

    base <- final_ranked[ , tail(current_base, 1)]
    challenger <- final_ranked[ , i + 1]
    if(degree == 1){
      sd.test <- NNS.FSD.uni(base, challenger, type = type)
    }else if(degree == 2){
      sd.test <- NNS.SSD.uni(base, challenger)
    }else if(degree == 3){
      sd.test <- NNS.TSD.uni(base, challenger)
    }
    if(sd.test == 1){
      Dominated_set[i] <- i + 1
    } else {
      I <- FALSE
      for (j in current_base){
        base <- final_ranked[ , j]
        if(degree == 1){
          new.base.sd.test <- NNS.FSD.uni(base, challenger, type = type)
        }else if(degree == 2){
          new.base.sd.test <- NNS.SSD.uni(base, challenger)
        }else if(degree == 3){
          new.base.sd.test <- NNS.TSD.uni(base, challenger)
        }
        if (new.base.sd.test == 0){
          I <- FALSE
          next
        } else {
          I <- TRUE
          Dominated_set[i] <- i + 1
          break
        }
      }
      if(!I){
        current_base <- c(current_base, i + 1)
      }
    }
  }
  if(length(Dominated_set) > 0){
    return(all_variables[-na.omit(Dominated_set)])
  }
  return(all_variables)
}
