#' NNS SD Efficient Set
#'
#' Determines the set of stochastic dominant variables for various degrees.
#' @param x a numeric matrix or data frame.
#' @param degree numeric options: (1, 2, 3); Degree of stochastic dominance test from (1, 2 or 3).
#' @param type options: ("discrete", "continuous"); \code{"discrete"} (default) selects the type of CDF.
#' @return Returns set of stochastic dominant variable names.
#' @keywords stochastic dominance
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2016) "LPM Density Functions for the Computation of the SD Efficient Set." Journal of Mathematical Finance, 6, 105-126. \url{http://www.scirp.org/Journal/PaperInformation.aspx?PaperID=63817}.
#'
#' Viole, F. (2017) "A Note on Stochastic Dominance." \url{https://ssrn.com/abstract=3002675}.
#' @examples
#' set.seed(123)
#' x <- rnorm(100) ; y<-rnorm(100) ; z<-rnorm(100)
#' A <- cbind(x, y, z)
#' NNS.SD.efficient.set(A, 1)
#' @export



NNS.SD.efficient.set <- function(x,degree,type="discrete") {
  n <- ncol(x)
  max_target <- max(x)
  LPM_order <- numeric()
  Dominated_set <- numeric()
  current_base <- numeric()


  LPM_order = sapply(1 : n,function(i) LPM(1, max_target, x[ , i]))

  final_ranked <- x[ , order(LPM_order)]

  current_base <- 1


  for (i in 1:(n-1)) {

      base <- final_ranked[ , current_base[length(current_base)]]

      challenger <- final_ranked[ , i + 1]

    if(degree == 1){
      sd.test = NNS.FSD.uni(base, challenger, type = type)
    }
    if(degree == 2){
      sd.test = NNS.SSD.uni(base, challenger)
    }
    if(degree == 3){
      sd.test = NNS.TSD.uni(base, challenger)
    }

      if (sd.test == 1){
        current_base[i] <- current_base[length(current_base)]
        Dominated_set[i] <- i + 1
      }


      if (sd.test == 0){
        for (j in current_base){
          base <- final_ranked[ , j]
          if(degree == 1){
            new.base.sd.test = NNS.FSD.uni(base, challenger, type = type)
          }
          if(degree == 2){
            new.base.sd.test = NNS.SSD.uni(base, challenger)
          }
          if(degree == 3){
            new.base.sd.test = NNS.TSD.uni(base, challenger)
          }

          if (new.base.sd.test == 0){ next
          } else {
            Dominated_set[i] <- i + 1
            }
        }

        current_base[i]<- i + 1
      }

  }


  if(length(Dominated_set) > 0){
      return(colnames(final_ranked[ , - na.omit(Dominated_set)]))
  } else {
      return(colnames(final_ranked))
    }


}
