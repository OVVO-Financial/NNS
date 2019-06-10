#' NNS Distance
#'
#' Internal function for NNS multivariate regression \link{NNS.reg} parallel instances.
#' @param rpm REGRESSION.POINT.MATRIX from \link{NNS.reg}
#' @param dist.estimate Vector to generate distances from.
#' @param type "L1" or "L2"
#' @param k \code{n.best} from \link{NNS.reg}
#'
#' @return Returns sum of weighted distances.
#'
#'
#' @export

NNS.distance <- function(rpm,dist.estimate,type,k){
  n <- length(dist.estimate)

  if(type=="L2"){
      row.sums <- rpm[,  `:=`(Sum= Reduce(`+`, lapply(1 : n,function(i) (rpm[[i]]-as.numeric(dist.estimate)[i])^2)))][,Sum]
  } else {
      row.sums <- rpm[,  `:=`(Sum= Reduce(`+`, lapply(1 : n,function(i) abs(rpm[[i]]-as.numeric(dist.estimate)[i]))))][,Sum]
  }

  row.sums[row.sums == 0] <- 1e-10
  total.row.sums <- sum(1 / row.sums)
  weights <- (1 / row.sums) / total.row.sums

  highest <- rev(order(weights))[1 : min(k, length(weights))]

  weights[-highest] <- 0
  weights.sum <- sum(weights)

  weights <- weights / weights.sum
  single.estimate <- sum(weights * rpm$y.hat)

  return(single.estimate)
}
