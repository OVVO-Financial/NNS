#' NNS Distance
#'
#' Internal function for NNS multivariate regression \link{NNS.reg} parallel instances.
#' @param rpm REGRESSION.POINT.MATRIX from \link{NNS.reg}
#' @param dist.estimate Vector to generate distances from.
#' @param type "L1", "L2" or "DTW"
#' @param k \code{n.best} from \link{NNS.reg}
#'
#' @return Returns sum of weighted distances.
#'
#'
#' @export

NNS.distance <- function(rpm, dist.estimate, type, k){
  type <- toupper(type)
  n <- length(dist.estimate)

  y.hat <- rpm$y.hat

  cols <- names(rpm)[names(rpm)!="y.hat"]

  if(type!="FACTOR"){
    rpm <- rbind(as.list(t(dist.estimate)), rpm[, .SD, .SDcols = cols])
    rpm[, names(rpm) := lapply(.SD, as.numeric)]
    rpm <- rpm[,lapply(.SD, function(b) (b - min(b)) / max(1e-10, (max(b) - min(b))))]
    dist.estimate <- as.numeric(rpm[1, ])
    rpm <- rpm[-1,]
  }

  rpm$y.hat <- y.hat


  if(type=="L2"){
    row.sums <- rpm[,  `:=` (Sum = Reduce(`+`, lapply(1 : n, function(i) (rpm[[i]]-as.numeric(dist.estimate)[i])^2)))][,Sum]
  }

  if(type=="L1"){
    row.sums <- rpm[,  `:=` (Sum = Reduce(`+`, lapply(1 : n, function(i) abs(rpm[[i]]-as.numeric(dist.estimate)[i]))))][,Sum]
  }

  if(type=="DTW"){
    row.sums <- rpm[,  `:=` (Sum = unlist(lapply(1 : nrow(rpm), function(i) dtw(as.numeric(rpm[i, ]), as.numeric(dist.estimate))$distance)))][,Sum]
  }

  if(type=="FACTOR"){
    row.sums <- rpm[,  `:=` (Sum = 1/Reduce(`+`, Map("==", rpm[, 1:n], as.numeric(dist.estimate))))][,Sum]
  }

  row.sums[row.sums == 0] <- 1e-10

  if(k==1){
    if(length(which(row.sums == min(row.sums)))>1){
      return(mode(rpm$y.hat[which(row.sums == min(row.sums))][1]))
    }  else {
      return(rpm$y.hat[which.min(row.sums)][1])
    }
  }

  total.row.sums <- sum(1 / row.sums)
  weights <- (1 / row.sums) / total.row.sums

  highest <- rev(order(weights))[1 : min(k, length(weights))]

  weights[-highest] <- 0

  weights.sum <- sum(weights)

  weights <- weights / weights.sum

  weights <- rowMeans(cbind(weights, rep(1/k, length(weights))))

  weights[-highest] <- 0

  weights.sum <- sum(weights)

  weights <- weights / weights.sum

  single.estimate <- sum(weights * rpm$y.hat)

  rpm[,"Sum":=NULL]

  return(mean(single.estimate))
}
