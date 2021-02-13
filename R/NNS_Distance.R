#' NNS Distance
#'
#' Internal kernel function for NNS multivariate regression \link{NNS.reg} parallel instances.
#' @param rpm REGRESSION.POINT.MATRIX from \link{NNS.reg}
#' @param dist.estimate Vector to generate distances from.
#' @param type "L1", "L2", "DTW" or "FACTOR"
#' @param k \code{n.best} from \link{NNS.reg}
#' @param n number of observations.
#'
#' @return Returns sum of weighted distances.
#'
#'
#' @export

NNS.distance <- function(rpm, dist.estimate, type, k, n){
  type <- toupper(type)
  n <- length(dist.estimate)
  l <- nrow(rpm)
  y.hat <- rpm$y.hat

  if(type!="FACTOR"){
    rpm <- rbind(as.list(t(dist.estimate)), rpm[, .SD, .SDcols = 1:n])
    rpm <- rpm[, names(rpm) := lapply(.SD, function(b) (b - min(b)) / max(1e-10, (max(b) - min(b)))), .SDcols = 1:n]
    dist.estimate <- unlist(rpm[1, ])
    rpm <- rpm[-1,]
  }

  rpm$y.hat <- y.hat

  rpm_mat <- t(rpm[, 1:n])


  if(type=="L2"){
    rpm$Sum <- Rfast::rowsums(t(rpm_mat - dist.estimate)^2 + 1/(1/l + t(ifelse(rpm_mat%%1 < .5, floor(rpm_mat), ceiling(rpm_mat)) == dist.estimate)), parallel = TRUE)
  }

  if(type=="L1"){
    rpm$Sum <- Rfast::rowsums(abs(t(rpm_mat - dist.estimate)) + 1/(1/l + t(ifelse(rpm_mat%%1 < .5, floor(rpm_mat), ceiling(rpm_mat))  == dist.estimate)), parallel = TRUE)
  }

  if(type=="DTW"){
    rpm[, "Sum" := ( unlist(lapply(1 : nrow(rpm), function(i) dtw::dtw(as.numeric(rpm[i, 1:n]), as.numeric(dist.estimate))$distance)))]
  }

  if(type=="FACTOR"){
    rpm$Sum <- 1/(1/l + ( Rfast::rowsums(t(ifelse(rpm_mat%%1 < .5, floor(rpm_mat), ceiling(rpm_mat)) == dist.estimate), parallel = TRUE)))
  }

  rpm$Sum[rpm$Sum == 0] <- 1e-10

  data.table::setkey(rpm, Sum)

  if(k==1){
    index <- which(rpm$Sum==min(rpm$Sum))
    if(length(index)>1){
      return(mode(rpm$y.hat[index]))
    }  else {
      return(rpm$y.hat[1])
    }
  }

  rpm <- rpm[1:min(k,l),]

  inv <- 1 / rpm$Sum
  emp_weights <- inv / sum(inv)

  norm_weights <- pmin(max(rpm$Sum), dnorm(rpm$Sum, mean(rpm$Sum), sd(rpm$Sum)))
  norm_weights <- norm_weights / sum(norm_weights)

  t_weights <- pmin(max(rpm$Sum), dt(x = rpm$y.hat, df = min(k,l)))
  t_weights <- t_weights/sum(t_weights)

  weights <- (emp_weights + norm_weights + t_weights)/3

  single.estimate <- weights%*%rpm$y.hat

  return(single.estimate)
}
