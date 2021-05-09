#' NNS Distance
#'
#' Internal kernel function for NNS multivariate regression \link{NNS.reg} parallel instances.
#' @param rpm REGRESSION.POINT.MATRIX from \link{NNS.reg}
#' @param rpm_class integer \code{rpm}.
#' @param dist.estimate Vector to generate distances from.
#' @param type "L1", "L2", "DTW" or "FACTOR"
#' @param k \code{n.best} from \link{NNS.reg}
#' @param n number of observations.
#'
#' @return Returns sum of weighted distances.
#'
#'
#' @export

NNS.distance <- function(rpm, rpm_class, dist.estimate, type, k, n){
  type <- tolower(type)
  n <- length(dist.estimate)
  l <- nrow(rpm)
  y.hat <- rpm$y.hat
  raw.dist.estimate <- unlist(dist.estimate)

  if(type!="factor"){
    rpm <- rbind(as.list(t(dist.estimate)), rpm[, .SD, .SDcols = 1:n])
    rpm <- rpm[, names(rpm) := lapply(.SD, function(b) (b - min(b)) / max(1e-10, (max(b) - min(b)))), .SDcols = 1:n]
    dist.estimate <- unlist(rpm[1, ])
    rpm <- rpm[-1,]
  }

  rpm$y.hat <- y.hat

  if(type=="l2"){
    rpm$Sum <- Rfast::rowsums(t(t(rpm[,1:n]) - dist.estimate)^2, parallel = TRUE)
  }

  if(type=="l1"){
    rpm$Sum <- Rfast::rowsums(abs(t(t(rpm[,1:n]) - dist.estimate)), parallel = TRUE)
  }

  if(type=="dtw"){
    rpm[, "Sum" := ( unlist(lapply(1 : nrow(rpm), function(i) dtw::dtw(as.numeric(rpm[i, 1:n]), as.numeric(dist.estimate))$distance)))]
  }


  rpm$Sum2 <-  Rfast::rowsums( (1/l + t(t(rpm_class == raw.dist.estimate)))^-1, parallel = TRUE)

  rpm$Sum[rpm$Sum == 0] <- 1e-10

  data.table::setorder(rpm, Sum, Sum2)
  rpm <- rpm[1:min(k,l),]

  if(k==1){
    index <- which(rpm$Sum==min(rpm$Sum))
    if(length(index)>1){
      return(mode(rpm$y.hat[index]))
    }  else {
      return(rpm$y.hat[1])
    }
  }

  unif_weights <- 1/min(k,l)

  inv <- 1 / rpm$Sum
  emp_weights <- inv / sum(inv)

  norm_weights <- pmin(max(rpm$Sum), dnorm(rpm$Sum, mean(rpm$Sum), sd(rpm$Sum)))
  norm_weights <- norm_weights / sum(norm_weights)

  t_weights <- dt(x = rpm$Sum, df = min(k,l)-1)
  t_weights <- t_weights/sum(t_weights)

  weights <- emp_weights + norm_weights + t_weights + unif_weights
  weights <- weights / sum(weights)

  single.estimate <- weights%*%rpm$y.hat

  inv <- 1 / rpm$Sum2
  emp_weights <- inv / sum(inv)

  norm_weights <- pmin(max(rpm$Sum2), dnorm(rpm$Sum2, mean(rpm$Sum2), sd(rpm$Sum2)))
  norm_weights <- norm_weights / sum(norm_weights)

  t_weights <- dt(x = rpm$Sum2, df = min(k,l)-1)
  t_weights <- t_weights/sum(t_weights)

  weights <- emp_weights + norm_weights + t_weights + unif_weights
  weights <- weights / sum(weights)

  single.estimate2 <- weights%*%rpm$y.hat

  return((single.estimate + single.estimate2)/2)
}
