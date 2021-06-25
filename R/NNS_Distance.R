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
  type <- toupper(type)
  l <- nrow(rpm)
  y.hat <- rpm$y.hat
  raw.dist.estimate <- unlist(dist.estimate)
  n <- length(raw.dist.estimate)



  if(type!="FACTOR"){
    rpm <- rbind(as.list(t(dist.estimate)), rpm[, .SD, .SDcols = 1:n])
    rpm <- rpm[, names(rpm) := lapply(.SD, function(b) (b - min(b)) / max(1e-10, (max(b) - min(b)))), .SDcols = 1:n]
    dist.estimate <- unlist(rpm[1, ])
    rpm <- rpm[-1,]
  }

  rpm$y.hat <- y.hat

  if(type=="L2"){
    rpm$Sum <- Rfast::rowsums( t((t(rpm[, 1:n]) - dist.estimate)^2) * ((l - (rpm_class == raw.dist.estimate))/l), parallel = TRUE)
  }

  if(type=="L1"){
    rpm$Sum <- Rfast::rowsums(abs(t(t(rpm[, 1:n]) - dist.estimate)) * ((l - (rpm_class == raw.dist.estimate))/l), parallel = TRUE)
  }

  if(type=="DTW"){
    rpm[, "Sum" := ( unlist(lapply(1 : nrow(rpm), function(i) dtw::dtw(as.numeric(rpm[i, 1:n]), as.numeric(dist.estimate))$distance)))]
  }

  if(type=="FACTOR"){
    rpm$Sum <- (1/l + ( Rfast::rowsums((rpm_class == raw.dist.estimate), parallel = TRUE)))^-1
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

  uni_weights <- rep(1/min(k,l), min(k,l))

  emp <- rpm$Sum^(-1/min(k,l))
  emp_weights <- emp / sum(emp)

  exp <- dexp(1:min(k,l), rate = 1/min(k,l))
  exp_weights <- exp / sum(exp)

  lnorm <- abs(rev(dlnorm(1:min(k, l), meanlog = min(k, l), sdlog = min(k, l), log = TRUE)))
  lnorm_weights <- lnorm / sum(lnorm)

  weights <- (emp_weights + exp_weights + lnorm_weights + uni_weights)/sum(emp_weights + exp_weights + lnorm_weights + uni_weights)

  single.estimate <- rpm$y.hat%*%weights

  return(single.estimate)
}
