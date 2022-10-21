#' NNS Distance
#'
#' Internal kernel function for NNS multivariate regression \link{NNS.reg} parallel instances.
#' @param rpm REGRESSION.POINT.MATRIX from \link{NNS.reg}
#' @param rpm_class integer \code{rpm}.
#' @param dist.estimate Vector to generate distances from.
#' @param type "L1", "L2" or "FACTOR"
#' @param k \code{n.best} from \link{NNS.reg}
#' @param class if classification problem.
#'
#' @return Returns sum of weighted distances.
#'
#'
#' @export

NNS.distance <- function(rpm, rpm_class, dist.estimate, type, k, class){
  type <- toupper(type)
  l <- nrow(rpm)
  y.hat <- rpm$y.hat
  raw.dist.estimate <- unlist(dist.estimate)
  n <- length(raw.dist.estimate)
  parallel <- FALSE

  if(type!="FACTOR"){
    rpm <- rbind(as.list(t(dist.estimate)), rpm[, .SD, .SDcols = 1:n])
    rpm <- rpm[, names(rpm) := lapply(.SD, function(b) ((((b - min(b))^2) / max(1e-10, (max(b) - min(b))^2)) + (((b - min(b))) / max(1e-10, (max(b) - min(b)))))/2),
               .SDcols = 1:n]
    dist.estimate <- unlist(rpm[1, ])
    rpm <- rpm[-1,]
  }

  rpm$y.hat <- y.hat

  
  if(type=="L2"){
    rpm$Sum <- Rfast::rowsums( t((t(rpm[, 1:n]) - dist.estimate)^2) * ((l - (rpm_class == raw.dist.estimate))/l), parallel = parallel)
  }

  if(type=="L1"){
    rpm$Sum <- Rfast::rowsums(abs(t(t(rpm[, 1:n]) - dist.estimate)) * ((l - (rpm_class == raw.dist.estimate))/l), parallel = parallel)
  }

  if(type=="FACTOR"){
    rpm$Sum <- (1/l + ( Rfast::rowsums((rpm_class == raw.dist.estimate), parallel = parallel)))^-1
  }

  rpm$Sum[rpm$Sum == 0] <- 1e-10

  data.table::setkey(rpm, Sum)

  rpm <- rpm[1:min(k,l),]

  if(k==1){
    index <- which(rpm$Sum==min(rpm$Sum))
    if(length(index)>1){
      return(mode(rpm$y.hat[index]))
    }  else {
      return(rpm$y.hat[1])
    }
  }

  uni_weights <- rep(1/min(k,l), min(k,l))

  emp <- rpm$Sum^(-1)
  emp_weights <- emp / sum(emp)
  if(any(is.na(emp_weights))) emp_weights <- 0

  exp <- dexp(1:min(k,l), rate = 1/min(k,l))
  exp_weights <- exp / sum(exp)
  if(any(is.na(exp_weights))) exp_weights <- 0

  lnorm <- abs(rev(dlnorm(1:min(k, l), meanlog = 0, sdlog = sd(1:min(k, l)), log = TRUE)))
  lnorm_weights <- lnorm / sum(lnorm)
  if(any(is.na(lnorm_weights))) lnorm_weights <- 0

  pl_weights <- (1:min(k, l)) ^ (-2)
  pl_weights <- pl_weights / sum(pl_weights)
  if(any(is.na(pl_weights))) pl_weights <- 0

  norm_weights <- dnorm(rpm$Sum, mean = 0, sd = sd(rpm$Sum))
  norm_weights <- norm_weights / sum(norm_weights)
  if(any(is.na(norm_weights))) norm_weights <- 0

  weights <- (emp_weights + exp_weights + lnorm_weights + norm_weights + pl_weights + uni_weights)/
    sum(emp_weights + exp_weights +  lnorm_weights + norm_weights + pl_weights + uni_weights)



  if(is.null(class)) single.estimate <- rpm$y.hat%*%weights else{
    single.estimate <- mode_class(rep(rpm$y.hat, ceiling(100*weights)))
  }

  return(single.estimate)
}
