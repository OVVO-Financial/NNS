#' NNS Distance
#'
#' Internal kernel function for NNS multivariate regression \link{NNS.reg} parallel instances.
#' @param rpm REGRESSION.POINT.MATRIX from \link{NNS.reg}
#' @param dist.estimate Vector to generate distances from.
#' @param k \code{n.best} from \link{NNS.reg}
#' @param class if classification problem.
#'
#' @return Returns sum of weighted distances.
#'
#'
#' @export

NNS.distance <- function(rpm, dist.estimate, k, class){
  l <- nrow(rpm)
  if(k=="all") k <- l
  y.hat <- rpm$y.hat
  raw.dist.estimate <- unlist(dist.estimate)
  raw.rpm <- rpm[ , -"y.hat"]
  n <- length(raw.dist.estimate)
  parallel <- FALSE

  
  rpm <- rbind(as.list(t(dist.estimate)), rpm[, .SD, .SDcols = 1:n])
  rpm <- rpm[, names(rpm) := lapply(.SD, function(b) NNS.rescale(b, 0, 1)), .SDcols = 1:n]
  dist.estimate <- unlist(rpm[1, ])
  rpm <- rpm[-1,]
  
  M <- matrix(rep(dist.estimate, l), byrow = T, ncol = n)

  rpm$Sum <- Rfast::rowsums( ((t(t(rpm)) - M)^2) + abs(t(t(rpm)) - M), parallel = parallel)

  rpm$Sum[rpm$Sum == 0] <- 1e-10
  rpm$y.hat <- y.hat
  
  data.table::setkey(rpm, Sum)
  
  ll <- min(k, l)

  rpm <- rpm[1:ll,]
  
  SUM = rpm$Sum

  if(k==1){
    index <- which(SUM==min(SUM))
    if(length(index)>1){
      return(mode(rpm$y.hat[index]))
    }  else {
      return(rpm$y.hat[1])
    }
  }
  
  

  uni_weights <- rep(1/ll, ll)
  
  
  t_weights <- dt(SUM, df = ll)
  t_weights <- t_weights/sum(t_weights)
  if(any(is.na(t_weights))) t_weights <- rep(0, ll)

  emp <- SUM^(-1)
  emp_weights <- emp / sum(emp)
  if(any(is.na(emp_weights))) emp_weights <- rep(0, ll)

  exp <- dexp(1:ll, rate = 1/ll)
  exp_weights <- exp / sum(exp)
  if(any(is.na(exp_weights))) exp_weights <- rep(0, ll)

  lnorm <- abs(rev(dlnorm(1:ll, meanlog = 0, sdlog = sd(1:ll), log = TRUE)))
  lnorm_weights <- lnorm / sum(lnorm)
  if(any(is.na(lnorm_weights))) lnorm_weights <- rep(0, ll)

  pl_weights <- (1:ll) ^ (-2)
  pl_weights <- pl_weights / sum(pl_weights)
  if(any(is.na(pl_weights))) pl_weights <- rep(0, ll)

  norm_weights <- dnorm(SUM, mean = 0, sd = sd(SUM))
  norm_weights <- norm_weights / sum(norm_weights)
  if(any(is.na(norm_weights))) norm_weights <- rep(0, ll)
  
  rbf_weights <- exp(- SUM / (2*var(SUM)))
  rbf_weights <- rbf_weights / sum(rbf_weights)
  if(any(is.na(rbf_weights))) rbf_weights <- rep(0, ll)

  weights <- (emp_weights + exp_weights + lnorm_weights + norm_weights + pl_weights + t_weights + uni_weights + rbf_weights)/
    sum(emp_weights + exp_weights + lnorm_weights + norm_weights + pl_weights + t_weights + uni_weights + rbf_weights)
  
  
  if(is.null(class)) single.estimate <- rpm$y.hat%*%weights else single.estimate <- mode_class(rep(rpm$y.hat, ceiling(100*weights)))
  

  return(single.estimate)
}
