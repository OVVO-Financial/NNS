#' NNS Distance
#'
#' Internal kernel function for NNS multivariate regression \link{NNS.reg} parallel instances.
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
  l <- nrow(rpm)

  y.hat <- rpm$y.hat

  cols <- colnames(rpm)[names(rpm)!="y.hat"]


  if(type!="FACTOR"){
    rpm <- rbind(as.list(t(dist.estimate)), rpm[, .SD, .SDcols = cols])
    rpm[, names(rpm) := lapply(.SD, as.numeric)]
    rpm <- rpm[,lapply(.SD, function(b) (b - min(b)) / max(1e-10, (max(b) - min(b))))]
    dist.estimate <- as.numeric(rpm[1, ])
    rpm <- rpm[-1,]
  }

  rpm$y.hat <- y.hat


  if(type=="L2"){
    rpm[, "Sum" := ( Reduce(`+`, lapply(1 : n, function(i) (rpm[[i]]-as.numeric(dist.estimate)[i])^2)        ))]#[,Sum]
    rpm$Sum <- rpm$Sum + 1/(1+Reduce(`+`, Map("==", rpm[, 1:n], as.numeric(dist.estimate))))
  }

  if(type=="L1"){
    rpm[, "Sum" := ( Reduce(`+`, lapply(1 : n, function(i) abs(rpm[[i]]-as.numeric(dist.estimate)[i]))))]#[,Sum]
    rpm$Sum <- rpm$Sum + 1/(1+Reduce(`+`, Map("==", rpm[, 1:n], as.numeric(dist.estimate))))
  }

  if(type=="DTW"){
    rpm[, "Sum" := ( unlist(lapply(1 : nrow(rpm), function(i) dtw::dtw(as.numeric(rpm[i, ]), as.numeric(dist.estimate))$distance)))]#[,Sum]
  }

  if(type=="FACTOR"){
    rpm[, "Sum"  := ( 1/(1+Reduce(`+`, Map("==", rpm[, 1:n], as.numeric(dist.estimate)))))]#[,Sum]
  }

  rpm$Sum[rpm$Sum == 0] <- 1e-10

  data.table::setkey(rpm, Sum)

  if(k==1){
    index <- which.min(rpm$Sum)
    if(length(index)>1){
        return(mode(rpm$y.hat[index]))
    }  else {
        return(rpm$y.hat[1])
    }
  }

  rpm <- rpm[1:min(k,l),]

  weights <- (1 / rpm$Sum) / sum(1 / rpm$Sum)

  norm <- dnorm(rpm$Sum)
  norm_weights <- norm / sum(norm)

  if(type!="FACTOR"){
    weights <- rowMeans(cbind(weights, norm_weights, rep(1/min(k,l), length(weights))))
    weights.sum <- sum(weights)
    weights <- weights / weights.sum
  }

  single.estimate <- sum(weights * rpm$y.hat)

  return(single.estimate)
}
