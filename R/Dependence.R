#' NNS Dependence
#'
#' Returns the dependence and nonlinear correlation between two variables based on higher order partial moment matrices measured by frequency or area.
#'
#' @param x a numeric vector, matrix or data frame.
#' @param y \code{NULL} (default) or a numeric vector with compatible dimensions to \code{x}.
#' @param asym logical; \code{FALSE} (default) Allows for asymmetrical dependencies.
#' @param p.value logical; \code{FALSE} (default) Generates 100 independent random permutations to test results against and plots 95 percent confidence intervals along with all results.
#' @param print.map logical; \code{FALSE} (default) Plots quadrant means, or p-value replicates.
#' @return Returns the bi-variate \code{"Correlation"} and \code{"Dependence"} or correlation / dependence matrix for matrix input.
#'
#' @note
#' \code{NNS.cor} has been deprecated \code{(NNS >= 0.5.4)} and can be called via \code{NNS.dep}.
#'
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{ https://www.amazon.com/dp/1490523995/ref=cm_sw_su_dp}
#' @examples
#' \dontrun{
#' set.seed(123)
#' x <- rnorm(100) ; y <- rnorm(100)
#' NNS.dep(x, y)
#'
#' ## Correlation / Dependence Matrix
#' x <- rnorm(100) ; y <- rnorm(100) ; z <- rnorm(100)
#' B <- cbind(x, y, z)
#' NNS.dep(B)
#' }
#' @export

NNS.dep = function(x,
                   y = NULL,
                   asym = FALSE,
                   p.value = FALSE,
                   print.map = FALSE){

  if(any(class(x)%in%c("tbl","data.table"))) x <- as.vector(unlist(x))

  if(sum(is.na(x)) > 0) stop("You have some missing values, please address.")

  if(p.value){
      y_p <- replicate(100, sample.int(length(y)))
      x <- cbind(x, y, matrix(y[y_p], ncol = dim(y_p)[2], byrow = F))
      y <- NULL
  }

  if(!is.null(y)){
    x <- as.numeric(x)
    l <- length(x)

    y <- as.numeric(y)
    obs <- max(10, l/5)

    # Define segments
    if(print.map) PART <- NNS.part(x, y, order = NULL, obs.req = obs, min.obs.stop = TRUE, type = "XONLY", Voronoi = TRUE) else PART <- NNS.part(x, y, order = NULL, obs.req = obs, min.obs.stop = TRUE, type = "XONLY", Voronoi = FALSE)

    if(dim(PART$regression.points)[1]==0) return(list("Correlation" = 0, "Dependence" = 0))

    PART <- PART$dt
    PART <- PART[complete.cases(PART),]

    PART[, weights := .N/l, by = prior.quadrant]
    weights <- PART[, weights[1], by = prior.quadrant]$V1

    ll <- expression(max(min(100, .N), 8))
    res <- suppressWarnings(tryCatch(PART[,  sign(cor(x[1:eval(ll)],y[1:eval(ll)]))*summary(lm(y[1:eval(ll)]~poly(x[1:eval(ll)], max(1, min(10, as.integer(sqrt(.N))-1)), raw = TRUE)))$r.squared, by = prior.quadrant],
                                     error = function(e) PART[, NNS.copula(cbind(x,y)), by = prior.quadrant]))

    if(sum(is.na(res))>0) res[is.na(res)] <- NNS.copula(cbind(x,y))

    # Compare each asymmetry
    res_xy <- suppressWarnings(tryCatch(PART[,  sign(cor(x[1:eval(ll)],(y[1:eval(ll)])))*summary(lm(abs(y[1:eval(ll)])~poly(x[1:eval(ll)], max(1, min(10, as.integer(sqrt(.N))-1)), raw = TRUE)))$r.squared, by = prior.quadrant],
                                        error = function(e) PART[, NNS.copula(cbind(x,y)), by = prior.quadrant]))
    res_yx <- suppressWarnings(tryCatch(PART[,  sign(cor(y[1:eval(ll)],(x[1:eval(ll)])))*summary(lm(abs(x[1:eval(ll)])~poly(y[1:eval(ll)], max(1, min(10, as.integer(sqrt(.N))-1)), raw = TRUE)))$r.squared, by = prior.quadrant],
                                        error = function(e) PART[, NNS.copula(cbind(x,y)), by = prior.quadrant]))

    if(sum(is.na(res_xy))>0) res_xy[is.na(res_xy)] <- NNS.copula(cbind(x,y))
    if(sum(is.na(res_yx))>0) res_yx[is.na(res_yx)] <- NNS.copula(cbind(x,y))


    if(asym) dependence <- sum(abs(res_xy$V1) * weights) else dependence <- max(sum(abs(res$V1) * weights),
                                                                                sum(abs(res_xy$V1) * weights),
                                                                                sum(abs(res_yx$V1) * weights))

    lx <- PART[, length(unique(x))]
    ly <- PART[, length(unique(y))]
    degree_x <- min(10, max(1,lx-1), max(1,ly-1))

    I_x <- lx > sqrt(l)
    I_y <- ly > sqrt(l)
    I <- I_x * I_y

    poly_base <- dependence

    if(I == 1) poly_base <- suppressWarnings(tryCatch(summary(lm(abs(y)~poly(x, degree_x), raw = TRUE))$r.squared,
                                                      warning = function(w) dependence,
                                                      error = function(e) dependence))

    dependence <- mean(c(rep(dependence,3), poly_base))

    corr <- mean(c(sum(res$V1 * weights),
                   sum(res_xy$V1 * weights),
                   sum(res_yx$V1 * weights)))


    return(list("Correlation" = corr,
                "Dependence" = dependence))

  } else {
    if(p.value){
      original.par <- par(no.readonly = TRUE)

      nns.mc <- apply(x, 2, function(g) NNS.dep(x[,1], g))

      ## Store results
      cors <- unlist(lapply(nns.mc, "[[", 1))
      deps <- unlist(lapply(nns.mc, "[[", 2))


          cor_lower_CI <- LPM.VaR(.025, 0, cors[-c(1,2)])
          cor_upper_CI <- UPM.VaR(.025, 0, cors[-c(1,2)])
          dep_lower_CI <- LPM.VaR(.025, 0, deps[-c(1,2)])
          dep_upper_CI <- UPM.VaR(.025, 0, deps[-c(1,2)])
      if(print.map){
          par(mfrow = c(1, 2))
          hist(cors[-c(1,2)], main = "NNS Correlation", xlab = NULL, xlim = c(min(cors), max(cors[-1])))
          abline(v = cors[2], col = "red", lwd = 2)
          mtext("Result", side = 3, col = "red", at = cors[2])
          abline(v =  cor_lower_CI, col = "red", lwd = 2, lty = 3)
          abline(v =  cor_upper_CI , col = "red", lwd = 2, lty = 3)
          hist(deps[-c(1,2)], main = "NNS Dependence", xlab = NULL, xlim = c(min(deps), max(deps[-1])))
          abline(v = deps[2], col = "red", lwd = 2)
          mtext("Result", side = 3, col = "red", at = deps[2])
          abline(v =  dep_lower_CI , col = "red", lwd = 2, lty = 3)
          abline(v =  dep_upper_CI , col = "red", lwd = 2, lty = 3)
          par(mfrow = original.par)
      }

      return(list("Correlation" = as.numeric((cors)[2]),
                  "Correlation p.value" = min(LPM(0, cors[2], cors[-c(1,2)]),
                                              UPM(0, cors[2], cors[-c(1,2)])),
                  "Correlation 95% CIs" = c(cor_lower_CI, cor_upper_CI),
                  "Dependence" = as.numeric((deps)[2]),
                  "Dependence p.value" = min(LPM(0, deps[2], deps[-c(1,2)]),
                                             UPM(0, deps[2], deps[-c(1,2)])),
                  "Dependence 95% CIs" = c(dep_lower_CI, dep_upper_CI)))
    } else return(NNS.dep.matrix(x, asym = asym))
  }

}
