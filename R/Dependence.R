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
#' For asymmetrical \code{(asym = TRUE)} matrices, directional dependence is returned as ([column variable] ---> [row variable]).
#'
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments" (ISBN: 1490523995)
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
  
  
  
  if (any(class(x) %in% c("tbl","data.table")) && !is.null(y)) {
    x <- if (!is.null(ncol(x)) && ncol(x) == 1) as.vector(unlist(x)) else as.numeric(x)
  }
  if (any(class(y) %in% c("tbl","data.table"))) {
    y <- if (!is.null(ncol(y)) && ncol(y) == 1) as.vector(unlist(y)) else as.numeric(y)
  }
  
  if(anyNA(x)) stop("You have some missing values, please address.")
  
  if(p.value){
    y_p <- replicate(100, sample.int(length(y)))
    x <- cbind(x, y, matrix(y[y_p], ncol = dim(y_p)[2], byrow = F))
    y <- NULL
  }
  
  if(!is.null(y)){
    l <- length(x)
    
    obs <- max(8, l/8)
    
    # Define segments
    if(print.map) PART_xy <- suppressWarnings(NNS.part(x, y, order = NULL, obs.req = obs, min.obs.stop = FALSE, type = NULL, Voronoi = TRUE)) else PART_xy <- suppressWarnings(NNS.part(x, y, order = NULL, obs.req = obs, min.obs.stop = FALSE, type = "XONLY", Voronoi = FALSE))
    
    PART_yx <- suppressWarnings(NNS.part(y, x, order = NULL, obs.req = obs, min.obs.stop = FALSE, type = "XONLY", Voronoi = FALSE))
    
    if(dim(PART_xy$regression.points)[1]==0) return(list("Correlation" = 0, "Dependence" = 0))
    
    PART_xy <- PART_xy$dt
    PART_xy <- PART_xy[complete.cases(PART_xy),]
    weights_xy <- PART_xy[, .N / l, by = quadrant]$V1
   
    PART_yx <- PART_yx$dt
    PART_yx <- PART_yx[complete.cases(PART_yx),]
    weights_yx <- PART_yx[, .N / l, by = quadrant]$V1

    
    dep_fn <- function(xx, yy) {
      NNS::NNS.copula(cbind(xx, yy)) * sign(fast_lm(xx, yy)$coef[2])
    }
  
    
    res_xy <- suppressWarnings(tryCatch(PART_xy[ , dep_fn(.SD[[1]], .SD[[2]]), by = quadrant, .SDcols = c(1,2)],
                                         error = function(e) PART_xy[ , dep_fn(.SD[[1]], .SD[[2]]), by = prior.quadrant, .SDcols = c(1,2)]))
     
    res_yx <- suppressWarnings(tryCatch(PART_yx[ , dep_fn(.SD[[1]], .SD[[2]]), by = quadrant, .SDcols = c(1,2)],
                                         error = function(e) PART_yx[ , dep_fn(.SD[[1]], .SD[[2]]), by = prior.quadrant, .SDcols = c(1,2)]))
    
    if(anyNA(res_xy)) res_xy[is.na(res_xy)] <- dep_fn(x, y)
    if(is.null(ncol(res_xy))) res_xy <- cbind(res_xy, res_xy)
    
    if(anyNA(res_yx)) res_yx[is.na(res_yx)] <- dep_fn(x, y)
    if(is.null(ncol(res_yx))) res_yx <- cbind(res_yx, res_yx)
    
    if(asym){
      dependence <- sum(abs(res_xy[,2]) * weights_xy)
    } else {
      dependence <- max(c(sum(abs(res_yx[,2]) * weights_yx),
                          sum(abs(res_xy[,2]) * weights_xy)))
    }
    
    lx <- PART_xy[, length(unique(x))]
    ly <- PART_xy[, length(unique(y))]
    degree_x <- min(10, max(1,lx-1), max(1,ly-1))
    
    I_x <- lx < sqrt(l)
    I_y <- ly < sqrt(l)
    I <- I_x * I_y
    
    if(I == 1){
      poly_base <- suppressWarnings(tryCatch(fast_lm_mult(poly(x, degree_x), abs(y))$r.squared,
                                             warning = function(w) dependence,
                                             error = function(e) dependence))
      
      dependence <- gravity(c(dependence, NNS.copula(cbind(x, y), plot = FALSE), poly_base))
    }
    
    if(asym){
      corr <- sum(res_xy[,2] * weights_xy)
    } else {
      corr <- max(c(sum(res_yx[,2] * weights_yx), sum(res_xy[,2] * weights_xy)))
    }
    
    
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