#' NNS Co-Partial Moments Higher Dimension Dependence
#'
#' Determines higher dimension dependence coefficients based on co-partial moment matrices ratios.
#'
#' @param X a numeric matrix or data frame.
#' @param target numeric; Typically the mean of Variable X for classical statistics equivalences, but does not have to be. (Vectorized)  \code{(target = NULL)} (default) will set the target as the mean of every variable.
#' @param continuous logical; \code{TRUE} (default) Generates a continuous measure using degree 1 \link{PM.matrix}, while discrete \code{FALSE} uses degree 0 \link{PM.matrix}.
#' @param plot logical; \code{FALSE} (default) Generates a 3d scatter plot with regression points.
#' @param independence.overlay logical; \code{FALSE} (default) Creates and overlays independent \link{Co.LPM} and \link{Co.UPM} regions to visually reference the difference in dependence from the data.frame of variables being analyzed.  Under independence, the light green and red shaded areas would be occupied by green and red data points respectively.
#'
#' @return Returns a multivariate dependence value [0,1]. Each co-partial moment is compared with its value under independence, the product of the marginal partial moments at \code{target}, so independent variables give a value near 0 for any target and marginal distributions.  At a target with every variable's observations on one side of it, each co-partial moment equals its independence value whatever the dependence, so the value is 0.
#'
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. (2016) "Beyond Correlation: Using the Elements of Variance for Conditional Means and Probabilities"  \doi{10.2139/ssrn.2745308}.
#' @examples
#' \dontrun{
#' set.seed(123)
#' x <- rnorm(1000) ; y <- rnorm(1000) ; z <- rnorm(1000)
#' A <- data.frame(x, y, z)
#' NNS.copula(A, target = colMeans(A), plot = TRUE, independence.overlay = TRUE)
#'
#' ### Target 0
#' NNS.copula(A, target = rep(0, ncol(A)), plot = TRUE, independence.overlay = TRUE)
#' }
#' @export


NNS.copula <- function (
    X,
    target = NULL,
    continuous = TRUE,
    plot = FALSE,
    independence.overlay = FALSE
){
  
  if(anyNA(X)) stop("You have some missing values, please address.")
  
  n <- ncol(X)
  
  if(any(class(X)%in%c("tbl","data.table"))) X <- as.data.frame(X)
  
  if(is.null(colnames(X))) colnames(X) <- paste0("Var ", seq_len(n))
  
  if((plot||independence.overlay) && n == 3){
    .nns_require_rgl()
    rgl::plot3d(x = X[ , 1], y = X[ , 2], z = X[ , 3], box = FALSE, size = 3,
                col=ifelse((X[ , 1] <= mean(X[ , 1])) & (X[ , 2] <= mean(X[ , 2])) & (X[ , 3] <= mean(X[ , 3])), 'red' ,
                           ifelse((X[ , 1] > mean(X[ , 1])) & (X[ , 2] > mean(X[ , 2])) & (X[ , 3] > mean(X[ , 3])), 'green',
                                  'steelblue')), xlab = colnames(X)[1], ylab = colnames(X)[2], zlab = colnames(X)[3])
    
    if(independence.overlay == TRUE){
      clpm.box <- rgl::cube3d(color = "red", alpha = 0.25)
      cupm.box <- rgl::cube3d(color = "green", alpha = 0.25)
      
      clpm.box$vb[1, ] <- replace(clpm.box$vb[1, ], clpm.box$vb[1, ] == -1, min(X[ , 1]))
      clpm.box$vb[2, ] <- replace(clpm.box$vb[2, ], clpm.box$vb[2, ] == -1, min(X[ , 2]))
      clpm.box$vb[3, ] <- replace(clpm.box$vb[3, ], clpm.box$vb[3, ] == -1, min(X[ , 3]))
      clpm.box$vb[1, ] <- replace(clpm.box$vb[1, ], clpm.box$vb[1, ] == 1, mean(X[, 1]))
      clpm.box$vb[2, ] <- replace(clpm.box$vb[2, ], clpm.box$vb[2, ] == 1, mean(X[, 2]))
      clpm.box$vb[3, ] <- replace(clpm.box$vb[3, ], clpm.box$vb[3, ] == 1, mean(X[, 3]))
      
      cupm.box$vb[1, ] <- replace(cupm.box$vb[1, ], cupm.box$vb[1, ] == 1, max(X[ , 1]))
      cupm.box$vb[2, ] <- replace(cupm.box$vb[2, ], cupm.box$vb[2, ] == 1, max(X[ , 2]))
      cupm.box$vb[3, ] <- replace(cupm.box$vb[3, ], cupm.box$vb[3, ] == 1, max(X[ , 3]))
      cupm.box$vb[1, ] <- replace(cupm.box$vb[1, ], cupm.box$vb[1, ] == -1, mean(X[, 1]))
      cupm.box$vb[2, ] <- replace(cupm.box$vb[2, ], cupm.box$vb[2, ] == -1, mean(X[, 2]))
      cupm.box$vb[3, ] <- replace(cupm.box$vb[3, ], cupm.box$vb[3, ] == -1, mean(X[, 3]))
      
      rgl::shade3d(clpm.box)
      rgl::shade3d(cupm.box)
    }
  }
  
  if(is.null(target)) target <- colMeans(X)
  
  # A constant margin is independent of the rest and has no degree-1 lower share
  vary <- vapply(seq_len(n), function(j) length(unique(X[ , j])) > 1L, logical(1))
  if(!all(vary)){
    if(sum(vary) < 2) return(0)
    X <- X[ , vary, drop = FALSE]
    target <- target[vary]
    n <- ncol(X)
  }
  
  # Pairwise  
  discrete_pm_cov <- PM.matrix(LPM_degree = 0, UPM_degree = 0, target = target, variable = X, pop_adj = FALSE)
  utr <- upper.tri(discrete_pm_cov$cupm, diag = FALSE)
  discrete_Co_pm <- sum(discrete_pm_cov$cupm[utr]) + sum(discrete_pm_cov$clpm[utr]) 
  # Perfect dependence: every pair fully concordant (1) or fully discordant (0).
  # A margin with every observation on one side of its target makes that trivially true.
  lower_share <- vapply(seq_len(n), function(j) LPM(0, target[j], X[ , j]), numeric(1))
  # Every co-partial moment then equals its independence product.
  if(all(lower_share == 0 | lower_share == 1)) return(0)
  pair_Co_pm <- discrete_pm_cov$cupm[utr] + discrete_pm_cov$clpm[utr]
  if(all(lower_share > 0 & lower_share < 1) &&
     all(abs(pair_Co_pm - 1) < 1e-9 | abs(pair_Co_pm) < 1e-9)) return(1)
  
  
  if(continuous){
    continuous_pm_cov <- PM.matrix(LPM_degree = 1, UPM_degree = 1, target = target, variable = X, pop_adj = TRUE, norm = TRUE)
  } else {
    continuous_pm_cov <- discrete_pm_cov
  }
  

  # Isolate the upper triangles from each of the partial moment matrices
  discrete_D_pm <- sum(discrete_pm_cov$dupm[utr]) + sum(discrete_pm_cov$dlpm[utr])
  
  continuous_Co_pm <- sum(continuous_pm_cov$cupm[utr]) + sum(continuous_pm_cov$clpm[utr]) 
  continuous_D_pm <- sum(continuous_pm_cov$dupm[utr]) + sum(continuous_pm_cov$dlpm[utr]) 
  
  # Independence anchors from each margin's own lower-tail share, so they hold at any target
  indep_anchors <- function(p) list(Co = sum((outer(p, p) + outer(1 - p, 1 - p))[utr]), D = 1 - prod(p) - prod(1 - p))
  discrete_indep <- indep_anchors(lower_share)
  continuous_indep <- if(continuous) indep_anchors(vapply(seq_len(n), function(j) LPM.ratio(1, target[j], X[ , j]), numeric(1))) else discrete_indep
  # A boundary target can zero an anchor; an observed 0 then matches it exactly
  rel_dev <- function(obs, indep) if(obs == indep) 0 else abs(obs - indep)/indep
  
  discrete_dep <- rel_dev(discrete_Co_pm, discrete_indep$Co)
  continuous_dep <- rel_dev(continuous_Co_pm, continuous_indep$Co)
  
   
  discrete_dep <- min(max(discrete_dep, 0), 1)
  continuous_dep <- min(max(continuous_dep, 0), 1)

  # n-dimensional
  discrete_D_pm <- DPM_nD(data = X, target = target, degree = 0, norm = TRUE)
  if(continuous) continuous_D_pm <- DPM_nD(data = X, target = target, degree = 1, norm = TRUE) else continuous_D_pm <- discrete_D_pm
  
  n_dim_discrete_dep <- rel_dev(discrete_D_pm, discrete_indep$D)
  n_dim_continuous_dep <- rel_dev(continuous_D_pm, continuous_indep$D)
  n_dim_discrete_dep <- min(n_dim_discrete_dep, 1)
  n_dim_continuous_dep <- min(n_dim_continuous_dep, 1)


  return(mean(c(discrete_dep, continuous_dep, n_dim_discrete_dep, n_dim_continuous_dep))^(1/2))
}