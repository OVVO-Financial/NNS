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
#' @return Returns a multivariate dependence value [0,1].
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
  
  # Pairwise  
  discrete_pm_cov <- PM.matrix(LPM_degree = 0, UPM_degree = 0, target = target, variable = X, pop_adj = FALSE)
  utr <- upper.tri(discrete_pm_cov$cupm, diag = FALSE)
  discrete_Co_pm <- sum(discrete_pm_cov$cupm[utr]) + sum(discrete_pm_cov$clpm[utr]) 
  if(discrete_Co_pm==1 || discrete_Co_pm==0) return(1)
  
  
  if(continuous){
    continuous_pm_cov <- PM.matrix(LPM_degree = 1, UPM_degree = 1, target = target, variable = X, pop_adj = TRUE, norm = TRUE)
  } else {
    continuous_pm_cov <- discrete_pm_cov
  }
  

  # Isolate the upper triangles from each of the partial moment matrices
  discrete_D_pm <- sum(discrete_pm_cov$dupm[utr]) + sum(discrete_pm_cov$dlpm[utr])
  
  continuous_Co_pm <- sum(continuous_pm_cov$cupm[utr]) + sum(continuous_pm_cov$clpm[utr]) 
  continuous_D_pm <- sum(continuous_pm_cov$dupm[utr]) + sum(continuous_pm_cov$dlpm[utr]) 
  
  indep_Co_pm <- .25 * (n^2 - n)
  
  discrete_dep <- abs(discrete_Co_pm-indep_Co_pm)/indep_Co_pm
  continuous_dep <- abs(continuous_Co_pm-indep_Co_pm)/indep_Co_pm 
  
   
  discrete_dep <- min(max(discrete_dep, 0), 1)
  continuous_dep <- min(max(continuous_dep, 0), 1)

  # n-dimensional
  discrete_D_pm <- DPM_nD(data = X, target = target, degree = 0, norm = TRUE)
  if(continuous) continuous_D_pm <- DPM_nD(data = X, target = target, degree = 1, norm = TRUE) else continuous_D_pm <- discrete_D_pm
  
  indep_D_pm <- 1-(0.5^n)
  
  n_dim_discrete_dep <- abs(discrete_D_pm - indep_D_pm)/indep_D_pm 
  n_dim_continuous_dep <- abs(continuous_D_pm - indep_D_pm)/indep_D_pm


  return(mean(c(discrete_dep, continuous_dep, n_dim_discrete_dep, n_dim_continuous_dep))^(1/2))
}