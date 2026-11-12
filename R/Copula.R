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
#' @return Returns a multivariate dependence value [0,1]. Each co-partial moment is compared with its value under independence, the product of the marginal partial moments at \code{target}, so independent variables give a value near 0 for any target and marginal distributions.  At a target with every variable's observations on one side, the normalized directional-region shares carry no contrast and the neutral value is 0. This target degeneracy is not evidence that the variables are independent. In general, 0 and 1 describe this target-specific summary, not global independence or deterministic dependence.
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
  
  if((plot||independence.overlay) && n == 3){
    .nns_require_rgl()
    if(is.null(colnames(X))) colnames(X) <- paste0("Var ", seq_len(n))
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
  
  # The dependence itself is computed in one compiled call (src/NNS_copula.cpp):
  # the pairwise degree-0 and degree-1 co-partial moments about the target and
  # the n-dimensional divergent moment, each measured against its independence
  # anchor from the margins' own lower-tail shares.
  if(is.data.frame(X)) X <- as.matrix(X)
  if(!is.double(X)) storage.mode(X) <- "double"
  NNS_copula_cpp(X, as.numeric(target), if(continuous) TRUE else FALSE)
}