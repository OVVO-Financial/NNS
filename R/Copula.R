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
#' @references Viole, F. (2016) "Beyond Correlation: Using the Elements of Variance for Conditional Means and Probabilities"  \url{https://www.ssrn.com/abstract=2745308}.
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
  if(sum(is.na(X)) > 0){
	stop("You have some missing values, please address.")
  }

  n <- ncol(X)
  l <- dim(X)[1]

  if(any(class(X)%in%c("tbl","data.table"))) X <- as.data.frame(X)

  if(is.null(colnames(X))){
    colnames.list <- list()
    for(i in 1 : n){
      colnames.list[i] <- paste0("Var ", i)
    }
    colnames(X) <- c(colnames.list)
  }
  
  discrete_pm_cov <- PM.matrix(0, 0, target = target, variable = X, pop_adj = FALSE)

  if(continuous){
    degree <- 1
    continuous_pm_cov <- PM.matrix(degree, degree, target = target, variable = X, pop_adj = TRUE)
  } else {
    degree <- 0
    continuous_pm_cov <- discrete_pm_cov
  }
  

  # Isolate the upper triangles from each of the partial moment matrices
  continuous_Co_pm <- sum(continuous_pm_cov$cupm[upper.tri(continuous_pm_cov$cupm, diag = FALSE)]) + sum(continuous_pm_cov$clpm[upper.tri(continuous_pm_cov$clpm, diag = FALSE)])
  continuous_D_pm <- sum(continuous_pm_cov$dupm[upper.tri(continuous_pm_cov$dupm, diag = FALSE)]) + sum(continuous_pm_cov$dlpm[upper.tri(continuous_pm_cov$dlpm, diag = FALSE)])

  discrete_Co_pm <- sum(discrete_pm_cov$cupm[upper.tri(discrete_pm_cov$cupm, diag = FALSE)]) + sum(discrete_pm_cov$clpm[upper.tri(discrete_pm_cov$clpm, diag = FALSE)])
  discrete_D_pm <- sum(discrete_pm_cov$dupm[upper.tri(discrete_pm_cov$dupm, diag = FALSE)]) + sum(discrete_pm_cov$dlpm[upper.tri(discrete_pm_cov$dlpm, diag = FALSE)])

 
  indep_Co_pm <- .25 * (n^2 - n)

  if(discrete_Co_pm > indep_Co_pm) discrete_dep <- (discrete_Co_pm-indep_Co_pm)/indep_Co_pm else discrete_dep <- (indep_Co_pm - discrete_Co_pm)/indep_Co_pm
  discrete_dep <- min(max(sqrt(discrete_dep), 0), 1)
  
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

  if(is.na(continuous_Co_pm) || is.null(continuous_Co_pm)) continuous_Co_pm <- 0
  if(is.na(continuous_D_pm)|| is.null(continuous_D_pm)) continuous_D_pm <- 0

  if(continuous_Co_pm == continuous_D_pm) return(mean(c(0, discrete_dep)))
  if(continuous_Co_pm==0 || continuous_D_pm==0) return(mean(c(1, discrete_dep)))
  
 

  if(continuous_Co_pm < continuous_D_pm) return(mean(c((1 - (continuous_Co_pm/continuous_D_pm)), discrete_dep)))
  if(continuous_Co_pm > continuous_D_pm) return(mean(c((1 - (continuous_D_pm/continuous_Co_pm)), discrete_dep)))

}
