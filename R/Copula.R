#' NNS Co-Partial Moments Higher Dimension Dependence
#'
#' Determines higher dimension dependence coefficients based on co-partial moment matrices ratios.
#'
#' @param X a numeric matrix or data frame.
#' @param continuous logical; \code{TRUE} (default) Generates a continuous measure using degree 1 \link{PM.matrix}, while discrete \code{FALSE} uses degree 0 \link{PM.matrix}.
#' @param plot logical; \code{FALSE} (default) Generates a 3d scatter plot with regression points using \link{plot3d}.
#' @param independence.overlay logical; \code{FALSE} (default) Creates and overlays independent \link{Co.LPM} and \link{Co.UPM} regions to visually reference the difference in dependence from the data.frame of variables being analyzed.  Under independence, the light green and red shaded areas would be occupied by green and red data points respectively.
#'
#' @return Returns a multivariate dependence value [0,1].
#'
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. (2016) "Beyond Correlation: Using the Elements of Variance for Conditional Means and Probabilities"  \url{https://www.ssrn.com/abstract=2745308}.
#' @examples
#' set.seed(123)
#' x <- rnorm(1000) ; y <- rnorm(1000) ; z <- rnorm(1000)
#' A <- data.frame(x, y, z)
#' NNS.copula(A, plot = TRUE, independence.overlay = TRUE)
#' @export


NNS.copula <- function (X,
                        continuous = TRUE,
                        plot = FALSE,
                        independence.overlay = FALSE){

    if(sum(is.na(X)) > 0) stop("You have some missing values, please address.")

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

    if(continuous) degree <- 1 else degree <- 0

    # Generate partial moment matrices
    pm_cov <- PM.matrix(degree, degree, variable = X, pop.adj = TRUE)

    # Isolate the upper triangles from each of the partial moment matrices
    Co_pm <- sum(pm_cov$cupm[upper.tri(pm_cov$cupm, diag = FALSE)]) + sum(pm_cov$clpm[upper.tri(pm_cov$clpm, diag = FALSE)])
    D_pm <- sum(pm_cov$dupm[upper.tri(pm_cov$dupm, diag = FALSE)]) + sum(pm_cov$dlpm[upper.tri(pm_cov$dlpm, diag = FALSE)])

    if(plot && n == 3){

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

    if(is.na(Co_pm) || is.null(Co_pm)) Co_pm <- 0
    if(is.na(D_pm)|| is.null(D_pm)) D_pm <- 0

    if(Co_pm == D_pm) return(0)
    if(Co_pm==0 || D_pm==0) return(1)
    if(Co_pm < D_pm) return(1 - Co_pm/D_pm)
    if(Co_pm > D_pm) return(1 - D_pm/Co_pm)


}
