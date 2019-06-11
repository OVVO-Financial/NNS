#' NNS Co-Partial Moments Higher Dimension Correlation
#'
#' Determines higher dimension correlation coefficients based on degree 0 co-partial moments.
#'
#' @param x a numeric matrix or data frame.
#' @param plot logical; \code{FALSE} (default) Generates a 3d scatter plot with regression points using \link{plot3d}.
#' @param independence.overlay logical; \code{FALSE} (default) Creates and overlays independent \link{Co.LPM} and \link{Co.UPM} regions to visually reference the difference in dependence from the data.frame of variables being analyzed.  Under independence, the light green and red shaded areas would be occupied by green and red data points respectively.
#' @return Returns multivariate nonlinear correlation coefficient
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. (2016) "Beyond Correlation: Using the Elements of Variance for Conditional Means and Probabilities"  \url{http://ssrn.com/abstract=2745308}.
#' @examples
#' set.seed(123)
#' x <- rnorm(1000) ; y <- rnorm(1000) ; z <- rnorm(1000)
#' A <- data.frame(x, y, z)
#' NNS.cor.hd(A, plot = TRUE, independence.overlay = TRUE)
#' @export


NNS.cor.hd <- function (x, plot = FALSE, independence.overlay = FALSE){
    A <- x
    n <- ncol(A)
    l <- length(A[ , 1])

    if(is.null(colnames(A))){
        colnames.list <- list()
        for(i in 1 : n){
            colnames.list[i] <- paste0("Var ", i)
        }
        colnames(A) <- c(colnames.list)
    }

    A_upm <- t(apply(A, 1, function(x) x > colMeans(A)))
    A_lpm <- t(apply(A, 1, function(x) x <= colMeans(A)))


    CO_upm <- sum(apply(A_upm, 1, prod)) / l
    CO_lpm <- sum(apply(A_lpm, 1, prod)) / l


    observed <- CO_upm + CO_lpm
    independence <- 2 * (.5 ^ n)

    if(plot && n == 3){

        plot3d(x = A[ , 1], y = A[ , 2], z = A[ , 3], box = FALSE, size = 3,
           col=ifelse((A[ , 1] <= mean(A[ , 1])) & (A[ , 2] <= mean(A[ , 2])) & (A[ , 3] <= mean(A[ , 3])), 'red' ,
                      ifelse((A[ , 1] > mean(A[ , 1])) & (A[ , 2] > mean(A[ , 2])) & (A[ , 3] > mean(A[ , 3])), 'green',
                      'steelblue')), xlab = colnames(A)[1], ylab = colnames(A)[2], zlab = colnames(A)[3])

        if(independence.overlay == TRUE){

            clpm.box <- cube3d(color = "red", alpha = 0.25)
            cupm.box <- cube3d(color = "green", alpha = 0.25)

            clpm.box$vb[1, ] <- replace(clpm.box$vb[1, ], clpm.box$vb[1, ] == -1, min(A[ , 1]))
            clpm.box$vb[2, ] <- replace(clpm.box$vb[2, ], clpm.box$vb[2, ] == -1, min(A[ , 2]))
            clpm.box$vb[3, ] <- replace(clpm.box$vb[3, ], clpm.box$vb[3, ] == -1, min(A[ , 3]))
            clpm.box$vb[1, ] <- replace(clpm.box$vb[1, ], clpm.box$vb[1, ] == 1, (max(A[ , 1]) + min(A[ , 1])) / 2)
            clpm.box$vb[2, ] <- replace(clpm.box$vb[2, ], clpm.box$vb[2, ] == 1, (max(A[ , 2]) + min(A[ , 2])) / 2)
            clpm.box$vb[3, ] <- replace(clpm.box$vb[3, ], clpm.box$vb[3, ] == 1, (max(A[ , 3]) + min(A[ , 3])) / 2)

            cupm.box$vb[1, ] <- replace(cupm.box$vb[1, ], cupm.box$vb[1, ] == 1, max(A[ , 1]))
            cupm.box$vb[2, ] <- replace(cupm.box$vb[2, ], cupm.box$vb[2, ] == 1, max(A[ , 2]))
            cupm.box$vb[3, ] <- replace(cupm.box$vb[3, ], cupm.box$vb[3, ] == 1, max(A[ , 3]))
            cupm.box$vb[1, ] <- replace(cupm.box$vb[1, ], cupm.box$vb[1, ] == -1, (max(A[ , 1]) + min(A[ , 1])) / 2)
            cupm.box$vb[2, ] <- replace(cupm.box$vb[2, ], cupm.box$vb[2, ] == -1, (max(A[ , 2]) + min(A[ , 2])) / 2)
            cupm.box$vb[3, ] <- replace(cupm.box$vb[3, ], cupm.box$vb[3, ] == -1, (max(A[ , 3]) + min(A[ , 3])) / 2)


            shade3d(clpm.box)
            shade3d(cupm.box)


        }

    }

    return((observed - independence) / (1 - independence))
}
