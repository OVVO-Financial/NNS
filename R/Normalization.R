#' NNS Normalization
#'
#' Normalizes a matrix of variables based on nonlinear scaling normalization method.
#'
#' @param X a numeric matrix or data frame.
#' @param linear logical; \code{FALSE} (default) Performs a linear scaling normalization, resulting in equal means for all variables.
#' @param chart.type  options: ("l", "b"); \code{NULL} (default).  Set \code{(chart.type = "l")} for line,
#' \code{(chart.type = "b")} for boxplot.
#' @param location Sets the legend location within the plot, per the \code{x} and \code{y} co-ordinates used in base graphics \link{legend}.
#' @return Returns a \link{data.frame} of normalized values.
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{https://www.amazon.com/dp/1490523995/ref=cm_sw_su_dp}
#' @examples
#' set.seed(123)
#' x <- rnorm(100) ; y<-rnorm(100)
#' A <- cbind(x, y)
#' NNS.norm(A)
#' @export

NNS.norm <- function(X,
                     linear = FALSE,
                     chart.type = NULL,
                     location = "topleft"){

  if(sum(is.na(X)) > 0) stop("You have some missing values, please address.")

  if(any(class(X)%in%c("tbl","data.table"))) X <- as.data.frame(X)

  X <- apply(X, 2, unlist)

  m  <- colMeans(X)
  m[m==0] <- 1e-10
  RG <- m %o% (1 / m)

  if(!linear){
      if(length(m) < 10){
        scale.factor <- abs(cor(X))
      } else {
        scale.factor <- abs(NNS.dep(X)$Dependence)
      }
    scales <- Rfast::colmeans(RG * scale.factor)
  } else {
      scales <- Rfast::colmeans(RG)
  }


  X_Normalized <- t(t(X) * scales)

  n <- ncol(X)
  i <- seq_len(n)

  if(is.null(colnames(X))){
      new.names <- list()
      for(i in 1 : n){
          new.names[[i]] <- paste0("x_", i)
      }
      colnames(X) <- unlist(new.names)
  }

  labels <- c(colnames(X), paste0(colnames(X), " Normalized"))

  colnames(X_Normalized) <- labels[(n + 1) : (2 * n)]

  rows <- rownames(X_Normalized)

if(!is.null(chart.type)){
    original.par <- par(no.readonly = TRUE)
        if(chart.type == 'b' ){
            par(mar = c(10, 4, 3, 2) + 0.1)
            boxplot(cbind(X, X_Normalized), las = 2, names = labels, col = c(rep("grey", n), rainbow(n)))
        }

        if(chart.type == 'l' ){
            par(mfrow = c(2, 1))
            par(mar = c(ifelse((class(rows)!="numeric" || !is.null(rows)),5,2), nchar(max(X_Normalized))*1/exp(1), 1, 2))


            matplot(X, type = 'l', col = c('steelblue', rainbow(n)), ylab = '', xaxt = 'n', lwd = 2, las = 1)
            legend(location, inset = c(0,0), c(colnames(X)), lty = 1, col = c('steelblue', rainbow(n)), bty = 'n', ncol = floor(n/sqrt(n)), lwd = 2, cex = n/sqrt(n)^exp(1))
            axis(1, at = seq(length(X_Normalized[ , 1]), 1, -floor(sqrt(length(X_Normalized[ , 1])))),
                 labels = rownames(X_Normalized[seq(length(X_Normalized[ , 1]), 1, -floor(sqrt(length(X_Normalized[ , 1])))),]), las = 1,
                 cex.axis = ifelse((class(rows)!="numeric" || !is.null(rows)),.75,1),
                 las = ifelse((class(rows)!="numeric" || !is.null(rows)),3,1),srt=45)

            matplot(X_Normalized, type = 'l', col = c('steelblue', rainbow(n)), ylab = '', xaxt = 'n', lwd = 2, las = 1)
            axis(1, at = seq(length(X_Normalized[ , 1]), 1, -floor(sqrt(length(X_Normalized[ , 1])))),
                 labels = rownames(X_Normalized[seq(length(X_Normalized[ , 1]), 1, -floor(sqrt(length(X_Normalized[ , 1])))),]), las = 1,
                 cex.axis = ifelse((class(rows)!="numeric" || !is.null(rows)),.75,1),
                 las = ifelse((class(rows)!="numeric" || !is.null(rows)),3,1),srt=45)

            legend(location, c(paste0(colnames(X), " Normalized")), lty = 1, col = c('steelblue', rainbow(n)), bty = 'n', ncol = ceiling(n/sqrt(n)), lwd = 2, cex = n/sqrt(n)^exp(1))
        }

  par(original.par)

  }



  return(X_Normalized)

}
