#' NNS Normalization
#'
#' Normalizes a matrix of variables based on nonlinear scaling normalization method.
#' @param A a numeric matrix or data frame.
#' @param chart.type  options: ("l", "b"); \code{NULL} (default).  Set \code{(chart.type = "l")} for line,
#' \code{(chart.type = "b")} for boxplot.
#' @param linear logical; \code{FALSE} (default) Performs a linear scaling normalization, resulting in equal means for all variables.
#' @return Returns a \link{data.frame} of normalized values.
#' @keywords normalization
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' set.seed(123)
#' x <- rnorm(100) ; y<-rnorm(100)
#' A <- cbind(x, y)
#' NNS.norm(A)
#' @export

NNS.norm <- function(A, chart.type = NULL, linear = FALSE) {
  m  <- colMeans(A)
  m <- pmax(1e-10, m)
  RG <- m %o% (1 / m)

  if(!linear){
      scale.factor = abs(cor(A))
      scales <- colMeans(RG * scale.factor)
  } else {
        scales <- colMeans(RG)
    }


  A_Normalized <- t(t(A) * scales)

  n <- ncol(A)
  i <- seq_len(n)

  if(is.null(colnames(A))){
    new.names = list()
    for(i in 1 : n){
      new.names[[i]] = paste0("x_", i)
    }
    colnames(A) = unlist(new.names)
    }

  labels <- c(colnames(A), paste0(colnames(A), " Normalized"))

  colnames(A_Normalized) <- labels[(n + 1) : (2 * n)]

if(!is.null(chart.type)){
    if(chart.type == 'b' ){
        par(mar = c(10, 4, 3, 2) + 0.1)
        boxplot(cbind(A, A_Normalized),
          las = 2, names = labels,
          col = c(rep("grey", n), rainbow(n)))
    }

    if(chart.type == 'l' ){
        par(mar = c(3, 2, 2, 2))
        par(mfrow = c(2, 1))

        matplot(A,type = 'l', col = c('steelblue', rainbow(n)), ylab = '', xaxt = 'n', lwd = 2)
        legend('top', inset = c(0,0), c(colnames(A)), lty = 1, col = c('steelblue', rainbow(n)), bty = 'n',horiz = TRUE, lwd = 2)
        axis(1, at = seq(length(A[ , 1]), 1, -round(sqrt(length(A[ , 1])))), labels = rownames(A[seq(length(A[ , 1]), 1,-round(sqrt(length(A[ , 1])))),]),las = 1,cex.axis = 1)

        matplot(A_Normalized, type = 'l', col = c('steelblue', rainbow(n)), ylab = '', xaxt = 'n', lwd = 2)
        axis(1, at = seq(length(A[ , 1]), 1, -round(sqrt(length(A[ , 1])))), labels = rownames(A[seq(length(A[ , 1]), 1, -round(sqrt(length(A[ , 1])))),]), las = 1, cex.axis = 1)
        legend('top', c(paste0(colnames(A), " Normalized")), lty = 1, col = c('steelblue', rainbow(n)), bty = 'n', horiz = TRUE, lwd = 2)
    }
  }

  par(mfrow=c(1, 1))

  return(A_Normalized)

}
