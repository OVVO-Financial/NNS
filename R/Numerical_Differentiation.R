#' NNS Numerical Differentiation
#'
#' Determines numerical derivative of a given function using projected secant lines on the y-axis.  These projected points infer finite steps \code{h}, in the finite step method.
#'
#' @param f an expression or call or a formula with no lhs.
#' @param point numeric; Point to be evaluated for derivative of a given function \code{f}.
#' @param h numeric [0, ...]; Initial step for secant projection.  Defaults to \code{(h = 0.1)}.
#' @param tol numeric; Sets the tolerance for the stopping condition of the inferred \code{h}.  Defualts to \code{(tol = 1e-10)}.
#' @param print.trace logical; \code{FALSE} (default) Displays each iteration, lower y-intercept, upper y-intercept and inferred \code{h}.
#' @return Returns a matrix of values, intercepts, derivatives, inferred step sizes for multiple methods of estimation.
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{https://www.amazon.com/dp/1490523995}
#' @examples
#' f <- function(x) sin(x) / x
#' NNS.diff(f, 4.1)
#'
#' @export


NNS.diff<- function(f, point, h = 0.1, tol = 1e-10, print.trace = FALSE){


  Finite.step <- function(f, point, h){

    f.x <- f(point)
    f.x.h.min <- f(point - h)
    f.x.h.pos <- f(point + h)

    neg.step <- (f.x - f.x.h.min) / h
    pos.step <- (f.x.h.pos - f.x) / h

    return((c("f(x-h)" = neg.step,
              "f(x+h)" = pos.step,
              mean(c(neg.step, pos.step)))))

  }




  orig.dig <- getOption("digits")
  options(digits = 20)
  Bs <- numeric()
  Bl <- numeric()
  Bu <- numeric()
  ### Step 1 initalize the boundaries for B

  ### Initial step size h
  f.x <- f(point)
  f.x.h <- f(point - h)

  ### Y = mX + B
  Y <- f.x

  m <- (f.x - f.x.h)/h

  mX <- ((f.x - f.x.h) / h) * point

  B <- f.x - ((f.x - f.x.h) / h) * point



  ### Initial interval for B given inputted h-step-value

  f.x.h.lower <- f(point - h)
  f.x.h.upper <- f(point + h)


  B1 <- f.x - ((f.x - f.x.h.lower) / h) * point
  B2 <- f.x - ((f.x.h.upper-f.x) / h) * point



  low.B <- min(c(B1, B2))
  high.B <- max(c(B1, B2))

  lower.B <- low.B
  upper.B <- high.B


  ## Return "Derivative Does Not Exist" if lower.B and upper.B are identical to 20 digits
  if(lower.B == upper.B){
      original.par=par(no.readonly = TRUE)
      par(mfrow = c(1, 2))
      plot(f, xlim = c(point - (100 * h), point + (100 * h)), col = 'blue', ylab = 'f(x)')
      points(point, f.x, pch = 19, col = 'red')
      plot(f, xlim = c(point - 1, point + 1), col = 'blue', ylab = 'f(x)')
      points(point, f.x, pch = 19, col = 'red')
      par(original.par)

      return(c("Derivative Does Not Exist"))}


      new.B <- mean(c(lower.B, upper.B))

      i <- 1L

      while(i >= 1L){
          Bl[i] <- lower.B
          Bu[i] <- upper.B
          new.f <- function(x) - f.x + ((f.x - f(point - x)) / x) * point + new.B

          ###  SOLVE FOR h, we just need the negative or positive sign from the tested B

          inferred.h <- uniroot(new.f, c(-2 * h, 2 * h), extendInt = 'yes')$root

          if(print.trace) {print(c("Iteration" = as.integer(i), "h" = inferred.h, "Lower B" = lower.B, "Upper B" = upper.B))}

          Bs[i] <- new.B

          ## Stop when the inferred h is within the tolerance level
          if(abs(inferred.h) < tol) {
              final.B <- mean(c(upper.B, lower.B))
              slope <- solve(point, f.x - final.B)

              z <- complex(real = point, imaginary = inferred.h)

              original.par <- par(no.readonly = TRUE)
              par(mfrow=c(1, 3))

              ## Plot #1
              plot(f, xlim = c(min(c(point - (100 * h), point + (100 * h)), 0), max(c(point - (100 * h), point + (100 * h)), 0)), col = 'azure4', ylab = 'f(x)', lwd = 2, ylim = c(min(c(min(c(B1, B2)), min(na.omit(f((point - (100 * h)) : (point + (100 * h))))))), max(c(max(na.omit(f((point - (100 * h)) : (point + (100 * h))))), max(c(B1, B2))))), main = 'f(x) and initial y-intercept range')
              abline(h = 0, v = 0, col = 'grey')
              points(point, f.x, pch = 19, col = 'green')
              points(point - h, f.x.h.lower, col = ifelse(B1 == high.B, 'blue', 'red'), pch = 19)
              points(point + h, f.x.h.upper, col = ifelse(B1 == high.B, 'red', 'blue'), pch = 19)
              points(x = rep(0, 2), y = c(B1, B2), col = c(ifelse(B1 == high.B, 'blue', 'red'), ifelse(B1 == high.B, 'red', 'blue')), pch = 1)
              segments(0, B1, point - h, f.x.h.lower, col = ifelse(B1 == high.B, 'blue','red'), lty = 2)
              segments(0, B2, point + h, f.x.h.upper, col = ifelse(B1 == high.B, 'red','blue'), lty = 2)

              ## Plot #2
              plot(f, col = 'azure4', ylab = 'f(x)', lwd = 3, main = 'f(x) narrowed range and secant lines', xlim = c(min(c(point - h, point + h,  0)), max(c(point + h,point - h, 0))), ylim= c(min(c(B1, B2, f.x.h.lower, f.x.h.upper)), max(c(B1, B2, f.x.h.lower, f.x.h.upper))))

              abline(h = 0, v = 0, col = 'grey')
              points(point,f.x, pch = 19, col = 'red')
              points(point - h, f.x.h.lower, col = ifelse(B1 == high.B, 'blue', 'red'), pch = 19)
              points(point + h, f.x.h.upper, col = ifelse(B1 == high.B, 'red', 'blue'), pch = 19)
              points(point, f.x, pch = 19, col = 'green')
              segments(0, B1, point - h, f.x.h.lower, col = ifelse(B1 == high.B, 'blue', 'red'), lty = 2)
              segments(0, B2, point + h, f.x.h.upper, col = ifelse(B1 == high.B, 'red', 'blue'), lty = 2)
              points(x = rep(0, 2), y = c(B1, B2), col = c(ifelse(B1 == high.B, 'blue', 'red'), ifelse(B1 == high.B, 'red', 'blue')), pch = 1)


              ## Plot #3
              plot(Bs, ylim = c(min(c(Bl, Bu)), max(c(Bl, Bu))), xlab = "Iterations", ylab = "y-inetercept", col = 'green', pch = 19, main = 'Iterated range of y-intercept')
              points(Bl, col = 'red', ylab = '')
              points(Bu, col = 'blue', ylab = '')

              legend('topright', c("Upper y-intercept", "Lower y-intercept", "Mean y-intercept"), col = c('blue', 'red', 'green'), pch = c(1, 1, 19), bty = 'n')

              par(original.par)

              return(as.matrix(c("Value of f(x) at point" = f(point),
                                "Final y-intercept (B)" = final.B,
                                "DERIVATIVE" = slope,
                                "Inferred h" = inferred.h,
                                "iterations" = i,
                                Finite.step(f, point, h)[1 : 2],
                                "Averaged Finite Step Initial h " = Finite.step(f, point, h)[3],
                                "Inferred h" = Finite.step(f, point, inferred.h)[1 : 2],
                                "Inferred h Averaged Finite Step" = Finite.step(f, point, inferred.h)[3],
                                "Complex Step Derivative (Initial h)" = Im(f(z)) / Im(z))))

          }


    ## NARROW THE RANGE OF B BASED ON SIGN OF INFERRED.H
    if(B1 == high.B){
        if(sign(inferred.h) < 0) {
            lower.B <- new.B
            upper.B <- upper.B
        } else {
            upper.B <- new.B
            lower.B <- lower.B
        }
    } else {
        if(sign(inferred.h) < 0) {
            lower.B <- lower.B
            upper.B <- new.B
        } else {
            upper.B <- upper.B
            lower.B <- new.B
        }
    }


    new.B = mean(c(lower.B, upper.B))

    i = i + 1
  }
  options(digits = orig.dig)
}


