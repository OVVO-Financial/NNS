#' NNS Seasonality Test
#'
#' Seasonality test based on the coefficient of variance for the variable and lagged component series.  A result of 1 signifies no seasonality present.
#'
#' @param variable a numeric vector.
#' @param plot logical; \code{TRUE} (default) Returns the plot of all periods exhibiting seasonality and the variable level reference.
#' @return Returns a matrix of all periods exhibiting less coefficient of variance than the variable with \code{"all.periods"}; and the single period exhibiting the least coefficient of variance versus the variable with \code{"best.period"}.  If no seasonality is detected, \code{NNS.seas} will return ("No Seasonality Detected").
#' @keywords seasonality
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' set.seed(123)
#' x <- rnorm(100)
#'
#' ## To call strongest period based on coefficient of variance:
#' NNS.seas(x)$best.period
#' @export


NNS.seas <- function(variable, plot = TRUE){

  variable_1 = variable[1 : (length(variable) - 1)]
  variable_2 = variable_1[1 : (length(variable_1) - 1)]


  output = numeric() ; output_1 = numeric() ; output_2 = numeric()
  instances = numeric() ; instances_1 = numeric() ; instances_2 = numeric()

  if(mean(variable) != 0){
  var.cov = abs(sd(variable) / mean(variable))
  } else {
    var.cov = Inf
  }

  for(i in 1 : (length(variable) / 4)){
    reverse.var = variable[seq(length(variable), 1, -i)]
    reverse.var_1 = variable_1[seq(length(variable_1), 1, -i)]
    reverse.var_2 = variable_2[seq(length(variable_2), 1, -i)]

    test = abs(sd(reverse.var) / mean(reverse.var))
    test_1 = abs(sd(reverse.var_1) / mean(reverse.var_1))
    test_2 = abs(sd(reverse.var_2) / mean(reverse.var_2))
      if (test <= var.cov){
          instances[i] <- i
          output[i] <- test
      } else {
            instances[i] <- 0
            output[i] <- 0
      }

    if (test_1 <= var.cov){
      instances_1[i] <- i
      output_1[i] <- test_1
    } else {
      instances_1[i] <- 0
      output_1[i] <- 0
    }

    if (test_2 <= var.cov){
      instances_2[i] <- i
      output_2[i] <- test_2
    } else {
      instances_2[i] <- 0
      output_2[i] <- 0
    }
  }

  ref.output = cbind(instances, output, output_1, output_2, output * output_1 * output_2 > 0)
  output = rowMeans(ref.output[ , 2 : 4]) * ref.output[ , 5]

  instances = ref.output[ , 1] * ref.output[ , 5]

  index = which(instances > 0 & output > 0)

  insts = sum(instances > 0) > 0

  if(insts){
      n = rep(var.cov, length(instances[index]))

      M = data.table("Period" = instances[index], "Coefficient.of.Variance" = output[index], "Variable.Coefficient.of.Variance" = n, key = "Coefficient.of.Variance")
  } else {
      M ="No Seasonality Detected"
  }


  if(insts){
      if(plot){
        plot(instances[index], output[index],
         xlab = "Period", ylab = "Coefficient of Variance", main = "Seasonality Test",
         ylim = c(0, 2 * abs(sd(variable) / mean(variable))))

    points(M[1, Period], M[1, Coefficient.of.Variance], pch = 19, col = 'red')

    abline(h = abs(sd(variable) / mean(variable)), col = "red", lty = 5)
    text(mean(instances[index]), abs(sd(variable) / mean(variable)), pos = 3,
         "Variable Coefficient of Variance", col = 'red')}

    return(list("all.periods" = M,
                "best.period" = M[1, Period]))
    } else {
      return(M)
    }

}
