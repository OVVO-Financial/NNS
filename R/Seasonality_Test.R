#' NNS Seasonality Test
#'
#' Seasonality test based on the coefficient of variation for the variable and lagged component series.  A result of 1 signifies no seasonality present.
#'
#' @param variable a numeric vector.
#' @param modulo integer(s); NULL (default) Used to find the nearest multiple(s) in the reported seasonal period.
#' @param mod.only logical; code{TRUE} (default) Limits the number of seasonal periods returned to the specified \code{modulo}.
#' @param plot logical; \code{TRUE} (default) Returns the plot of all periods exhibiting seasonality and the variable level reference.
#' @return Returns a matrix of all periods exhibiting less coefficient of variation than the variable with \code{"all.periods"}; and the single period exhibiting the least coefficient of variation versus the variable with \code{"best.period"}; as well as a vector of \code{"periods"} for easy call into \link{NNS.ARMA.optim}.  If no seasonality is detected, \code{NNS.seas} will return ("No Seasonality Detected").
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{https://www.amazon.com/dp/1490523995/ref=cm_sw_su_dp}
#' @examples
#' set.seed(123)
#' x <- rnorm(100)
#'
#' ## To call strongest period based on coefficient of variation:
#' NNS.seas(x, plot = FALSE)$best.period
#'
#' ## Using modulos for logical seasonal inference:
#' NNS.seas(x, modulo = c(2,3,5,7), plot = FALSE)
#' @export


NNS.seas <- function(variable,
                     modulo = NULL,
                     mod.only = TRUE,
                     plot = TRUE){

  if(any(class(variable)=="tbl")) variable <- as.vector(unlist(variable))

  if(length(variable) < 5){
    return(data.table::data.table("Period" = 0, "Coefficient.of.Variation" = 0, "Variable.Coefficient.of.Variation" = 0, key = "Coefficient.of.Variation"))
  }

  if(is.null(modulo)){
      mod.only <- FALSE
  }

  variable_1 <- variable[1 : (length(variable) - 1)]
  variable_2 <- variable_1[1 : (length(variable_1) - 1)]


  output <- numeric() ; output_1 = numeric() ; output_2 = numeric()
  instances <- numeric() ; instances_1 = numeric() ; instances_2 = numeric()

  if(mean(variable) != 0){
    var.cov <- abs(sd(variable) / mean(variable))
  } else {
    var.cov <- abs(acf(variable, lag.max = 1, plot = FALSE)$acf[2])^-1
  }

  for(i in 1 : (length(variable) / 2)){
    reverse.var <- variable[seq(length(variable), 1, -i)]
    reverse.var_1 <- variable_1[seq(length(variable_1), 1, -i)]
    reverse.var_2 <- variable_2[seq(length(variable_2), 1, -i)]

    if(mean(variable) != 0){
        test <- abs(sd(reverse.var) / mean(reverse.var)); test <- ifelse(is.na(test), var.cov, test)
        test_1 <- abs(sd(reverse.var_1) / mean(reverse.var_1)); test_1 <- ifelse(is.na(test_1), var.cov, test_1)
        test_2 <- abs(sd(reverse.var_2) / mean(reverse.var_2)); test_2 <- ifelse(is.na(test_2), var.cov, test_2)
    } else {
        test <- abs(acf(reverse.var, lag.max = 1, plot = FALSE)$acf[2])^-1; test <- ifelse(is.na(test), var.cov, test)
        test_1 <- abs(acf(reverse.var_1, lag.max = 1, plot = FALSE)$acf[2])^-1; test_1 <- ifelse(is.na(test_1), var.cov, test_1)
        test_2 <- abs(acf(reverse.var_2, lag.max = 1, plot = FALSE)$acf[2])^-1; test_2 <- ifelse(is.na(test_2), var.cov, test_2)
    }

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

  ref.output <- cbind(instances, output, output_1, output_2, output * output_1 * output_2 > 0)
  output <- rowMeans(ref.output[ , 2 : 4]) * ref.output[ , 5]

  instances <- ref.output[ , 1] * ref.output[ , 5]

  index <- which(instances > 0 & output > 0)

  insts <- sum(instances > 0) > 0

  if(insts){
    n <- rep(var.cov, length(instances[index]))

    M <- data.table::data.table("Period" = instances[index], "Coefficient.of.Variation" = output[index], "Variable.Coefficient.of.Variation" = n, key = "Coefficient.of.Variation")
  } else {
    M <- data.table::data.table("Period" = 1, "Coefficient.of.Variation" = var.cov, "Variable.Coefficient.of.Variation" = var.cov, key = "Coefficient.of.Variation")
  }





    if(!is.null(modulo)){
        a <- M$Period
        plus <- a+(modulo-a%%modulo)
        minus <- a-a%%modulo

        periods <- unique(c(rbind(minus,plus)))

        if(mod.only){
            periods <- c(periods[!is.na(periods) & periods>0])
            mod_index <- which(unlist(M[, 1])%in%periods)
        } else {
            if(!1%in%unlist(M[,1])){
                periods <- c(periods[!is.na(periods) & periods>0], 1)
            } else {
                periods <- c(periods[!is.na(periods) & periods>0])
            }
            mod_index <- seq_along(unlist(M[,1]))
        }

        periods <- unique(periods[!periods%in%unlist(M[, 1])])

        mod_cv <- data.table::data.table(cbind(periods,
                                   rep(M[1, 3], length(periods)),
                                   rep(M[1, 3], length(periods))))

        M <- data.table::rbindlist(list(M[mod_index, ], mod_cv), use.names = FALSE)
    }

    M <- M[Period < length(variable)/2,]

    if(plot){
        plot(unlist(M[, 1]), unlist(M[, 2]), xlab = "Period", ylab = "Coefficient of Variation", main = "Seasonality Test", ylim = c(0, 2 * abs(sd(variable) / mean(variable))))

        points(unlist(M[, 1])[1], unlist(M[, 2])[1], pch = 19, col = 'red')

        abline(h = abs(sd(variable) / mean(variable)), col = "red", lty = 5)
        text(mean(unlist(M[, 1])), abs(sd(variable) / mean(variable)), pos = 3, "Variable Coefficient of Variation", col = 'red')
    }

    return(list("all.periods" = M,
                "best.period" = unlist(M[1, 1]),
                "periods" = as.vector(unlist(M[, 1]))))

}



