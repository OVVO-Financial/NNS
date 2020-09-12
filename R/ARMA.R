#' NNS ARMA
#'
#' Autoregressive model incorporating nonlinear regressions of component series.
#'
#' @param variable a numeric vector.
#' @param h integer; 1 (default) Number of periods to forecast.
#' @param training.set numeric; \code{NULL} (default) Sets the number of variable observations
#'
#'  \code{(variable[1 : training.set])} to monitor performance of forecast over in-sample range.
#' @param seasonal.factor logical or integer(s); \code{TRUE} (default) Automatically selects the best seasonal lag from the seasonality test.  To use weighted average of all seasonal lags set to \code{(seasonal.factor = FALSE)}.  Otherwise, directly input known frequency integer lag to use, i.e. \code{(seasonal.factor = 12)} for monthly data.  Multiple frequency integers can also be used, i.e. \code{(seasonal.factor = c(12, 24, 36))}
#' @param weights numeric; \code{NULL} (default) sets the weights of the \code{seasonal.factor} vector when specified as integers.  If \code{(weights = NULL)} each \code{seasonal.factor} is weighted on its \link{NNS.seas} result and number of observations it contains.
#' @param best.periods integer; [2] (default) used in conjunction with \code{(seasonal.factor = FALSE)}, uses the \code{best.periods} number of detected seasonal lags instead of \code{ALL} lags when
#'
#' \code{(seasonal.factor = FALSE)}.
#' @param negative.values logical; \code{FALSE} (default) If the variable can be negative, set to
#' \code{(negative.values = TRUE)}.  If there are negative values within the variable, \code{negative.values} will automatically be detected.
#' @param method options: ("lin", "nonlin", "both"); \code{"nonlin"} (default)  To select the regression type of the component series, select \code{(method = "both")} where both linear and nonlinear estimates are generated.  To use a nonlinear regression, set to
#' \code{(method = "nonlin")}; to use a linear regression set to \code{(method = "lin")}.
#' @param dynamic logical; \code{FALSE} (default) To update the seasonal factor with each forecast point, set to \code{(dynamic = TRUE)}.  The default is \code{(dynamic = FALSE)} to retain the original seasonal factor from the inputted variable for all ensuing \code{h}.
#' @param plot logical; \code{TRUE} (default) Returns the plot of all periods exhibiting seasonality and the \code{variable} level reference in upper panel.  Lower panel returns original data and forecast.
#' @param seasonal.plot logical; \code{TRUE} (default) Adds the seasonality plot above the forecast.  Will be set to \code{FALSE} if no seasonality is detected or \code{seasonal.factor} is set to an integer value.
#' @param conf.intervals numeric [0, 1]; \code{NULL} (default) Plots and returns the associated confidence intervals for the final estimate.  Constructed using the maximum entropy bootstrap \link{meboot} on the final estimates.
#' @param ncores integer; value specifying the number of cores to be used in the parallelized  procedure. If NULL (default), the number of cores to be used is equal to half the number of cores of the machine - 1.
#' @return Returns a vector of forecasts of length \code{(h)} if no \code{conf.intervals} specified.  Else, returns a \link{data.table} with the forecasts as well as lower and upper confidence intervals per forecast point.
#' @note
#' For monthly data series, increased accuracy may be realized from forcing seasonal factors to multiples of 12.  For example, if the best periods reported are: \{37, 47, 71, 73\}  use
#' \code{(seasonal.factor = c(36, 48, 72))}.
#'
#' \code{(seasonal.factor = FALSE)} can be a very computationally expensive exercise due to the number of seasonal periods detected.
#'
#' If error encountered when \code{(seasonal.factor = TRUE)}:
#'
#' \code{"NaNs produced Error in seq.default(length(variable)+1, 1, -lag[i]) :
#'  wrong sign in 'by' argument"}
#'
#' use the combination of \code{(seasonal.factor = FALSE, best.periods = 1)}.
#'
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{https://www.amazon.com/dp/1490523995/ref=cm_sw_su_dp}
#'
#' Viole, F. (2019) "Forecasting Using NNS"
#' \url{https://www.ssrn.com/abstract=3382300}
#'
#' @examples
#'
#' ## Nonlinear NNS.ARMA using AirPassengers monthly data and 12 period lag
#' \dontrun{
#' NNS.ARMA(AirPassengers, h = 45, training.set = 100, seasonal.factor = 12, method = "nonlin")
#'
#' ## Linear NNS.ARMA using AirPassengers monthly data and 12, 24, and 36 period lags
#' NNS.ARMA(AirPassengers, h = 45, training.set = 120, seasonal.factor = c(12, 24, 36), method = "lin")
#'
#' ## Nonlinear NNS.ARMA using AirPassengers monthly data and 2 best periods lag
#' NNS.ARMA(AirPassengers, h = 45, training.set = 120, seasonal.factor = FALSE, best.periods = 2)}
#'
#' @export



# Autoregressive Model
NNS.ARMA <- function(variable,
                     h = 1,
                     training.set = NULL,
                     seasonal.factor = TRUE,
                     weights = NULL,
                     best.periods = 2,
                     negative.values = FALSE,
                     method = "nonlin",
                     dynamic = FALSE,
                     plot = TRUE,
                     seasonal.plot = TRUE,
                     conf.intervals = NULL,
                     ncores = NULL){


  if(is.numeric(seasonal.factor) && dynamic){
      stop('Hmmm...Seems you have "seasonal.factor" specified and "dynamic = TRUE".  Nothing dynamic about static seasonal factors!  Please set "dynamic = FALSE" or "seasonal.factor = FALSE"')
  }

  oldw <- getOption("warn")
  options(warn = -1)


  if (is.null(ncores)) {
      num_cores <- as.integer(detectCores() / 2) - 1
  } else {
      num_cores <- ncores
  }

  if(num_cores>1){
      cl <- makeCluster(num_cores)
      registerDoParallel(cl)
  } else { cl <- NULL }

  if(!is.null(best.periods) && !is.numeric(seasonal.factor)){
      seasonal.factor <- FALSE
  }

  label <- deparse(substitute(variable))
  variable <- as.numeric(variable)
  OV <- variable

  if(min(variable)<0) negative.values <- TRUE

  if(!is.null(training.set)){
      variable <- variable[1 : training.set]
      FV <- variable[1 : training.set]
  } else {
      training.set <- length(variable)
      variable <- variable
      FV <- variable
  }

  Estimates <- numeric()


  if(is.numeric(seasonal.factor)){
      seasonal.plot = FALSE
      M <- matrix(seasonal.factor, ncol=1)
      colnames(M) <- "Period"
      lag <- seasonal.factor
      output <- numeric(length(seasonal.factor))
      for(i in 1 : length(seasonal.factor)){
          rev.var <- variable[seq(length(variable), 1, -i)]
          output[i] <- abs(sd(rev.var) / mean(rev.var))
      }

      if(is.null(weights)){
          Relative.seasonal <- output / abs(sd(variable)/mean(variable))
          Seasonal.weighting <- 1 / Relative.seasonal
          Observation.weighting <- 1 / sqrt(seasonal.factor)
          Weights <- (Seasonal.weighting * Observation.weighting) / sum(Observation.weighting * Seasonal.weighting)
          seasonal.plot <- FALSE
      } else {
          Weights <- weights
      }

  } else {
    M <- NNS.seas(variable, plot=FALSE)
    if(!is.list(M)){
        M <- t(1)
    } else {
        if(is.null(best.periods)){
            M <- M$all.periods
        } else {
            if(!seasonal.factor && is.numeric(best.periods) && (length(M$all.periods$Period) < best.periods)){
                best.periods <- length(M$all.periods$Period)
            }
            if(!seasonal.factor && is.null(best.periods)){
                best.periods <- length(M$all.periods$Period)
            }
        M <- M$all.periods[1 : best.periods, ]
      }
    }


    ASW <- ARMA.seas.weighting(seasonal.factor, M)
    lag <- ASW$lag

    if(is.null(weights)){
        Weights <- ASW$Weights
    } else {
        Weights <- weights
    }
  }



  # Regression for each estimate in h
  for (j in 1 : h){
      ## Regenerate seasonal.factor if dynamic
      if(dynamic){
          seas.matrix <- NNS.seas(variable, plot = FALSE)
          if(!is.list(seas.matrix)){
              M <- t(1)
          } else {
              if(is.null(best.periods)){
                  M <- seas.matrix$all.periods
                  best.periods <- length(M$all.periods$Period)
              } else {
                  if(length(M$all.periods$Period) < best.periods){
                      best.periods <- length(M$all.periods$Period)
                  }
                  M <- seas.matrix$all.periods[1 : best.periods, ]
              }
          }

      ASW <- ARMA.seas.weighting(seasonal.factor, M)
      lag <- ASW$lag
      Weights <- ASW$Weights
    }

    ## Generate vectors for 1:lag
    GV <- generate.vectors(variable,lag)
    Component.index <- GV$Component.index
    Component.series <- GV$Component.series

    ## Regression on Component Series
    #Regression.Estimates <- numeric()

    if(method == 'nonlin' | method == 'both'){
      Regression.Estimates <- list(length(lag))

      Regression.Estimates <- foreach(i = 1 : length(lag),.packages = "NNS")%dopar%{
        x <- Component.index[[i]] ; y <- Component.series[[i]]
        last.y <- tail(y, 1)

        ## Skeleton NNS regression for NNS.ARMA
        reg.points <- tail(NNS.reg(x, y, return.values = FALSE , plot = FALSE, multivariate.call = TRUE), 4)
        reg.points <- reg.points[complete.cases(reg.points),]

        run <- mean(diff(reg.points$x))
        rise <- mean(diff(reg.points$y))

        last.y + (rise / run)
      }

      Regression.Estimates <- unlist(Regression.Estimates)

      NL.Regression.Estimates <- Regression.Estimates
      Nonlin.estimates <- sum(Regression.Estimates * Weights)

    }#Linear == F

    if(method == 'lin' | method == 'both'){

      Regression.Estimates <- list(length(lag))

      Regression.Estimates <- foreach(i = 1 : length(lag))%dopar%{
          last.x <- tail(Component.index[[i]], 1)
          coefs <- coef(lm(Component.series[[i]] ~ Component.index[[i]]))

          coefs[1] + (coefs[2] * (last.x + 1))
      }

      Regression.Estimates <- unlist(Regression.Estimates)

      L.Regression.Estimates <- Regression.Estimates
      Lin.estimates <- sum(Regression.Estimates * Weights)

    }#Linear == T


    if(!negative.values){
        Regression.Estimates <- pmax(0, Regression.Estimates)
    }



    if(method == 'both'){
        Estimates[j] <- mean(c(Lin.estimates, Nonlin.estimates))
    } else {
        Estimates[j] <- sum(Regression.Estimates * Weights)
    }

    variable <- c(variable, Estimates[j])
    FV <- variable

  } # j loop

if(!is.null(cl)){
    stopCluster(cl)
    registerDoSEQ()
}


  if(!is.null(conf.intervals)){
      CIs <- NNS.meboot(Estimates, reps=399)$replicates
  }


  #### PLOTTING
  if(plot){
    original.par = par(no.readonly = TRUE)
    if(seasonal.plot){
      par(mfrow = c(2, 1))
      if(ncol(M) > 1){
        plot(M[, Period], M[, Coefficient.of.Variation],
             xlab = "Period", ylab = "Coefficient of Variation", main = "Seasonality Test", ylim = c(0, 2 * M[1, Variable.Coefficient.of.Variation]))
        points(M[ , Period], M[ , Coefficient.of.Variation], pch = 19, col = 'red')
        abline(h = M[1, Variable.Coefficient.of.Variation], col = "red", lty = 5)
        text((M[ , min(Period)] + M[ , max(Period)]) / 2, M[1, Variable.Coefficient.of.Variation], pos = 3, "Variable Coefficient of Variation", col = 'red')
      } else {
        plot(1,1, pch = 19, col = 'blue', xlab = "Period", ylab = "Coefficient of Variation", main = "Seasonality Test",
             ylim = c(0, 2 * abs(sd(FV) / mean(FV))))
        text(1, abs(sd(FV) / mean(FV)), pos = 3, "NO SEASONALITY DETECTED", col = 'red')
      }
    }


    if(is.null(label)){
      label <- "Variable"
    }


    if(!is.null(conf.intervals)){
      plot(OV, type = 'l', lwd = 2, main = "NNS.ARMA Forecast", col = 'steelblue',
           xlim = c(1, max((training.set + h), length(OV))),
           ylab = label, ylim = c(min(Estimates, OV,  unlist(CIs) ), max(OV, Estimates, unlist(CIs) )) )

      for(i in 1 : 399){
        lines((training.set+1) : (training.set+h), CIs[,i],  col = rgb(0.75,0.75,0.75, 0.05))
      }

      lines((training.set + 1) : (training.set + h), Estimates, type = 'l', lwd = 2, lty = 1, col = 'red')
      segments(training.set, FV[training.set], training.set + 1, Estimates[1],lwd = 2,lty = 1,col = 'red')
      legend('topleft', bty = 'n', legend = c("Original", paste0("Forecast ", h, " period(s)")), lty = c(1, 1), col = c('steelblue', 'red'), lwd = 2)
    } else {
      plot(OV, type = 'l', lwd = 2, main = "NNS.ARMA Forecast", col = 'steelblue',
           xlim = c(1, max((training.set + h), length(OV))),
           ylab = label, ylim = c(min(Estimates, OV), max(OV, Estimates)))

      if(training.set[1] < length(OV)){
        lines((training.set + 1) : (training.set + h), Estimates, type = 'l',lwd = 2, lty = 3, col = 'red')
        segments(training.set, FV[training.set], training.set + 1, Estimates[1], lwd = 2, lty = 3, col = 'red')
        legend('topleft', bty = 'n', legend = c("Original", paste0("Forecast ", h, " period(s)")), lty = c(1, 2), col = c('steelblue', 'red'), lwd = 2)
      } else {
        lines((training.set + 1) : (training.set + h), Estimates, type = 'l', lwd = 2, lty = 1, col = 'red')
        segments(training.set, FV[training.set], training.set + 1, Estimates[1], lwd = 2, lty = 1, col = 'red')
        legend('topleft', bty = 'n', legend = c("Original", paste0("Forecast ", h, " period(s)")),lty = c(1, 1), col = c('steelblue', 'red'), lwd = 2)
      }


    }
    points(training.set, OV[training.set], col = "green", pch = 18)
    points(training.set + h, tail(FV, 1), col = "green", pch = 18)

    par(original.par)
  }

  options(warn = oldw)
  if(!is.null(conf.intervals)){
      upper_CIs <- apply(CIs, 1, function(z) UPM.VaR(1-conf.intervals, 0, z))
      lower_CIs <- apply(CIs, 1, function(z) LPM.VaR(1-conf.intervals, 0, z))
      results <- cbind.data.frame(Estimates,  lower_CIs,  upper_CIs)
      colnames(results) = c("Estimates",
                            paste0("Lower ", round(conf.intervals*100,2), "% CI"),
                            paste0("Upper ", round(conf.intervals*100,2), "% CI"))
      return(data.table::data.table(results))
  } else {
      return(Estimates)
  }
}
