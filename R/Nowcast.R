#' NNS Nowcast
#'
#' Wrapper function for NNS nowcasting method using \link{NNS.VAR} as detailed in Viole (2020), \url{https://www.ssrn.com/abstract=3586658}.
#'
#' @param h integer; \code{(h = 1)} (default) Number of periods to forecast. \code{(h = 0)} will return just the interpolated and extrapolated values.
#' @param additional.regressors character; \code{NULL} (default) add more regressors to the base model.  The format must utilize the \code{\link[quantmod]{getSymbols}} format for FRED data.
#' @param specific.regressors integer; \code{NULL} (default) Select individual regressors from the base model per Viole (2020) listed in the References.
#' @param start.date character; \code{"2000-01-03"} (default) Starting date for all data series download.
#' @param keep.data logical; \code{FALSE} (default) Keeps downloaded variables in a new environment \code{NNSdata}.
#' @param status logical; \code{TRUE} (default) Prints status update message in console.
#' @param ncores integer; value specifying the number of cores to be used in the parallelized subroutine \link{NNS.ARMA.optim}. If NULL (default), the number of cores to be used is equal to the number of cores of the machine - 1.
#'
#' @return Returns the following matrices of forecasted variables:
#' \itemize{
#'  \item{\code{"interpolated_and_extrapolated"}} Returns a \code{data.frame} of the linear interpolated and \link{NNS.ARMA} extrapolated values to replace \code{NA} values in the original \code{variables} argument.  This is required for working with variables containing different frequencies, e.g. where \code{NA} would be reported for intra-quarterly data when indexed with monthly periods.
#'  \item{\code{"relevant_variables"}} Returns the relevant variables from the dimension reduction step.
#'
#'  \item{\code{"univariate"}} Returns the univariate \link{NNS.ARMA} forecasts.
#'
#'  \item{\code{"multivariate"}} Returns the multi-variate \link{NNS.reg} forecasts.
#'
#'  \item{\code{"ensemble"}} Returns the ensemble of both \code{"univariate"} and \code{"multivariate"} forecasts.
#'  }
#'
#'
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{https://www.amazon.com/dp/1490523995/ref=cm_sw_su_dp}
#'
#' Viole, F. (2019) "Multi-variate Time-Series Forecasting: Nonparametric Vector Autoregression Using NNS"
#' \url{https://www.ssrn.com/abstract=3489550}
#'
#' Viole, F. (2020) "NOWCASTING with NNS"
#' \url{https://www.ssrn.com/abstract=3589816}
#'
#'
#' @examples
#'
#'  \dontrun{
#'  NNS.nowcast(h = 12)
#'  }
#'
#' @export


NNS.nowcast <- function(h = 1,
                        additional.regressors = NULL,
                        specific.regressors = NULL,
                        start.date = "2000-01-03",
                        keep.data = FALSE,
                        status = TRUE,
                        ncores = NULL){


  variables <- c("PAYEMS", "JTSJOL",  "CPIAUCSL", "DGORDER", "RSAFS",
                 "UNRATE", "HOUST", "INDPRO", "DSPIC96", "BOPTEXP",
                 "BOPTIMP", "TTLCONS", "IR", "CPILFESL", "PCEPILFE",
                 "PCEPI", "PERMIT", "TCU", "BUSINV", "ULCNFB",
                 "IQ", "GACDISA066MSFRBNY", "GACDFSA066MSFRBPHI", "PCEC96", "GDPC1",
                 "ICSA", "DGS10", "T10Y2Y", "WALCL")


  if(is.null(specific.regressors)) variable_list <- as.character(c(variables, additional.regressors)) else variable_list <- as.character(c(variables, additional.regressors)[specific.regressors])

  NNSdata <- new.env()
  quantmod::getSymbols(variable_list, src = "FRED", env = NNSdata)

  if(sum(ls(envir = NNSdata)%in%variable_list) < length(variable_list)){
    quantmod::getSymbols(variable_list[!variable_list%in%ls(envir = NNSdata)[(ls(envir = NNSdata)%in%variable_list)]], src = "FRED")
  }
 
  raw_econ_variables <- lapply(mget(variable_list, envir = NNSdata), function(x) xts::to.monthly(x)[,4])
  
  oldw <- getOption("warn")
  options(warn = -1)
  
  if(!keep.data) rm(list = ls(), envir = NNSdata)
  
  econ_variables <- Reduce(function(...) merge(..., all=TRUE), raw_econ_variables)[paste0(start.date,"::")]

  colnames(econ_variables) <- variable_list
  
  options(warn = oldw)
  
  NNS.VAR(econ_variables, h = h, tau = 12, status = status, ncores = ncores, nowcast = TRUE)
}
