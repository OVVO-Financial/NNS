#' NNS Nowcast
#'
#' Wrapper function for NNS nowcasting method using \link{NNS.VAR} as detailed in Viole (2020), \url{https://www.ssrn.com/abstract=3586658}.
#'
#' @param h integer; 1 (default) Number of periods to forecast. \code{(h = 0)} will return just the interpolated and extrapolated values.
#' @param additional.regressors character; \code{NULL} (default) add more regressors to the base model.  The format must utilize the Quandl exchange format as described in \url{https://docs.data.nasdaq.com/docs/data-organization}.  For example, the 10-year US Treasury yield using the St. Louis Federal Reserve data is \code{"FRED/DGS10"}.
#' @param start.date character; \code{"2000-01-03"} (default) Starting date for all data series download.
#' @param Quandl.key character; \code{NULL} (default) User provided \link{Quandl} API key WITH QUOTES.
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
#' @note
#' \itemize{
#' \item This function requires an API key from Quandl.  Sign up via \url{https://data.nasdaq.com/}.
#' }
#'
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{https://www.amazon.com/dp/1490523995/ref=cm_sw_su_dp}
#'
#' Viole, F. (2019) "Multi-variate Time-Series Forecasting: Nonparametric Vector Autoregression Using NNS"
#' \url{https://www.ssrn.com/abstract=3489550}
#'
#' Viole, F. (2020) "NOWCASTING with NNS"
#' \url{https://www.ssrn.com/abstract=3586658}
#'
#'
#' @examples
#'
#'  \dontrun{
#'  NNS.nowcast(h = 12)
#'  }
#'
#' @export


NNS.nowcast <- function(h = 12,
                        additional.regressors = NULL,
                        start.date = "2000-01-03",
                        Quandl.key = NULL,
                        status = TRUE,
                        ncores = NULL){


  if(is.null(Quandl.key)){
    message("Please enter your Quandl API key WITHOUT QUOTES and press enter:")

    key <- readline(": ")
    Quandl::Quandl.api_key(key)

  } else { Quandl::Quandl.api_key(Quandl.key) }

  variables <- c("PAYEMS", "JTSJOL",  "CPIAUCSL", "DGORDER", "RSAFS",
                 "UNRATE", "HOUST", "INDPRO", "DSPIC96", "BOPTEXP",
                 "BOPTIMP", "TTLCONS", "IR", "CPILFESL", "PCEPILFE",
                 "PCEPI", "PERMIT", "TCU", "BUSINV", "ULCNFB",
                 "IQ", "GACDISA066MSFRBNY", "GACDFSA066MSFRBPHI", "PCEC96", "GDPC1",
                 "DGS10", "T10Y2Y", "ICSA", "DGS10", "T10Y2Y", "ICSA")

  variable_list <- character()

  for(i in 1:length(variables)){
    variable_list[i] <- paste0("FRED/", variables[i])
  }

  variable_list <- c(variable_list, additional.regressors)

  econ_variables <- Quandl::Quandl(variable_list, type = 'ts',order = "asc",
                                   collapse = "monthly", start_date = start.date)

  if(!is.null(attributes(econ_variables)$errors)) return(attributes(econ_variables)$errors)

  econ_variables <- stats::.preformat.ts(econ_variables)

  colnames(econ_variables) <- variable_list

  NNS.VAR(econ_variables, h = h, tau = 12, status = status, ncores = ncores, nowcast = TRUE)
}
