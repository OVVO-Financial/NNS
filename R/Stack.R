#' NNS stack
#'
#' Prediction model using the predictions of the NNS base models \link{NNS.reg} as features (i.e. meta-features) for the stacked model.
#'
#' @param IVs.train a vector, matrix or data frame of variables of numeric or factor data types.
#' @param IVs.test a vector, matrix or data frame of variables of numeric or factor data types.
#' @param DV.train a numeric or factor vector with compatible dimsensions to \code{(IVs.train)}.
#' @param CV.size numeric [0, 1]; \code{NULL} (default) Sets the cross-validation size if \code{(IVs.test = NULL)}.  Defaults to 0.25 for a 25 percent random sampling of the training set under \code{(CV.size = NULL)}.
#' @param weight options: ("SSE", "Features") method for selecting model output weight; Set \code{(weight = "SSE")} for optimum parameters and weighting based on each base model's sum of squared errors.  \code{(weight = "Feautures")} uses a weighting based on the number of features present, whereby logistic \link{NNS.reg} receives higher relative weights for more regressors.  Defaults to \code{"SSE"}.
#' @param order integer; \code{NULL} (default) Sets the order for \link{NNS.reg}, where \code{(order = 'max')} is the k-nearest neighbors equivalent.
#' @param norm options: ("std", "NNS", NULL); \code{NULL} (default) 3 settings offered: \code{NULL}, \code{"std"}, and \code{"NNS"}.  Selects the \code{norm} parameter in \link{NNS.reg}.
#' @param method numeric options: (1, 2); Select the NNS method to include in stack.  \code{(method = 1)} selects \link{NNS.reg}; \code{(method = 2)} selects \link{NNS.reg} dimension reduction regression.  Defaults to \code{method = c(1, 2)}, including both NNS regression methods in the stack.
#' @param dim.red.method options: ("cor", "NNS.cor", "NNS.caus", "all") method for determining synthetic X* coefficients.  \code{(dim.red.method = "cor")} (default) uses standard linear correlation for weights.  \code{(dim.red.method = "NNS.cor")} uses \link{NNS.cor} for nonlinear correlation weights, while \code{(dim.red.method = "NNS.caus")} uses \link{NNS.caus} for causal weights.  \code{(dim.red.method = "all")} averages all methods for further feature engineering.
#' @param seed numeric; 123 (default) Sets seed for CV sampling.
#' @return Returns a vector of fitted values for the dependent variable test set for all models.
#' \itemize{
#' \item{\code{"NNS.reg.n.best"}} returns the optimum \code{"n.best"} paramater for the \link{NNS.reg} multivariate regression.  \code{"SSE.reg"} returns the SSE for the \link{NNS.reg} multivariate regression.
#' \item{\code{"NNS.dim.red.threshold"}} returns the optimum \code{"threshold"} from the \link{NNS.reg} dimension reduction regression.
#' \item{\code{"SSE.dim.red"}} returns the SSE for the \link{NNS.reg} dimension reduction regression.
#' \item{\code{"reg"}} returns \link{NNS.reg} output.
#' \item{\code{"dim.red"}} returns \link{NNS.reg} dimension reduction regression output.
#' \item{\code{"stack"}} returns the output of the stacked model.
#' }
#'
#' @keywords classifier
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. (2016) "Classification Using NNS Clustering Analysis"
#' \url{https://ssrn.com/abstract=2864711}
#' @note If character variables are used, transform them first to factors using \link{as.factor}, or \link{data.matrix} to ensure overall dataset is numeric.  A multifunction \link{sapply} can also be applied to the overall dataset: \code{data <- sapply(data,function(x){as.factor(x) ; as.numeric(x)})}.  Then run \code{NNS.stack} with transormed variables.
#'
#' Missing data should be handled prior as well using \link{na.omit} or \link{complete.cases} on the full dataset.
#'
#' If error received:
#'
#' \code{"Error in is.data.frame(x) : object 'RP' not found"}
#'
#' reduce the \code{CV.size}.
#'
#'
#' @examples
#'  ## Using 'iris' dataset where test set [IVs.test] is 'iris' rows 141:150.
#'  \dontrun{
#'  NNS.stack(iris[1:140, 1:4], iris[1:140, 5], IVs.test = iris[141:150, 1:4])}
#'
#'  ## Using 'iris' dataset to determine [n.best] and [threshold] with no test set.
#'  \dontrun{
#'  NNS.stack(iris[ , 1:4], iris[ , 5])}
#'
#'  ## Selecting NNS.reg and dimension reduction techniques.
#'  \dontrun{
#'  NNS.stack(iris[1:140, 1:4], iris[1:140, 5], iris[141:150, 1:4], method = c(1, 2))}
#' @export

NNS.stack <- function(IVs.train,
                      DV.train,
                      IVs.test = NULL,
                      CV.size = NULL,
                      weight = "SSE",
                      order = NULL,
                      norm = NULL,
                      method = c(1, 2),
                      dim.red.method = "cor",
                      seed = 123){



  IVs.train <- apply(IVs.train, 2, as.numeric)
  DV.train <- as.numeric(DV.train)

  n <- ncol(IVs.train)

  l = length(IVs.train[ , 1])

  if(is.null(CV.size)){
    if(is.null(IVs.test)){
      CV.size = 0.25
    } else {
        CV.size = mean(c(.2, min(length(IVs.test[ , 1]) / l, .5)))
    }
  }

  #IV test provided...
  if(!is.null(IVs.test)){
    IVs.test <- data.matrix(IVs.test)
  }

  #No IV test provided...

  set.seed(seed * l)
  test.set = sample(1 : l, as.integer(CV.size * l), replace = FALSE)

  CV.IVs.train <- IVs.train[c(-test.set), ]

  CV.IVs.test <- IVs.train[c(test.set), ]
  CV.DV.train <- DV.train[c(-test.set)]
  CV.DV.test <- DV.train[c(test.set)]

  IVs.train <- CV.IVs.train
  DV.train <- CV.DV.train

  ### NORMALIZATION OF VARIABLES and SELECTION OF ORDER:
  np = nrow(CV.IVs.test)
  points.norm = rbind(CV.IVs.test, CV.IVs.train)
  colnames(points.norm) = colnames(CV.IVs.test)
  if(!is.null(norm)){
    if(norm == "std"){
      CV.IVs.train = apply(CV.IVs.train, 2, function(b) (b - min(b)) / (max(b) - min(b)))
      CV.IVs.test = apply(points.norm, 2, function(b) (b - min(b)) / (max(b) - min(b)))[1 : np, ]
    }

    if(norm == "NNS"){
      CV.IVs.train = NNS.norm(CV.IVs.train)
      CV.IVs.test = NNS.norm(points.norm)[1 : np, ]
    }
  }

  if(1 %in% method){

    nns.cv = sapply(1 : (2 * n), function(i) sum((NNS.reg(CV.IVs.train, CV.DV.train, point.est = CV.IVs.test, plot = FALSE, n.best = i, order=order)$Point.est - CV.DV.test) ^ 2))

    nns.cv = c(nns.cv, sum((NNS.reg(CV.IVs.train, CV.DV.train, point.est = CV.IVs.test, plot = FALSE, n.best = 'all', order = order)$Point.est - CV.DV.test) ^ 2))

    k = which.min(na.omit(nns.cv))
    if(k == length(nns.cv)){
      k = 'all'
    }

    nns.1 = NNS.reg(IVs.train, DV.train, point.est = IVs.test, plot = FALSE, n.best = k)$Point.est
    nns.cv.sse = min(nns.cv)

  } else {
    k = 'N/A'
    nns.1 = 0
    nns.cv = 1e-10
    nns.cv.sse = 'N/A'
  }

  # Dimension Reduction Regression Output
  if(2 %in% method){
    var.cutoffs = abs(round(NNS.reg(CV.IVs.train, CV.DV.train, dim.red.method = dim.red.method, plot = FALSE)$equation$Coefficient, digits = 2))

    var.cutoffs = unique(var.cutoffs[var.cutoffs < max(var.cutoffs)])

    nns.ord = sapply(1 : length(var.cutoffs), function(i) sum((NNS.reg(CV.IVs.train, CV.DV.train, point.est = CV.IVs.test, plot = FALSE, dim.red.method = dim.red.method, threshold = var.cutoffs[i])$Point.est - CV.DV.test) ^ 2))

    nns.2 = NNS.reg(IVs.train, DV.train,point.est = IVs.test, dim.red.method = dim.red.method, plot = FALSE, threshold = var.cutoffs[which.min(nns.ord)])$Point.est

    nns.ord.sse = min(na.omit(nns.ord))
  } else {
    nns.2 = 0
    nns.ord = 1e-10
    nns.ord.sse = NA
    var.cutoffs = NA
  }



  ### Weights for combining NNS techniques
  if(weight == "Features"){
    weights = c(n, n - 1)
  }

  nns.cv[nns.cv == 0] <- 1e-10
  nns.ord[nns.ord == 0] <- 1e-10


  if(weight == "SSE"){
    weights = c(max(1e-10, 1 / min(na.omit(nns.cv))), max(1e-10, 1 / min(na.omit(nns.ord))))
  }



  weights = pmax(weights, c(0, 0))
  weights[!(c(1, 2) %in% method)] <- 0
  weights = weights / sum(weights)


  estimates = (weights[1] * nns.1 + weights[2] * nns.2)

  return(list(NNS.reg.n.best = k,
              sse.reg = nns.cv.sse,
              NNS.dim.red.threshold = var.cutoffs[which.min(nns.ord)],
              SSE.dim.red = nns.ord.sse,
              reg = nns.1,
              dim.red = nns.2,
              stack=estimates))

}
