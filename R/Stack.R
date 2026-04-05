#' NNS Stack
#'
#' Prediction model using the predictions of the NNS base models \link{NNS.reg} as features (i.e. meta-features) for the stacked model.
#'
#' @param IVs.train a vector, matrix or data frame of variables of numeric or factor data types.
#' @param DV.train a numeric or factor vector with compatible dimensions to \code{(IVs.train)}.
#' @param IVs.test a vector, matrix or data frame of variables of numeric or factor data types with compatible dimensions to \code{(IVs.train)}.  If NULL, will use \code{(IVs.train)} as default.
#' @param type \code{NULL} (default).  To perform a classification of discrete integer classes from factor target variable \code{(DV.train)} with a base category of 1, set to \code{(type = "CLASS")}, else for continuous \code{(DV.train)} set to \code{(type = NULL)}.   Like a logistic regression, this setting is not necessary for target variable of two classes e.g. [0, 1].
#' @param obj.fn expression; \code{expression(sum((predicted - actual)^2))} (default) Sum of squared errors is the default objective function.  Any \code{expression()} using the specific terms \code{predicted} and \code{actual} can be used.
#' @param objective options: ("min", "max") \code{"min"} (default) Select whether to minimize or maximize the objective function \code{obj.fn}.
#' @param optimize.threshold logical; \code{TRUE} (default) Will optimize the probability threshold value for rounding in classification problems.  If \code{FALSE}, returns 0.5.
#' @param dist options:("L1", "L2", "DTW", "FACTOR") the method of distance calculation; Selects the distance calculation used. \code{dist = "L2"} (default) selects the Euclidean distance and \code{(dist = "L1")} selects the Manhattan distance; \code{(dist = "DTW")} selects the dynamic time warping distance; \code{(dist = "FACTOR")} uses a frequency.
#' @param CV.size numeric [0, 1]; \code{NULL} (default) Sets the cross-validation size if \code{(IVs.test = NULL)}.  Defaults to a random value between 0.2 and 0.33 for a random sampling of the training set.
#' @param balance logical; \code{FALSE} (default) Uses both up and down sampling to balance the classes.  \code{type="CLASS"} required.
#' @param ts.test integer; NULL (default) Sets the length of the test set for time-series data; typically \code{2*h} parameter value from \link{NNS.ARMA} or double known periods to forecast.
#' @param folds integer; \code{folds = 5} (default) Select the number of cross-validation folds.
#' @param order options: (integer, "max", NULL); \code{NULL} (default) Sets the order for \link{NNS.reg}, where \code{(order = "max")} is the k-nearest neighbors equivalent, which is suggested for mixed continuous and discrete (unordered, ordered) data.
#' @param method numeric options: (1, 2); Select the NNS method to include in stack.  \code{(method = 1)} selects \link{NNS.reg}; \code{(method = 2)} selects \link{NNS.reg} dimension reduction regression.  Defaults to \code{method = c(1, 2)}, which will reduce the dimension first, then find the optimal \code{n.best}.
#' @param stack logical; \code{TRUE} (default) Uses dimension reduction output in \code{n.best} optimization, otherwise performs both analyses independently.
#' @param dim.red.method options: ("cor", "NNS.dep", "NNS.caus", "equal", "all") method for determining synthetic X* coefficients.  \code{(dim.red.method = "cor")} uses standard linear correlation for weights.  \code{(dim.red.method = "NNS.dep")} (default) uses \link{NNS.dep} for nonlinear dependence weights, while \code{(dim.red.method = "NNS.caus")} uses \link{NNS.caus} for causal weights.  \code{(dim.red.method = "all")} averages all methods for further feature engineering.
#' @param pred.int numeric [0,1]; \code{NULL} (default) Returns the associated prediction intervals with each \code{method}.
#' @param status logical; \code{TRUE} (default) Prints status update message in console.
#' @param ncores integer; value specifying the number of cores to be used in the parallelized subroutine \link{NNS.reg}. If NULL (default), the number of cores to be used is equal to the number of cores of the machine - 1.
#'
#' @return Returns a vector of fitted values for the dependent variable test set for all models.
#' \itemize{
#' \item{\code{"NNS.reg.n.best"}} returns the optimum \code{"n.best"} parameter for the \link{NNS.reg} multivariate regression.  \code{"SSE.reg"} returns the SSE for the \link{NNS.reg} multivariate regression.
#' \item{\code{"OBJfn.reg"}} returns the \code{obj.fn} for the \link{NNS.reg} regression.
#' \item{\code{"NNS.dim.red.threshold"}} returns the optimum \code{"threshold"} from the \link{NNS.reg} dimension reduction regression.
#' \item{\code{"OBJfn.dim.red"}} returns the \code{obj.fn} for the \link{NNS.reg} dimension reduction regression.
#' \item{\code{"probability.threshold"}} returns the optimum probability threshold for classification, else 0.5 when set to \code{FALSE}.
#' \item{\code{"reg"}} returns \link{NNS.reg} output.
#' \item{\code{"reg.pred.int"}} returns the prediction intervals for the regression output.
#' \item{\code{"dim.red"}} returns \link{NNS.reg} dimension reduction regression output.
#' \item{\code{"dim.red.pred.int"}} returns the prediction intervals for the dimension reduction regression output.
#' \item{\code{"stack"}} returns the output of the stacked model.
#' \item{\code{"pred.int"}} returns the prediction intervals for the stacked model.
#' }
#'
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. (2016) "Classification Using NNS Clustering Analysis"  \doi{10.2139/ssrn.2864711}
#'
#' @note
#' \itemize{
#' \item Incorporate any objective function from external packages (such as \code{Metrics::mape}) via \code{NNS.stack(..., obj.fn = expression(Metrics::mape(actual, predicted)), objective = "min")}
#' 
#' \item Like a logistic regression, the \code{(type = "CLASS")} setting is not necessary for target variable of two classes e.g. [0, 1].  The response variable base category should be 1 for multiple class problems.
#'
#' \item Missing data should be handled prior as well using \link{na.omit} or \link{complete.cases} on the full dataset.
#' }
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
#'  NNS.stack(iris[1:140, 1:4], iris[1:140, 5], IVs.test = iris[141:150, 1:4], type = "CLASS", 
#'  balance = TRUE)
#'
#'  ## Using 'iris' dataset to determine [n.best] and [threshold] with no test set.
#'  NNS.stack(iris[ , 1:4], iris[ , 5], type = "CLASS")
#'  }
#' @export

NNS.stack <- function(IVs.train,
                      DV.train,
                      IVs.test = NULL,
                      type = NULL,
                      obj.fn = expression( sum((predicted - actual)^2) ),
                      objective = "min",
                      optimize.threshold = TRUE,
                      dist = "L2",
                      CV.size = NULL,
                      balance = FALSE,
                      ts.test = NULL,
                      folds = 5,
                      order = NULL,
                      method = c(1, 2),
                      stack = TRUE,
                      dim.red.method = "cor",
                      pred.int = NULL,
                      status = TRUE,
                      ncores = NULL){
  
  if(anyNA(cbind(IVs.train,DV.train))) stop("You have some missing values, please address.")
  if(is.null(obj.fn)) stop("Please provide an objective function")
  
  if(balance && is.null(type)) warning("type = 'CLASS' selected due to balance = TRUE.")
  if(balance) type <- "CLASS"
  
  if(!is.null(type) && min(as.numeric(DV.train))==0) warning("Base response variable category should be 1, not 0.")
  
  if(any(class(IVs.train)%in%c("tbl","data.table"))) IVs.train <- as.data.frame(IVs.train)
  if(any(class(DV.train)%in%c("tbl","data.table"))) DV.train <- as.vector(unlist(DV.train))
  
  if(is.vector(IVs.train) || is.null(dim(IVs.train)) || ncol(IVs.train)==1){
    IVs.train <- data.frame(IVs.train)
    method <- 1
    order <- NULL
  }
  
  if(!is.null(type)){
    type <- tolower(type)
    if(type == "class" && identical(obj.fn,expression( sum((predicted - actual)^2) ))){
      obj.fn <- expression(mean( predicted == as.numeric(actual)))
      objective <- "max"
    }
  }
  
  objective <- tolower(objective)
  
  if(!is.null(type) && type=="class"){
    DV.train <- as.numeric(factor(DV.train))
    smoothness <- FALSE
  } else {
    smoothness <- FALSE
    DV.train <- as.numeric(DV.train)
  }
  
  n <- ncol(IVs.train)
  l <- floor(sqrt(length(IVs.train[ , 1])))
  
  if(is.null(IVs.test)){
    IVs.test <- IVs.train
  } else {
    if(any(class(IVs.test)%in%c("tbl","data.table"))) IVs.test <- as.data.frame(IVs.test)
  }
  
  if(is.null(dim(IVs.test))) IVs.test <- data.frame(t(IVs.test)) else IVs.test <- data.frame(IVs.test)
  
  dist <- tolower(dist)
  
  i_s <- numeric()
  THRESHOLDS <- vector(mode = "list", folds)
  best.k <- vector(mode = "list", folds)
  best.nns.cv <- vector(mode = "list", folds)
  best.nns.ord <- vector(mode = "list", folds)
  
  if(is.null(colnames(IVs.train))){
    colnames.list <- lapply(1 : ncol(IVs.train), function(i) paste0("X", i))
    colnames(IVs.test) <- colnames(IVs.train) <- as.character(colnames.list)
  } else {
    # FIX: if names exist on training and dimensions match, mirror them on the test matrix
    if(!is.null(IVs.test) && ncol(IVs.test) == ncol(IVs.train))
      colnames(IVs.test) <- colnames(IVs.train)
  }
  
  # var.cutoffs_1 (full-data importance scores) removed: computing importance on the
  # entire training set before any CV split leaks held-out information into every
  # fold's threshold grid.  The grid is now built solely from fold-local scores
  # (var.cutoffs_2) computed inside the fold loop on CV.IVs.train only.
  
  # Balance applied ONCE before any fold splitting -- fixes the bug where
  # IVs.train/DV.train were overwritten on each fold iteration, causing the same
  # balanced dataset to be reused across folds and leaking information between them.
  if (balance) {
    y_train  <- as.factor(DV.train)
    ycol <- "Class"
    training_1 <- downSample(IVs.train, y_train, list = FALSE, yname = ycol)
    training_2 <- upSample(IVs.train,   y_train, list = FALSE, yname = ycol)
    training_bal <- rbind.data.frame(training_1, training_2)
    IVs.train <- training_bal[, setdiff(names(training_bal), ycol), drop = FALSE]
    DV.train  <- as.numeric(as.factor(training_bal[[ycol]]))
    colnames(IVs.test) <- colnames(IVs.train)
  }
  
  if(is.null(CV.size)) new.CV.size <- round(runif(1, .2, 1/3), 3) else new.CV.size <- CV.size
  
  for(b in 1 : folds){
    if(status) message("Folds Remaining = " , folds-b," ","\r",appendLF=TRUE)
    
    set.seed(123 * b)
    
    test.set <- as.integer(seq(b, length(unlist(IVs.train[ , 1])), length.out = as.integer(new.CV.size * length(unlist(IVs.train[ , 1])))))
    
    if(!is.null(ts.test)){
      test.set <- 1:(length(DV.train) - ts.test)
    }
    
    test.set <- unlist(test.set)
    
    CV.IVs.train <- data.frame(IVs.train[c(-test.set), ])
    
    if(dim(CV.IVs.train)[2]!=dim(IVs.train)[2]) CV.IVs.train <- t(CV.IVs.train)
    if(dim(CV.IVs.train)[2]!=dim(IVs.train)[2]) CV.IVs.train <- t(CV.IVs.train)
    
    CV.IVs.test <- data.frame(IVs.train[test.set, ])
    if(dim(CV.IVs.test)[2]!=dim(IVs.train)[2]) CV.IVs.test <- t(CV.IVs.test)
    if(dim(CV.IVs.test)[2]!=dim(IVs.train)[2]) CV.IVs.test <- t(CV.IVs.test)
    
    CV.DV.train <- DV.train[c(-test.set)]
    CV.DV.test <- DV.train[c(test.set)]
    
    training <- cbind(IVs.train[c(-test.set),], DV.train[c(-test.set)])
    training <- training[complete.cases(training),]
    
    CV.IVs.train <- data.frame(training[, -(ncol(training))])
    CV.DV.train <- as.numeric(training[, ncol(training)])
    
    
    # Dimension Reduction Regression Output
    if (2 %in% method && ncol(IVs.train) > 1) {
      actual <- CV.DV.test
      
      # --- compute per-variable scores for threshold grid ---
      if (dim.red.method == "cor") {
        var.cutoffs_2 <- abs(round(suppressWarnings(
          cor(data.matrix(cbind(CV.DV.train, CV.IVs.train)), method = "spearman")
        )[-1, 1], digits = 2))
      } else {
        var.cutoffs_2 <- abs(round(suppressWarnings(
          NNS.reg(CV.IVs.train, CV.DV.train,
                  dim.red.method = dim.red.method,
                  plot = FALSE, residual.plot = FALSE,
                  order = order, ncores = ncores,
                  type = type, point.only = TRUE, smooth = smoothness)$equation$Coefficient[-(n + 1)]
        ), digits = 2))
      }
      
      # Threshold grid built solely from fold-local scores -- no full-data leakage.
      var.cutoffs <- var.cutoffs_2
      var.cutoffs <- var.cutoffs[var.cutoffs < 1 & var.cutoffs >= 0]
      var.cutoffs[is.na(var.cutoffs)] <- 0
      var.cutoffs <- rev(sort(unique(var.cutoffs)))[-1]
      if (length(var.cutoffs) == 0 || is.null(var.cutoffs)) var.cutoffs <- 0
      if (n == 2) var.cutoffs <- unique(c(var.cutoffs, 0))
      if (dist == "factor" && length(var.cutoffs) > 1) var.cutoffs <- var.cutoffs[-1]
      if (dim.red.method == "equal") var.cutoffs <- 0
      
      # --- evaluate ALL thresholds (no early stopping) ---
      threshold_results_2 <- vector(mode = "list", length = length(var.cutoffs))
      nns.ord <- rep(NA_real_, length(var.cutoffs))
      
      for (i in seq_along(var.cutoffs)) {
        
        predicted <- suppressWarnings(
          NNS.reg(CV.IVs.train, CV.DV.train,
                  point.est = CV.IVs.test,
                  plot = FALSE,
                  dim.red.method = dim.red.method,
                  threshold = var.cutoffs[i],
                  order = order, ncores = ncores,
                  type = NULL, dist = dist,
                  point.only = TRUE, smooth = smoothness)$Point.est
        )
        
        # fill NA predictions with gravity of non-NA (original behavior)
        predicted[is.na(predicted)] <- gravity(na.omit(predicted))
        
        # per-threshold classification rounding (if needed)
        if (!is.null(type)) {
          if (length(unique(predicted)) == 1) {
            pred_matrix <- matrix(replicate(100, predicted), nrow = length(predicted))
          } else {
            pred_matrix <- sapply(seq(.01, .99, .01),
                                  function(z) ifelse(predicted %% 1 < z,
                                                     as.integer(floor(predicted)),
                                                     as.integer(ceiling(predicted))))
          }
          z <- apply(pred_matrix, 2, function(z) mean(z == as.numeric(actual)))
          threshold_results_2[[i]] <- seq(.01, .99, .01)[as.integer(median(which(z == max(z))))]
          predicted <- ifelse(predicted %% 1 < threshold_results_2[[i]],
                              floor(predicted), ceiling(predicted))
          end_if <- TRUE
        } # end if classification
        
        # objective at this threshold
        nns.ord[i] <- eval(obj.fn)
        
        # print threshold + eval(obj.fn) when status = TRUE
        if (status) {
          message(sprintf(
            "Current NNS.reg(... , threshold = %.4f ) | eval(obj.fn) = %.6f | MAX Iterations Remaining = %d",
            var.cutoffs[i],
            nns.ord[i],
            length(var.cutoffs) - i
          ))
        }
      } # end for each threshold
      
      # --- pick best threshold across ALL tested ---
      if (objective == "min") {
        best.idx <- which.min(na.omit(nns.ord))
        best.nns.ord[[b]] <- min(na.omit(nns.ord))
      } else {
        best.idx <- which.max(na.omit(nns.ord))
        best.nns.ord[[b]] <- max(na.omit(nns.ord))
      }
      if (length(best.idx) == 0) best.idx <- 1L  # fallback if all NA
      best.threshold <- var.cutoffs[best.idx]
      THRESHOLDS[[b]] <- best.threshold
      
      # --- downstream: finalize relevant vars and fit once using the chosen threshold ---
      relevant_vars <- colnames(IVs.train)
      if (is.null(relevant_vars)) relevant_vars <- 1:n
      
      # Compute per-fold X* (synthetic dim-red predictor) for use in Method 1 CV
      # when stack = TRUE and both methods are requested.
      # Replicates NNS.reg internals exactly (Regression.R lines 379-419):
      #   norm.x    <- apply(original.variable, 2, function(b) (b-min(b))/(max(b)-min(b)))
      #   X*_train  <- rowSums(eachrow(norm.x, coef, "*")) / sum(abs(coef)>0)
      # For point.est NNS.reg does: points.norm <- apply(rbind(point.est, x), 2, rescale)
      # then new.point.est <- points.norm[test_rows,] %*% coef / xn
      if (stack && identical(sort(method), c(1, 2))) {
        xstar_cv_fit <- suppressWarnings(
          NNS.reg(CV.IVs.train, CV.DV.train,
                  point.est = CV.IVs.test,
                  dim.red.method = dim.red.method,
                  plot = FALSE, residual.plot = FALSE,
                  order = order, threshold = best.threshold,
                  ncores = ncores, type = NULL,
                  dist = dist, point.only = FALSE, smooth = smoothness)
        )
        # X* for training rows -- returned directly
        xstar_CV_train <- as.numeric(unlist(xstar_cv_fit$x.star))
        
        # X* for test rows -- replicate NNS.reg point.est projection exactly:
        # coefficients from equation (all rows except last DENOMINATOR row)
        eq        <- xstar_cv_fit$equation
        coef_vals <- as.numeric(eq$Coefficient[-nrow(eq)])
        xn        <- sum(abs(coef_vals) > 0)
        if (xn == 0) xn <- 1L
        train_mat <- data.matrix(CV.IVs.train)
        test_mat  <- data.matrix(CV.IVs.test)
        np_cv     <- nrow(test_mat)
        # joint normalisation: rbind(test, train) then rescale each column
        joint     <- rbind(test_mat, train_mat)
        joint_norm <- apply(joint, 2, function(b) {
          rng <- max(b) - min(b)
          (b - min(b)) / ifelse(rng == 0, 1, rng)
        })
        test_norm <- joint_norm[seq_len(np_cv), , drop = FALSE]
        # guard: coef_vals must match ncol of test_norm
        if (length(coef_vals) == ncol(test_norm)) {
          xstar_CV_test <- as.numeric(test_norm %*% coef_vals / xn)
        } else {
          # fallback: use training X* mean
          xstar_CV_test <- rep(mean(xstar_CV_train, na.rm = TRUE), np_cv)
        }
        xstar_CV_test[is.na(xstar_CV_test)] <- gravity(na.omit(xstar_CV_test))
      }
      
      if (b == folds) {
        threshold.table <- sort(table(unlist(THRESHOLDS)), decreasing = TRUE)
        nns.ord.threshold <- gravity(as.numeric(names(threshold.table[threshold.table == max(threshold.table)])))
        if (is.na(nns.ord.threshold)) nns.ord.threshold <- 0
        
        nns.method.2 <- NNS.reg(IVs.train, DV.train,
                                point.est = IVs.test,
                                dim.red.method = dim.red.method,
                                plot = FALSE,
                                order = order, threshold = nns.ord.threshold,
                                ncores = ncores,
                                type = type, point.only = FALSE,
                                confidence.interval = pred.int,
                                smooth = smoothness)
        
        actual    <- nns.method.2$Fitted.xy$y
        predicted <- nns.method.2$Fitted.xy$y.hat
        pred.int.2 <- nns.method.2$pred.int
        best.nns.ord <- eval(obj.fn)
        
        # Capture full-data X* for Method 1's final fit (when stacking)
        # Use $x.star for training rows; replicate NNS.reg's joint normalisation for test rows
        if (stack && identical(sort(method), c(1, 2))) {
          xstar_full_train <- as.numeric(unlist(nns.method.2$x.star))
          
          eq        <- nns.method.2$equation
          coef_vals <- as.numeric(eq$Coefficient[-nrow(eq)])
          xn        <- sum(abs(coef_vals) > 0)
          if (xn == 0) xn <- 1L
          train_mat <- data.matrix(IVs.train)
          test_mat  <- data.matrix(IVs.test)
          np_full   <- nrow(test_mat)
          joint     <- rbind(test_mat, train_mat)
          joint_norm <- apply(joint, 2, function(b) {
            rng <- max(b) - min(b)
            (b - min(b)) / ifelse(rng == 0, 1, rng)
          })
          test_norm <- joint_norm[seq_len(np_full), , drop = FALSE]
          if (length(coef_vals) == ncol(test_norm)) {
            xstar_full_test <- as.numeric(test_norm %*% coef_vals / xn)
          } else {
            xstar_full_test <- rep(mean(xstar_full_train, na.rm = TRUE), np_full)
          }
          xstar_full_test[is.na(xstar_full_test)] <- gravity(na.omit(xstar_full_test))
        }
        
        rel_vars <- nns.method.2$equation
        rel_vars <- which(rel_vars$Coefficient > 0)
        rel_vars <- rel_vars[rel_vars <= n]
        if (length(rel_vars) == 0 || is.null(rel_vars)) rel_vars <- 1:n
        
        if (!stack) relevant_vars <- 1:n else relevant_vars <- rel_vars
        if (all(relevant_vars == "FALSE")) relevant_vars <- 1:n
        
        if (!is.null(type) && !is.null(nns.method.2$Point.est)) {
          threshold_results_2 <- mean(unlist(threshold_results_2))
          nns.method.2 <- ifelse(nns.method.2$Point.est %% 1 < threshold_results_2,
                                 floor(nns.method.2$Point.est), ceiling(nns.method.2$Point.est))
          nns.method.2 <- pmin(nns.method.2, max(as.numeric(DV.train)))
          nns.method.2 <- pmax(nns.method.2, min(as.numeric(DV.train)))
        } else {
          nns.method.2 <- nns.method.2$Point.est
        }
      }
      
    } else {
      THRESHOLDS <- NA
      test.set.2 <- NULL
      nns.method.2 <- NA
      if (objective == "min") { best.nns.ord <- Inf } else { best.nns.ord <- -Inf }
      nns.ord.threshold <- NA
      threshold_results_2 <- NA
      relevant_vars <- 1:n
    } # 2 %in% method
    
    
    
    # --- Method 1 (NNS.reg / k-NN path) — optimized: only 1..l plus q, safe C++ calls ---
    if (1 %in% method) {
      actual <- CV.DV.test
      
      # When stacking with Method 2, replace the CV design matrices with cbind(X*, X*)
      # so that n.best is cross-validated over the synthetic dim-red predictor.
      if (stack && identical(sort(method), c(1, 2)) &&
          exists("xstar_CV_train") && !anyNA(xstar_CV_train)) {
        CV.IVs.train <- data.frame(Xstar = xstar_CV_train,
                                   Xstar2 = xstar_CV_train)
        CV.IVs.test  <- data.frame(Xstar = xstar_CV_test,
                                   Xstar2 = xstar_CV_test)
      } else {
        if (is.character(relevant_vars)) relevant_vars <- relevant_vars != ""
        if (is.logical(relevant_vars)) {
          CV.IVs.train <- data.frame(CV.IVs.train[, relevant_vars, drop = FALSE])
          CV.IVs.test  <- data.frame(CV.IVs.test[,  relevant_vars, drop = FALSE])
        }
        if (ncol(CV.IVs.train) != n) CV.IVs.train <- t(CV.IVs.train)
        if (ncol(CV.IVs.train) != n) CV.IVs.train <- t(CV.IVs.train)
        if (ncol(CV.IVs.test)  != n) CV.IVs.test  <- t(CV.IVs.test)
        if (ncol(CV.IVs.test)  != n) CV.IVs.test  <- t(CV.IVs.test)
      }
      
      threshold_results_1 <- vector(mode = "list", length = length(c(1:l, length(IVs.train[, 1]))))
      nns.cv.1 <- numeric()
      
      q     <- length(IVs.train[, 1])
      Kcand <- c(1:l, q)
      
      # build *aligned* dummy matrices for TRAIN and TEST in one shot
      build_design_pair <- function(train_df, test_df) {
        tr <- as.data.frame(train_df, stringsAsFactors = TRUE)
        te <- as.data.frame(test_df,  stringsAsFactors = TRUE)
        
        # If either has no names, synthesize consistent names
        if (is.null(names(tr)) || anyNA(names(tr))) names(tr) <- paste0("X", seq_len(ncol(tr)))
        if (is.null(names(te)) || anyNA(names(te))) names(te) <- paste0("X", seq_len(ncol(te)))
        
        # 1) take the UNION of names
        alln <- union(names(tr), names(te))
        
        # 2) add any missing columns as NA (they’ll dummy to zeros after factor -> dummy)
        add_missing <- function(df, alln) {
          miss <- setdiff(alln, names(df))
          for (m in miss) df[[m]] <- NA
          # reorder to the common order
          df[, alln, drop = FALSE]
        }
        tr <- add_missing(tr, alln)
        te <- add_missing(te, alln)
        
        # proceed with factor_2_dummy_FR on the combined columns
        pieces_tr <- list(); pieces_te <- list()
        for (nm in names(tr)) {
          combo <- c(tr[[nm]], te[[nm]])
          block <- factor_2_dummy_FR(combo)
          if (is.null(dim(block))) block <- matrix(as.numeric(block), ncol = 1L)
          ntr <- NROW(tr)
          pieces_tr[[nm]] <- block[seq_len(ntr), , drop = FALSE]
          pieces_te[[nm]] <- block[(ntr + 1L):(ntr + NROW(te)), , drop = FALSE]
        }
        Xtr <- do.call(cbind, pieces_tr); storage.mode(Xtr) <- "double"
        Xte <- do.call(cbind, pieces_te); storage.mode(Xte) <- "double"
        list(Xtr = Xtr, Xte = Xte)
      }
      
      pred_path_small <- NULL   # |Xtest| x l (k = 1..l)
      pred_q          <- NULL   # |Xtest| vector (k = q)
      
      for (i in Kcand) {
        index <- which(Kcand == i)[1L]
        
        if (index == 1L) {
          # One NNS.reg call per fold to get fitted y.hat and
          # a baseline prediction for threshold optimisation
          setup <- suppressWarnings(
            NNS.reg(
              CV.IVs.train, CV.DV.train,
              point.est = CV.IVs.test,
              plot = FALSE, residual.plot = FALSE,
              n.best = 1, order = order,
              type = type, factor.2.dummy = TRUE,
              dist = dist, ncores = ncores,
              point.only = FALSE, smooth = smoothness
            )
          )
          
          # y.hat for each TRAINING point – used as the value aggregated by k-NN
          yhat_vec <- as.numeric(setup$Fitted.xy$y.hat)
          if (length(yhat_vec) != NROW(CV.IVs.train)) {
            stop("Internal: length(yhat_vec) must equal nrow(CV.IVs.train).")
          }
          
          # Build aligned numeric design matrices directly from TRAIN & TEST IVs
          design   <- build_design_pair(CV.IVs.train, CV.IVs.test)
          RPM_num   <- design$Xtr   # one row per training observation
          Xtest_num <- design$Xte   # one row per test observation
          
          # Sanity guards (fail fast instead of crashing in C++)
          stopifnot(
            is.matrix(RPM_num), is.double(RPM_num),
            is.matrix(Xtest_num), is.double(Xtest_num),
            nrow(RPM_num) == length(yhat_vec),
            ncol(RPM_num) == ncol(Xtest_num)
          )
          
          # Index == 1: original point estimates from NNS.reg for thresholding
          predicted <- setup$Point.est
          predicted[is.na(predicted)] <- mean(predicted, na.rm = TRUE)
          
          if (!is.null(type)) {
            pred_matrix <- if (length(unique(predicted)) == 1L) {
              matrix(replicate(100L, predicted), nrow = length(predicted))
            } else {
              sapply(seq(.01, .99, .01),
                     function(z)
                       ifelse(predicted %% 1 < z,
                              floor(predicted),
                              ceiling(predicted)))
            }
            threshold_results_1[[index]] <-
              seq(.01, .99, .01)[
                which.max(apply(pred_matrix, 2L,
                                function(z) mean(z == as.numeric(actual))))
              ]
            predicted <- ifelse(predicted %% 1 < threshold_results_1[[index]],
                                floor(predicted), ceiling(predicted))
          }
          
          # Precompute only what we need via C++ for k = 1..min(l, n_train)
          kmax_use <- min(l, nrow(RPM_num))
          pred_path_small <- NNS_distance_path_cpp(
            RPM   = RPM_num,
            yhat  = yhat_vec,
            Xtest = Xtest_num,
            kmax  = kmax_use,
            is_class = !is.null(type)
          )
          
          # Cache the "all training points" prediction for k = min(q, n_train)
          n_train_fold <- nrow(RPM_num)
          if (q > ncol(pred_path_small) || q > kmax_use) {
            pred_q <- NNS_distance_bulk_cpp(
              RPM   = RPM_num,
              yhat  = yhat_vec,
              Xtest = Xtest_num,
              k     = min(q, n_train_fold),
              is_class = !is.null(type)
            )
          } else {
            pred_q <- pred_path_small[, q, drop = TRUE]
          }
          
        } else {
          if (!is.null(dim(CV.IVs.train)) && ncol(CV.IVs.train) > 1) {
            if (i <= ncol(pred_path_small)) {
              predicted <- pred_path_small[, i, drop = TRUE]
            } else if (i == q) {
              predicted <- pred_q
            } else {
              predicted <- NNS_distance_bulk_cpp(
                RPM   = RPM_num,
                yhat  = yhat_vec,
                Xtest = Xtest_num,
                k     = min(i, nrow(RPM_num)),
                is_class = !is.null(type)
              )
            }
          } else {
            predicted <- suppressWarnings(
              NNS.reg(
                CV.IVs.train, CV.DV.train,
                point.est = if (is.null(dim(CV.IVs.test))) unlist(CV.IVs.test) else CV.IVs.test,
                plot = FALSE, residual.plot = FALSE,
                n.best = i, order = order, ncores = ncores,
                type = type, factor.2.dummy = TRUE,
                dist = dist, point.only = TRUE, smooth = smoothness
              )$Point.est
            )
          }
          
          if (!is.null(type)) {
            pred_matrix <- if (length(unique(predicted)) == 1) {
              matrix(replicate(100, predicted), nrow = length(predicted))
            } else {
              sapply(seq(.01, .99, .01),
                     function(z) ifelse(predicted %% 1 < z, floor(predicted), ceiling(predicted)))
            }
            z <- apply(pred_matrix, 2, function(z) mean(z == as.numeric(actual)))
            threshold_results_1[[index]] <- seq(.01, .99, .01)[as.integer(median(which(z == max(z))))]
            predicted <- ifelse(predicted %% 1 < threshold_results_1[[index]],
                                floor(predicted), ceiling(predicted))
          }
        }
        
        # objective at this n.best
        nns.cv.1[index] <- eval(obj.fn)
        
        # print n.best + eval(obj.fn) when status = TRUE
        if (status) {
          message(sprintf(
            "Current NNS.reg(. , n.best = %d ) | eval(obj.fn) = %.6f | MAX Iterations Remaining = %d",
            i,
            nns.cv.1[index],
            length(Kcand) - index
          ))
        }
        
        if (length(na.omit(nns.cv.1)) > 3) {
          if (objective == 'min') nns.cv.1[is.na(nns.cv.1)] <- max(na.omit(nns.cv.1)) else
            nns.cv.1[is.na(nns.cv.1)] <- min(na.omit(nns.cv.1))
          if (objective == 'min' && nns.cv.1[index] >= nns.cv.1[index - 1] && nns.cv.1[index] >= nns.cv.1[index - 2]) break
          if (objective == 'max' && nns.cv.1[index] <= nns.cv.1[index - 1] && nns.cv.1[index] <= nns.cv.1[index - 2]) break
        }
      }
      
      ks <- Kcand[!is.na(nns.cv.1)]
      if (objective == 'min') {
        k <- ks[which.min(na.omit(nns.cv.1))]; nns.cv.1 <- min(na.omit(nns.cv.1))
      } else {
        k <- ks[which.max(na.omit(nns.cv.1))]; nns.cv.1 <- max(na.omit(nns.cv.1))
      }
      
      best.k[[b]]      <- k
      best.nns.cv[[b]] <- if (!is.null(type)) min(max(nns.cv.1, 0), 1) else nns.cv.1
      
      if (b == folds) {
        ks_tab <- table(unlist(best.k))
        best.k <- mode_class(as.numeric(rep(names(ks_tab), as.numeric(unlist(ks_tab)))))
        best.k <- ifelse(best.k %% 1 < 0.5, floor(best.k), ceiling(best.k))
        
        # When stacking with Method 2, fit Method 1 on cbind(X*, X*) using full data
        if (stack && identical(sort(method), c(1, 2)) &&
            exists("xstar_full_train") && !anyNA(xstar_full_train)) {
          IVs.train.m1 <- data.frame(Xstar  = xstar_full_train,
                                     Xstar2 = xstar_full_train)
          IVs.test.m1  <- data.frame(Xstar  = xstar_full_test,
                                     Xstar2 = xstar_full_test)
          nns.method.1 <- suppressWarnings(
            NNS.reg(
              IVs.train.m1, DV.train,
              point.est = IVs.test.m1,
              plot = FALSE, n.best = best.k, order = order, ncores = ncores,
              type = type, point.only = FALSE, confidence.interval = pred.int, smooth = smoothness
            )
          )
        } else if (length(relevant_vars) > 1) {
          nns.method.1 <- suppressWarnings(
            NNS.reg(
              IVs.train[, relevant_vars], DV.train,
              point.est = IVs.test[, relevant_vars],
              plot = FALSE, n.best = best.k, order = order, ncores = ncores,
              type = type, point.only = FALSE, confidence.interval = pred.int, smooth = smoothness
            )
          )
        } else {
          nns.method.1 <- suppressWarnings(
            NNS.reg(
              IVs.train[, relevant_vars], DV.train,
              point.est = unlist(IVs.test[, relevant_vars]),
              plot = FALSE, n.best = best.k, order = order, ncores = ncores,
              type = type, point.only = FALSE, confidence.interval = pred.int, smooth = smoothness
            )
          )
        }
        
        actual    <- nns.method.1$Fitted.xy$y
        predicted <- nns.method.1$Fitted.xy$y.hat
        
        best.nns.cv <- eval(obj.fn)
        
        pred.int.1   <- nns.method.1$pred.int
        nns.method.1 <- nns.method.1$Point.est
        
        if (!is.null(type) && !is.null(nns.method.1)) {
          threshold_results_1 <- mean(unlist(threshold_results_1))
          nns.method.1 <- ifelse(nns.method.1 %% 1 < threshold_results_1,
                                 floor(nns.method.1), ceiling(nns.method.1))
          nns.method.1 <- pmin(nns.method.1, max(as.numeric(DV.train)))
          nns.method.1 <- pmax(nns.method.1, min(as.numeric(DV.train)))
        }
      }
      
    } else {
      test.set.1          <- NULL
      best.k              <- NA
      nns.method.1        <- NA
      threshold_results_1 <- NA
      if (objective == 'min') { best.nns.cv <- Inf } else { best.nns.cv <- -Inf }
    } # end: 1 %in% method
    
  } # errors (b) loop
  
  
  ### Weights for combining NNS techniques
  best.nns.cv[best.nns.cv == 0] <- 1e-10
  best.nns.ord[best.nns.ord == 0] <- 1e-10
  
  if(objective=="min"){
    weights <- c(max(1e-10, 1 / best.nns.cv^2), max(1e-10, 1 / best.nns.ord^2))
  } else {
    weights <- c(max(1e-10, best.nns.cv^2), max(1e-10, best.nns.ord^2))
  }
  
  
  weights <- pmax(weights, c(0, 0))
  weights[!(c(1, 2) %in% method)] <- 0
  weights[is.nan(weights)] <- 0
  weights[is.infinite(weights)] <- 0
  
  if(sum(weights)>0)  weights <- weights / sum(weights) else weights <- c(.5, .5)
  
  if(!is.null(type)) probability.threshold <-  mean(c(threshold_results_1, threshold_results_2), na.rm = TRUE) else probability.threshold <- .5
  
  if(identical(sort(method),c(1,2))){
    if(anyNA(nns.method.1)){
      na.1.index <- which(is.na(nns.method.1))
      nns.method.1[na.1.index] <- nns.method.2[na.1.index]
    }
    if(anyNA(nns.method.2)){
      na.2.index <- which(is.na(nns.method.2))
      nns.method.2[na.2.index] <- nns.method.1[na.2.index]
    }
    
    estimates <- (weights[1] * nns.method.1 + weights[2] * nns.method.2)
    if(!is.null(pred.int)) stacked.pred.int <- (weights[1] * pred.int.1 + weights[2] * pred.int.2) else stacked.pred.int <- NULL
    
    if(!is.null(type)){
      estimates <- ifelse(estimates%%1 < probability.threshold, floor(estimates), ceiling(estimates))
      estimates <- pmin(estimates, max(as.numeric(DV.train)))
      estimates <- pmax(estimates, min(as.numeric(DV.train)))
      
      if(!is.null(pred.int)) stacked.pred.int <- data.table::data.table(apply(stacked.pred.int, 2, function(x) ifelse(x%%1 <0.5, floor(x), ceiling(x))))
    }
  } else {
    if(method==1){
      estimates <- nns.method.1
      pred.int.2 <- NULL
      stacked.pred.int <- pred.int.1
    } else {
      if(method==2){
        estimates <- nns.method.2
        pred.int.1 <- NULL
        stacked.pred.int <- pred.int.2
      }
    }
  }
  
  
  if(is.null(probability.threshold)) probability.threshold <- .5
  
  return(list(OBJfn.reg = best.nns.cv,
              NNS.reg.n.best = best.k,
              probability.threshold = probability.threshold,
              OBJfn.dim.red = best.nns.ord,
              NNS.dim.red.threshold = nns.ord.threshold,
              reg = nns.method.1,
              reg.pred.int = pred.int.1,
              dim.red = nns.method.2,
              dim.red.pred.int = pred.int.2,
              stack = estimates,
              pred.int = stacked.pred.int))
  
}