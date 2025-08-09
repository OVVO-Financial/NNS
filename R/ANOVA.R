#' NNS ANOVA: Nonparametric Analysis of Variance
#' 
#' Performs a distribution-free ANOVA using partial moment statistics to evaluate differences between control and treatment groups. Returns a certainty metric (0-1) indicating separation likelihood rather than traditional p-values. Includes bootstrapped effect size bounds.
#' 
#' @param control Numeric vector of control group observations
#' @param treatment Numeric vector of treatment group observations
#' @param means.only Logical; \code{FALSE} (default) uses full distribution analysis. Set \code{TRUE} for mean-only comparison
#' @param medians Logical; \code{FALSE} (default) uses means. Set \code{TRUE} for median-based analysis
#' @param confidence.interval Numeric [0,1]; confidence level for effect size bounds (e.g., 0.95)
#' @param tails Character; specifies CI tail(s): "both", "left", or "right"
#' @param pairwise logical; \code{FALSE} (default) Returns pairwise certainty tests when set to \code{pairwise = TRUE}.
#' @param robust logical; \code{FALSE} (default) Generates 100 independent random permutations to test results, and returns / plots 95 percent confidence intervals along with robust central tendency of all results for pairwise analysis only.
#' @param plot Logical; \code{TRUE} (default) generates distribution plot
#' 
#' @return Returns a list containing:
#' \itemize{
#'   \item \code{Control_Statistic}: Mean/median of control group
#'   \item \code{Treatment_Statistic}: Mean/median of treatment group
#'   \item \code{Grand_Statistic}: Grand mean/median
#'   \item \code{Control_CDF}: CDF value at grand statistic (control)
#'   \item \code{Treatment_CDF}: CDF value at grand statistic (treatment)
#'   \item \code{Certainty}: Separation certainty (0-1)
#'   \item \code{Effect_Size_LB}: Lower bound of treatment effect (if CI requested)
#'   \item \code{Effect_Size_UB}: Upper bound of treatment effect (if CI requested)
#'   \item \code{Confidence_Level}: Confidence level used (if CI requested)
#' }
#' 
#'
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments" (ISBN: 1490523995)
#'
#' Viole, F. (2017) "Continuous CDFs and ANOVA with NNS"  \doi{10.2139/ssrn.3007373}
#'
#' @examples
#'  \dontrun{
#' ### Binary analysis and effect size
#' set.seed(123)
#' x <- rnorm(100) ; y <- rnorm(100)
#' NNS.ANOVA(control = x, treatment = y)
#'
#' ### Two variable analysis with no control variable
#' A <- cbind(x, y)
#' NNS.ANOVA(A)
#' 
#' ### Medians test
#' NNS.ANOVA(A, means.only = TRUE, medians = TRUE)
#'
#' ### Multiple variable analysis with no control variable
#' set.seed(123)
#' x <- rnorm(100) ; y <- rnorm(100) ; z <- rnorm(100)
#' A <- cbind(x, y, z)
#' NNS.ANOVA(A)
#' 
#' ### Different length vectors used in a list
#' x <- rnorm(30) ; y <- rnorm(40) ; z <- rnorm(50)
#' A <- list(x, y, z)
#' NNS.ANOVA(A)
#' }
#' @export


NNS.ANOVA <- function(
    control,
    treatment,
    means.only = FALSE,
    medians = FALSE,
    confidence.interval = 0.95,
    tails = "Both",
    pairwise = FALSE,
    plot = TRUE,
    robust = FALSE
){
  tails <- tolower(tails)
  if(!any(tails%in%c("left","right","both"))){
    stop("Please select tails from 'left', 'right', or 'both'")
  }
  if(!missing(treatment)){
    # with treatment
    if(any(class(control)%in%c("tbl","data.table"))) control <- as.vector(unlist(control))
    
    if(any(class(treatment)%in%c("tbl","data.table"))) treatment <- as.vector(unlist(treatment))
    
    if(robust){
      l <- min(length(treatment), length(control))
      treatment_p <- replicate(100, sample.int(l, replace = TRUE))
      treatment_matrix <- matrix(treatment[treatment_p], ncol = dim(treatment_p)[2], byrow = F)
      treatment_matrix <- cbind(treatment[treatment_p[,1]], treatment_matrix[, rev(seq_len(ncol(treatment_matrix)))])
      control_matrix <- matrix(control[treatment_p], ncol = dim(treatment_p)[2], byrow = F)
      full_matrix <- cbind(control[treatment_p[,1]], control_matrix, treatment[treatment_p[,1]], treatment_matrix)
      
      nns.certainties <- sapply(
        1:ncol(control_matrix), 
        function(g) NNS.ANOVA.bin(control_matrix[,g], treatment_matrix[,g], means.only = means.only, medians = medians, plot = FALSE)$Certainty
      )
      
      cer_lower_CI <- LPM.VaR(.025, 1, nns.certainties)
      cer_upper_CI <- UPM.VaR(.025, 1, nns.certainties)
      
      robust_estimate <- gravity(nns.certainties)
      
      if(plot){
        original.par <- par(no.readonly = TRUE)
        par(mfrow = c(1, 2))
        hist(nns.certainties, main = "NNS Certainty")
        abline(v = robust_estimate, col = 'red', lwd = 3)
        mtext("Robust Certainty Estimate", side = 3, col = "red", at = robust_estimate)
        abline(v =  cer_lower_CI, col = "red", lwd = 2, lty = 3)
        abline(v =  cer_upper_CI , col = "red", lwd = 2, lty = 3)
      }
      return(
        as.list(c(unlist(
          NNS.ANOVA.bin(
            control, 
            treatment,
            means.only = means.only,
            medians = medians,
            confidence.interval = confidence.interval, 
            plot = plot, 
            tails = tails, 
            par = original.par
          )
        ),
        "Robust Certainty Estimate" = robust_estimate,
        "Lower 95% CI" = cer_lower_CI,
        "Upper 95% CI" = cer_upper_CI))
      )
    } else {
      return(
        NNS.ANOVA.bin(
          control, 
          treatment, 
          means.only = means.only,
          medians = medians,
          confidence.interval = confidence.interval, 
          plot = plot, 
          tails = tails
        )
      )
    }
  }
  
  # without treatment
  if(is.list(control)) n <- length(control) else n <- ncol(control)
  
  if(is.null(n)) stop("supply both 'control' and 'treatment' or a matrix-like 'control', or a list 'control'")
  
  if(n == 1) stop("supply both 'control' and 'treatment' or a matrix-like 'control', or a list 'control'")
  
  if(n >= 2){
    if(any(class(control) %in% c("tbl","data.table"))){
      A <- as.data.frame(control)
    } else {
      if(any(class(control) %in% "list")) A <- do.call(cbind, lapply(control, `length<-`, max(lengths(control)))) else A <- control
    }
  } else {
    A <- control
  }
  
  if(medians) mean.of.means <- mean(apply(A, 2, function(i) median(i, na.rm = TRUE))) else mean.of.means <- mean(colMeans(A, na.rm = T))
  
  if(!pairwise){
    #Continuous CDF for each variable from Mean of Means
    if(medians){
      LPM_ratio <- sapply(1:n, function(b) LPM.ratio(0, mean.of.means, na.omit(unlist(A[ , b]))))
    } else {
      LPM_ratio <- sapply(1:n, function(b) LPM.ratio(1, mean.of.means, na.omit(unlist(A[ , b]))))
    }
    
    lower.25.target  <- mean(sapply(1:n, function(i) LPM.VaR(.25,  1, na.omit(unlist(A[,i])))))
    upper.25.target  <- mean(sapply(1:n, function(i) UPM.VaR(.25,  1, na.omit(unlist(A[,i])))))
    lower.125.target <- mean(sapply(1:n, function(i) LPM.VaR(.125, 1, na.omit(unlist(A[,i])))))
    upper.125.target <- mean(sapply(1:n, function(i) UPM.VaR(.125, 1, na.omit(unlist(A[,i])))))
    
    raw.certainties <- list(n - 1)
    for(i in 1:(n - 1)){
      raw.certainties[[i]] <- sapply(
        (i + 1) : n, 
        function(b) NNS.ANOVA.bin(
          na.omit(unlist(A[ , i])),
          na.omit(unlist(A[ , b])),
          means.only = means.only,
          medians = medians,
          mean.of.means = mean.of.means,
          upper.25.target = upper.25.target,
          lower.25.target = lower.25.target,
          upper.125.target = upper.125.target,
          lower.125.target = lower.125.target,
          plot = FALSE
        )$Certainty
      )
    }
    
    #Certainty associated with samples
    NNS.ANOVA.rho <- mean(unlist(raw.certainties))
    
    #Graphs
    if(plot){
      boxplot(
        A, 
        las = 2,  
        ylab = "Variable", 
        horizontal = TRUE, 
        main = "NNS ANOVA", 
        col = c('steelblue', rainbow(n - 1))
      )
      #For ANOVA Visualization
      abline(v = mean.of.means, col = "red", lwd = 4)
      if(medians) mtext("Grand Median", side = 3,col = "red", at = mean.of.means) else mtext("Grand Mean", side = 3,col = "red", at = mean.of.means)
    }
    return(c("Certainty" = NNS.ANOVA.rho))
  }
  
  raw.certainties <- list(n - 1)
  for(i in 1:(n - 1)){
    raw.certainties[[i]] <- sapply(
      (i + 1) : n, 
      function(b) NNS.ANOVA.bin(na.omit(unlist(A[ , i])), na.omit(unlist(A[ , b])), means.only = means.only, medians = medians, plot = FALSE)$Certainty
    )
  }
  
  certainties <- matrix(NA, n, n)
  certainties[lower.tri(certainties, diag = FALSE)] <- unlist(raw.certainties)
  diag(certainties) <- 1
  certainties <- pmax(certainties, t(certainties), na.rm = TRUE)
  colnames(certainties) <- rownames(certainties) <- colnames(A)
  
  if(plot){
    boxplot(
      A, 
      las = 2, 
      ylab = "Variable", 
      horizontal = TRUE, 
      main = "ANOVA", 
      col = c('steelblue', rainbow(n - 1))
    )
    abline(v = mean.of.means, col = "red", lwd = 4)
    if(medians) mtext("Grand Median", side = 3,col = "red", at = mean.of.means) else mtext("Grand Mean", side = 3,col = "red", at = mean.of.means)
  }
  return(certainties)
}