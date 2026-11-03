#' NNS ANOVA: Nonparametric Analysis of Variance
#'
#' Compares control and treatment groups through benchmark-relative partial
#' moments. The chosen mode compares means, medians, or a distribution profile
#' with central and tail benchmarks, without selecting a parametric family.
#'
#' The key output, \code{Certainty}, is a descriptive agreement index in
#' \eqn{[0, 1]} for the chosen target:
#' \itemize{
#'   \item \code{means.only = TRUE, medians = FALSE}: compares group means
#'         using degree-one deviation-area ratios at the pooled mean.
#'   \item \code{means.only = TRUE, medians = TRUE}: compares medians using
#'         the distance of 0.5 from each interval \eqn{[F(t-), F(t)]} at the
#'         pooled median \eqn{t}. This interval contains 0.5 exactly when
#'         \eqn{t} is a median, including tied and discrete data.
#'   \item \code{means.only = FALSE}: adds the 12.5\% and 25\% lower and
#'         upper deviation-area benchmarks to the chosen central comparison.
#'         This examines distributional shape through a finite profile;
#'         agreement at these benchmarks does not prove equality of laws.
#' }
#'
#' At a constant group's only value both deviation areas vanish. For agreement
#' comparisons this satisfies the expectile balance at every level, so its
#' area discrepancy is zero there. Returned CDF values retain the public
#' partial-moment step convention; other benchmarks and nonconstant samples
#' use the ordinary ratio discrepancy.
#'
#' The raw score is 1 at maximal benchmark agreement. The returned score is
#' multiplied by \eqn{((n_1 + n_2 - 2)/(n_1 + n_2))^2}, so the finite-sample
#' maximum is below 1. \code{Certainty} is not a calibrated probability,
#' a p-value, or statistical power, and has no universal significance cutoff.
#' Equal population means need not produce the same score distribution when
#' sample sizes or spreads differ. \code{robust = TRUE} describes sampling
#' variability of the score; it does not generate a null reference distribution.
#'
#' Random relabelings of the pooled observations provide a reference under
#' identical distributions only when observations are exchangeable. They are
#' not a general unequal-variance equal-means test. Any inferential calibration
#' must match the chosen target and the sampling design; dependent observations
#' require a suitable resampling design.
#'
#' Effect-size bounds retain the conservative difference of the groups' separate
#' bootstrap endpoints. They are not quantiles of the bootstrap difference and
#' need not have exactly the requested nominal coverage. The independent-group
#' bootstrap assumes independent observations within and between groups.
#'
#' @param control Numeric vector of control group observations
#' @param treatment Numeric vector of treatment group observations
#' @param means.only Logical; \code{FALSE} (default) uses full distribution analysis. Set \code{TRUE} for the central comparison only (means, or medians when \code{medians = TRUE})
#' @param medians Logical; \code{FALSE} (default) uses means. Set \code{TRUE} for median-based analysis
#' @param confidence.interval Numeric [0,1]; confidence level for effect size bounds (e.g., 0.95)
#' @param tails Character; specifies CI tail(s): "both", "left", or "right"
#' @param pairwise logical; \code{FALSE} (default) Returns pairwise certainty tests when set to \code{pairwise = TRUE}.
#' @param robust logical; \code{FALSE} (default) Two-sample (\code{control}, \code{treatment}) path only: draws 100 independent bootstrap resamples of each group (with replacement, from its own full sample), and returns / plots the robust central tendency (\link{NNS.gravity}) and the \code{confidence.interval} bounds of the resampled \code{Certainty} values.  It is ignored when \code{control} is a matrix or list (multi-group, including \code{pairwise = TRUE}).
#' @param plot Logical; \code{TRUE} (default) generates distribution plot
#' 
#' @return Returns a list containing:
#' \itemize{
#'   \item \code{Control}: Mean/median of control group
#'   \item \code{Treatment}: Mean/median of treatment group
#'   \item \code{Grand_Statistic}: Mean/median of all observations pooled
#'   \item \code{Control_CDF}: CDF value at grand statistic (control)
#'   \item \code{Treatment_CDF}: CDF value at grand statistic (treatment)
#'   \item \code{Certainty}: Agreement index for the chosen central or
#'         distribution-profile comparison, with the finite-sample adjustment;
#'         not a calibrated probability.
#'   \item \code{Effect_Size_LB}: Lower bound of treatment effect (if confidence.interval requested)
#'   \item \code{Effect_Size_UB}: Upper bound of treatment effect (if confidence.interval requested)
#'   \item \code{Confidence_Level}: Confidence level used (if confidence.interval requested)
#' }
#' 
#' 
#'
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments" (ISBN: 1490523995, 2nd edition: \url{https://ovvo-financial.github.io/NNS/book/})
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
  # Standardize and validate tail-selection input used by CI calculations.
  # normalize tails to lower case for downstream functions
  tails <- tolower(tails)
  if (!any(tails %in% c("left","right","both"))) {
    stop("Please select tails from 'left', 'right', or 'both'")
  }
  
  if (!missing(treatment) && !is.null(treatment)) {
    # Two-sample path: directly compare control vs. treatment.
    # with treatment
    if (any(class(control)   %in% c("tbl","data.table"))) control   <- as.vector(unlist(control))
    if (any(class(treatment) %in% c("tbl","data.table"))) treatment <- as.vector(unlist(treatment))

    # Missing values are omitted, as in the multi-group path (na.omit per column).
    control   <- control[!is.na(control)]
    treatment <- treatment[!is.na(treatment)]

    if (robust) {
      # ---------------------------
      # Robust path: bootstrap each group
      # ---------------------------
      control_p        <- replicate(100, sample.int(length(control),   replace = TRUE))
      treatment_p      <- replicate(100, sample.int(length(treatment), replace = TRUE))
      control_matrix   <- matrix(control[control_p],     ncol = 100L, byrow = FALSE)
      treatment_matrix <- matrix(treatment[treatment_p], ncol = 100L, byrow = FALSE)
      
      nns.certainties <- sapply(
        1:ncol(control_matrix),
        function(g) NNS.ANOVA.bin(
          control_matrix[, g],
          treatment_matrix[, g],
          means.only = means.only,
          medians = medians,
          plot = FALSE
        )$Certainty
      )
      
      # Robust point estimate of certainty
      robust_estimate <- gravity(nns.certainties)
      
      # --- CI for robust certainty, honoring confidence.interval and tails
      alpha <- if (tails == "both") (1 - confidence.interval) / 2 else 1 - confidence.interval
      
      # Use degree = 0 to obtain empirical quantiles of the certainty samples
      cer_lower_CI <- if (tails %in% c("both","left"))  LPM.VaR(alpha, 0, nns.certainties) else NA_real_
      cer_upper_CI <- if (tails %in% c("both","right")) UPM.VaR(alpha, 0, nns.certainties) else NA_real_
      
      # optional plotting (restore red lines and red label)
      if (plot) {
        original.par <- par(no.readonly = TRUE); on.exit(par(original.par), add = TRUE)
        par(mfrow = c(1, 2))
        hist(nns.certainties, main = "NNS Certainty")
        abline(v = robust_estimate, col = "red", lwd = 3)
        if (tails %in% c("both","left"))  abline(v = cer_lower_CI, col = "red", lwd = 2, lty = 3)
        if (tails %in% c("both","right")) abline(v = cer_upper_CI, col = "red", lwd = 2, lty = 3)
        mtext("Robust Certainty Estimate", side = 3, at = robust_estimate, col = "red")
      }
      
      # Compute main ANOVA (non-robust) results but DO NOT pass 'par'
      base <- NNS.ANOVA.bin(
        control,
        treatment,
        means.only = means.only,
        medians = medians,
        confidence.interval = confidence.interval,
        plot = plot,
        tails = tails
      )
      
      lvl <- round(100 * confidence.interval)
      lower_label <- if (tails == "both") paste0("Lower ", lvl, "% CI") else paste0("Lower ", lvl, "% bound")
      upper_label <- if (tails == "both") paste0("Upper ", lvl, "% CI") else paste0("Upper ", lvl, "% bound")
      
      # Build return object to match prior shape (flatten base then append)
      ci_vals  <- numeric(0); ci_names <- character(0)
      if (tails %in% c("both","left"))  { ci_vals <- c(ci_vals, cer_lower_CI); ci_names <- c(ci_names, lower_label) }
      if (tails %in% c("both","right")) { ci_vals <- c(ci_vals, cer_upper_CI); ci_names <- c(ci_names, upper_label) }
      if (length(ci_vals)) names(ci_vals) <- ci_names
      
      return(
        as.list(c(
          unlist(base),
          "Robust Certainty Estimate" = robust_estimate,
          ci_vals
        ))
      )
      
    } else {
      # Fast/default path delegates to binary ANOVA implementation.
      # non-robust, binary ANOVA (already honors confidence.interval & tails)
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
  
  # ------------------------------
  # Without treatment (k >= 2 vars)
  # ------------------------------
  # Multi-sample path: treat each column/list element as a separate group.
  if (is.list(control)) n <- length(control) else n <- ncol(control)
  if (is.null(n) || n == 1)
    stop("supply both 'control' and 'treatment' or a matrix-like 'control', or a list 'control'")
  
  if (n >= 2) {
    if (any(class(control) %in% c("tbl","data.table"))) {
      # Keep tabular inputs as data.frame for stable column-wise operations.
      A <- as.data.frame(control)
    } else {
      if (any(class(control) %in% "list")) {
        # Pad unequal-length vectors to a rectangular matrix for pairwise scans.
        A <- do.call(cbind, lapply(control, `length<-`, max(lengths(control))))
      } else {
        A <- control
      }
    }
  } else {
    A <- control
  }
  
  if (medians) {
    mean.of.means <- median(unlist(A), na.rm = TRUE)
  } else {
    mean.of.means <- mean(unlist(A), na.rm = TRUE)
  }
  
  if (!pairwise) {
    # Aggregate mode: compute one global certainty across all group pairs.
    # Continuous CDF for each variable from grand statistic
    if (medians) {
      LPM_ratio <- sapply(1:n, function(b) LPM.ratio(0, mean.of.means, na.omit(unlist(A[, b]))))
    } else {
      LPM_ratio <- sapply(1:n, function(b) LPM.ratio(1, mean.of.means, na.omit(unlist(A[, b]))))
    }
    
    lower.25.target  <- mean(sapply(1:n, function(i) LPM.VaR(.25,  1, na.omit(unlist(A[,i])))))
    upper.25.target  <- mean(sapply(1:n, function(i) UPM.VaR(.25,  1, na.omit(unlist(A[,i])))))
    lower.125.target <- mean(sapply(1:n, function(i) LPM.VaR(.125, 1, na.omit(unlist(A[,i])))))
    upper.125.target <- mean(sapply(1:n, function(i) UPM.VaR(.125, 1, na.omit(unlist(A[,i])))))
    
    raw.certainties <- vector("list", n - 1)
    for (i in 1:(n - 1)) {
      # Collect upper-triangle pairwise certainties, then average below.
      raw.certainties[[i]] <- sapply(
        (i + 1):n,
        function(b) NNS.ANOVA.bin(
          na.omit(unlist(A[, i])),
          na.omit(unlist(A[, b])),
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
    
    # Certainty associated with samples
    NNS.ANOVA.rho <- mean(unlist(raw.certainties))
    
    # Graphs
    if (plot) {
      boxplot(
        A,
        las = 2,
        ylab = "Variable",
        horizontal = TRUE,
        main = "NNS ANOVA",
        col = c('steelblue', rainbow(n - 1))
      )
      abline(v = mean.of.means, col = "red", lwd = 4)
      if (medians) mtext("Grand Median", side = 3, col = "red", at = mean.of.means) else mtext("Grand Mean", side = 3, col = "red", at = mean.of.means)
    }
    return(c("Certainty" = NNS.ANOVA.rho))
  }
  
  # pairwise = TRUE: return symmetric matrix of certainties
  raw.certainties <- vector("list", n - 1)
  for (i in 1:(n - 1)) {
    raw.certainties[[i]] <- sapply(
      (i + 1):n,
      function(b) NNS.ANOVA.bin(
        na.omit(unlist(A[, i])),
        na.omit(unlist(A[, b])),
        means.only = means.only,
        medians = medians,
        plot = FALSE
      )$Certainty
    )
  }
  
  certainties <- matrix(NA_real_, n, n)
  certainties[lower.tri(certainties, diag = FALSE)] <- unlist(raw.certainties)
  diag(certainties) <- 1
  certainties <- pmax(certainties, t(certainties), na.rm = TRUE)
  colnames(certainties) <- rownames(certainties) <- colnames(A)
  
  if (plot) {
    boxplot(
      A,
      las = 2,
      ylab = "Variable",
      horizontal = TRUE,
      main = "ANOVA",
      col = c('steelblue', rainbow(n - 1))
    )
    abline(v = mean.of.means, col = "red", lwd = 4)
    if (medians) mtext("Grand Median", side = 3, col = "red", at = mean.of.means) else mtext("Grand Mean", side = 3, col = "red", at = mean.of.means)
  }
  return(certainties)
}
