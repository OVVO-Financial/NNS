NNS.ANOVA.bin <- function(control, treatment,
                          means.only = FALSE,
                          medians = FALSE,
                          mean.of.means = NULL,
                          upper.25.target = NULL,
                          lower.25.target = NULL,
                          upper.125.target = NULL,
                          lower.125.target = NULL,
                          confidence.interval = NULL, 
                          tails = NULL, 
                          plot = TRUE, 
                          par = NULL,
                          n_boot = 1000) {
  
  # Calculate grand statistic if not provided
  if(is.null(mean.of.means)) {
    if(medians) {
      mean.of.means <- median(c(median(control), median(treatment)))
    } else {
      mean.of.means <- mean(c(mean(control), mean(treatment)))
    }
  }
  
  # Calculate partial moment targets if not provided
  if(is.null(upper.25.target) && is.null(lower.25.target)) {
    upper.25.target <- mean(c(UPM.VaR(0.25, 1, control), UPM.VaR(0.25, 1, treatment)))
    lower.25.target <- mean(c(LPM.VaR(0.25, 1, control), LPM.VaR(0.25, 1, treatment)))
    upper.125.target <- mean(c(UPM.VaR(0.125, 1, control), UPM.VaR(0.125, 1, treatment)))
    lower.125.target <- mean(c(LPM.VaR(0.125, 1, control), LPM.VaR(0.125, 1, treatment)))
  }
  
  # Calculate partial moment ratios
  if(medians) {
    # Median: Use degree 0 (frequency-based)
    LPM_ratio.1 <- LPM.ratio(0, mean.of.means, control)
    LPM_ratio.2 <- LPM.ratio(0, mean.of.means, treatment)
  } else {
    # Mean: Use degree 1 (moment-based)
    LPM_ratio.1 <- LPM.ratio(1, mean.of.means, control) / 
      (LPM.ratio(1, mean.of.means, control) + UPM.ratio(1, mean.of.means, control))
    LPM_ratio.2 <- LPM.ratio(1, mean.of.means, treatment) / 
      (LPM.ratio(1, mean.of.means, treatment) + UPM.ratio(1, mean.of.means, treatment))
  }
  
  # Calculate partial moment ratios at thresholds
  Upper_25_ratio.1 <- UPM.ratio(1, upper.25.target, control)
  Upper_25_ratio.2 <- UPM.ratio(1, upper.25.target, treatment)
  Upper_25_ratio <- mean(c(Upper_25_ratio.1, Upper_25_ratio.2))
  
  Lower_25_ratio.1 <- LPM.ratio(1, lower.25.target, control)
  Lower_25_ratio.2 <- LPM.ratio(1, lower.25.target, treatment)
  Lower_25_ratio <- mean(c(Lower_25_ratio.1, Lower_25_ratio.2))
  
  Upper_125_ratio.1 <- UPM.ratio(1, upper.125.target, control)
  Upper_125_ratio.2 <- UPM.ratio(1, upper.125.target, treatment)
  Upper_125_ratio <- mean(c(Upper_125_ratio.1, Upper_125_ratio.2))
  
  Lower_125_ratio.1 <- LPM.ratio(1, lower.125.target, control)
  Lower_125_ratio.2 <- LPM.ratio(1, lower.125.target, treatment)
  Lower_125_ratio <- mean(c(Lower_125_ratio.1, Lower_125_ratio.2))
  
  # Calculate CDF deviations
  MAD.CDF <- min(0.5, max(c(abs(LPM_ratio.1 - 0.5), abs(LPM_ratio.2 - 0.5))))
  upper.25.CDF <- min(0.25, max(c(abs(Upper_25_ratio.1 - 0.25), abs(Upper_25_ratio.2 - 0.25))))
  lower.25.CDF <- min(0.25, max(c(abs(Lower_25_ratio.1 - 0.25), abs(Lower_25_ratio.2 - 0.25))))
  upper.125.CDF <- min(0.125, max(c(abs(Upper_125_ratio.1 - 0.125), abs(Upper_125_ratio.2 - 0.125))))
  lower.125.CDF <- min(0.125, max(c(abs(Lower_125_ratio.1 - 0.125), abs(Lower_125_ratio.2 - 0.125))))
  
  # Calculate certainty statistic
  if(means.only) {
    NNS.ANOVA.rho <- ((0.5 - MAD.CDF)^2) / 0.25 
  } else {
    NNS.ANOVA.rho <- sum(
      c( ((0.5 - MAD.CDF)^2) / 0.25,
         0.5 * (((0.25 - upper.25.CDF)^2) / (0.25^2)),
         0.5 * (((0.25 - lower.25.CDF)^2) / (0.25^2)),
         0.25 * (((0.125 - upper.125.CDF)^2) / (0.125^2)),
         0.25 * (((0.125 - lower.125.CDF)^2) / (0.125^2))
      )) / 2.5
  }
  
  # Population size adjustment
  pop.adjustment <- ((length(control) + length(treatment) - 2) / 
                       (length(control) + length(treatment)))^2
  
  # Plotting
  if(plot) {
    if(is.null(par)) {
      original.par <- par(no.readonly = TRUE)
      on.exit(par(original.par))
    }
    
    boxplot(list(control, treatment), 
            names = c("Control", "Treatment"), 
            horizontal = TRUE, 
            main = "NNS ANOVA and Effect Size", 
            col = c("grey", "white"),
            cex.axis = 0.75)
    
    abline(v = mean.of.means, col = "red", lwd = 4)
    if(medians) {
      mtext("Grand Median", side = 3, col = "red", at = mean.of.means)
    } else {
      mtext("Grand Mean", side = 3, col = "red", at = mean.of.means)
    }
  }
  
  # Confidence interval and effect size calculation
  if(!is.null(confidence.interval)) {
    # Validate tails parameter
    if(is.null(tails)) stop("tails must be specified with confidence.interval")
    tails <- match.arg(tails, c("both", "left", "right"))
    
    # Bootstrap both groups
    control_boot <- matrix(sample(control, size = n_boot * length(control), replace = TRUE), 
                           nrow = length(control))
    treatment_boot <- matrix(sample(treatment, size = n_boot * length(treatment), replace = TRUE), 
                             nrow = length(treatment))
    
    if(medians) {
      control_stats <- apply(control_boot, 2, median)
      treatment_stats <- apply(treatment_boot, 2, median)
    } else {
      control_stats <- colMeans(control_boot)
      treatment_stats <- colMeans(treatment_boot)
    }
    
    # Calculate confidence bounds
    alpha <- if(tails == "both") (1 - confidence.interval)/2 else 1 - confidence.interval
    
    if(tails %in% c("both", "right")) {
      control_upper <- UPM.VaR(alpha, 0, control_stats)
      treatment_upper <- UPM.VaR(alpha, 0, treatment_stats)
    }
    
    if(tails %in% c("both", "left")) {
      control_lower <- LPM.VaR(alpha, 0, control_stats)
      treatment_lower <- LPM.VaR(alpha, 0, treatment_stats)
    }
    
    # Calculate conservative effect size bounds
    if(tails == "both") {
      min_effect <- treatment_lower - control_upper  # Minimum plausible effect
      max_effect <- treatment_upper - control_lower  # Maximum plausible effect
    } else if(tails == "left") {
      min_effect <- treatment_lower - control_upper
      max_effect <- Inf
    } else if(tails == "right") {
      min_effect <- -Inf
      max_effect <- treatment_upper - control_lower
    }
    
    
    # Add confidence bounds to plot
    if(plot) {
      # Add CI bounds with solid lines
      if(tails %in% c("both", "right")) {
        abline(v = control_upper, col = "blue", lwd = 2)
        abline(v = treatment_upper, col = "darkblue", lwd = 2)
      }
      
      if(tails %in% c("both", "left")) {
        abline(v = control_lower, col = "green", lwd = 2)
        abline(v = treatment_lower, col = "darkgreen", lwd = 2)
      }
      
      # Add separate legends for lower and upper bounds
      ci_label <- paste0(confidence.interval * 100, "% CI")
      
      if(tails %in% c("both", "left")) {
        legend("topleft", 
               legend = c(paste("Control Lower", ci_label), 
                          paste("Treatment Lower", ci_label)),
               col = c("green", "darkgreen"),
               lty = 1,
               lwd = 2,
               cex = 0.8,
               bty = "n")
      }
      
      if(tails %in% c("both", "right")) {
        legend("topright", 
               legend = c(paste("Control Upper", ci_label), 
                          paste("Treatment Upper", ci_label)),
               col = c("blue", "darkblue"),
               lty = 1,
               lwd = 2,
               cex = 0.8,
               bty = "n")
      }
    }
  
    
    # Return results with effect size bounds
    result <- list(
      Control = if(medians) median(control) else mean(control),
      Treatment = if(medians) median(treatment) else mean(treatment),
      Grand_Statistic = mean.of.means,
      Control_CDF = LPM_ratio.1,
      Treatment_CDF = LPM_ratio.2,
      Certainty = min(1, NNS.ANOVA.rho * pop.adjustment),
      Effect_Size_LB = min_effect,
      Effect_Size_UB = max_effect,
      Confidence_Level = confidence.interval
    )
    return(result)
    
  } else {
    # Return basic results without effect size bounds
    result <- list(
      Control = if(medians) median(control) else mean(control),
      Treatment = if(medians) median(treatment) else mean(treatment),
      Grand_Statistic = mean.of.means,
      Control_CDF = LPM_ratio.1,
      Treatment_CDF = LPM_ratio.2,
      Certainty = min(1, NNS.ANOVA.rho * pop.adjustment)
    )
    return(result)
  }
}