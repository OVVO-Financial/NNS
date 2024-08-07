NNS.ANOVA.bin <- function(control, treatment,
                         means.only = FALSE,
                         mean.of.means = NULL,
                         upper.25.target = NULL,
                         lower.25.target = NULL,
                         upper.125.target = NULL,
                         lower.125.target = NULL,
                         confidence.interval = NULL, tails = NULL, plot = TRUE, par = NULL){

  if(is.null(upper.25.target) && is.null(lower.25.target)){
        mean.of.means <- mean(c(mean(control), mean(treatment)))
        upper.25.target <- mean(c(UPM.VaR(.25, 1, control), UPM.VaR(.25, 1, treatment)))
        lower.25.target <- mean(c(LPM.VaR(.25, 1, control), LPM.VaR(.25, 1, treatment)))
        upper.125.target <- mean(c(UPM.VaR(.125, 1, control), UPM.VaR(.125, 1, treatment)))
        lower.125.target <- mean(c(LPM.VaR(.125, 1, control), LPM.VaR(.125, 1, treatment)))
  }



  #Continuous CDF for each variable from Mean of Means
        LPM_ratio.1 <- LPM.ratio(1, mean.of.means, control)
        LPM_ratio.2 <- LPM.ratio(1, mean.of.means, treatment)

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


  #Continuous CDF Deviation from 0.5
        MAD.CDF <- min(0.5, max(c(abs(LPM_ratio.1 - 0.5), abs(LPM_ratio.2 - 0.5))))
        upper.25.CDF <- min(0.25, max(c(abs(Upper_25_ratio.1 - 0.25), abs(Upper_25_ratio.2 - 0.25))))
        lower.25.CDF <- min(0.25, max(c(abs(Lower_25_ratio.1 - 0.25), abs(Lower_25_ratio.2 - 0.25))))
        upper.125.CDF <- min(0.125, max(c(abs(Upper_125_ratio.1 - 0.125), abs(Upper_125_ratio.2 - 0.125))))
        lower.125.CDF <- min(0.125, max(c(abs(Lower_125_ratio.1 - 0.125), abs(Lower_125_ratio.2 - 0.125))))


  #Certainty associated with samples
        if(means.only) NNS.ANOVA.rho <- ((.5 - MAD.CDF)^2) / .25 
        else {
        NNS.ANOVA.rho <- sum(c( ((.5 - MAD.CDF)^2) / .25,
                                .5 * (( (.25 - upper.25.CDF)^2) / .25^2),
                                .5 * (( (.25 - lower.25.CDF)^2) / .25^2),
                                .25 * (( (.125 - upper.125.CDF)^2) / .125^2),
                                .25 * (( (.125 - lower.125.CDF)^2) / .125^2)
                             )) / 2.5

        }
        
    
        
        
        
        pop.adjustment <- ((length(control) + length(treatment) - 2) / (length(control)  + length(treatment) )) ^ 2

  #Graphs

        if(plot){
            if(is.null(par)) original.par <- par(no.readonly = TRUE) else original.par <- par

            boxplot(list(control, treatment), las = 2, names = c("Control", "Treatment"), horizontal = TRUE, main = "NNS ANOVA and Effect Size", col = c("grey", "white"), cex.axis = 0.75)

            #For ANOVA Visualization
            abline(v = mean.of.means, col = "red", lwd = 4)
            mtext("Grand Mean", side = 3, col = "red", at = mean.of.means)
        }

    if(is.null(confidence.interval)){

          return(list("Control Mean" = mean(control),
                "Treatment Mean" = mean(treatment),
                "Grand Mean" = mean.of.means,
                "Control CDF" = LPM_ratio.1,
                "Treatment CDF" = LPM_ratio.2,
                "Certainty" = min(1, NNS.ANOVA.rho * pop.adjustment)))
    } else {

        #Upper end of CDF confidence interval for control mean
        if(tails == "both") CI <- (1-confidence.interval)/2
        
        if(tails == "left" || tails == "right") CI <- 1 - confidence.interval
        
            # Resample control means
            y_p <- replicate(1000, sample.int(length(control), replace = TRUE))
            control_means <- Rfast::colmeans(matrix(control[y_p], ncol = ncol(y_p), byrow = T))

            a <- UPM.VaR(CI, 0, control_means)
            b <- mean(control_means)
            
            if(plot){
                if(tails == "both" | tails == "right"){
                    abline(v = max(a, b), col = "green", lwd = 2, lty = 3)
                    text(max(a, b), pos = 4, 0.5, paste0("<--- ", "ctl mu+ ", CI * 100, "%" ), col = "green")
                }
            }

            #Lower end of CDF confidence interval for control mean
            c <- LPM.VaR(CI, 0, control_means)
            d <- mean(control_means)

            if(plot){
                if(tails == "both" | tails == "left"){
                    abline(v = min(c, d), col = "blue", lwd = 2, lty = 3)
                    text(min(c, d), pos = 2, 0.5, paste0("ctl mu- ", paste0(CI * 100, "% --->")) , col = "blue")
                }

                par(original.par)
            }

            #Effect Size Lower Bound
            if(tails == "both") Lower.Bound.Effect <- mean(treatment) - max(a, b)
            
            if(tails == "left") Lower.Bound.Effect <- mean(treatment) - max(c, d)
            
            if(tails == "right") Lower.Bound.Effect <- mean(treatment) - max(a, b)
            


            #Effect Size Upper Bound
            if(tails == "both") Upper.Bound.Effect <- mean(treatment) - min(c, d)
            
            if(tails == "left") Upper.Bound.Effect <- mean(treatment) - min(c, d)
            
            if(tails == "right") Upper.Bound.Effect <- mean(treatment) - min(a, b)
            



        #Certainty Statistic and Effect Size Given Confidence Interval
        return(list("Control Mean" = mean(control),
                    "Treatment Mean" = mean(treatment),
                    "Grand Mean" = mean.of.means,
                    "Control CDF" = LPM_ratio.1,
                    "Treatment CDF" = LPM_ratio.2,
                    "Certainty" = min(1, NNS.ANOVA.rho * pop.adjustment),
                    "Lower Bound Effect" = Lower.Bound.Effect,
                    "Upper Bound Effect" = Upper.Bound.Effect))
  }
}
