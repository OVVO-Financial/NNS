NNS.ANOVA.bin<- function(control, treatment,
                         mean.of.means = NULL,
                         upper.25.target = NULL,
                         lower.25.target = NULL,
                         upper.125.target = NULL,
                         lower.125.target = NULL,
                         confidence.interval = NULL, tails = NULL, plot = TRUE){

  if(is.null(upper.25.target) && is.null(lower.25.target)){
        mean.of.means <- mean(c(mean(control), mean(treatment)))
        upper.25.target <- mean(c(UPM.VaR(.75, 1, control), UPM.VaR(.75, 1, treatment)))
        lower.25.target <- mean(c(LPM.VaR(.75, 1, control), LPM.VaR(.75, 1, treatment)))
        upper.125.target <- mean(c(UPM.VaR(.875, 1, control), UPM.VaR(.875, 1, treatment)))
        lower.125.target <- mean(c(LPM.VaR(.875, 1, control), LPM.VaR(.875, 1, treatment)))
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
        MAD.CDF <- mean(c(abs(LPM_ratio.1 - 0.5), abs(LPM_ratio.2 - 0.5)))
        upper.25.CDF <- mean(c(Upper_25_ratio.1, Upper_25_ratio.2))
        lower.25.CDF <- mean(c(Lower_25_ratio.1, Lower_25_ratio.2))
        upper.125.CDF <- mean(c(Upper_125_ratio.1, Upper_125_ratio.2))
        lower.125.CDF <- mean(c(Lower_125_ratio.1, Lower_125_ratio.2))


  #Certainty associated with samples
        NNS.ANOVA.rho <- sum(c( ((.5- MAD.CDF)^2) / .25 ,
                                ifelse(Upper_25_ratio==0, 0, .5 * ( abs(Upper_25_ratio)/ upper.25.CDF)),
                                ifelse(Lower_25_ratio==0, 0, .5 * ( abs(Lower_25_ratio)/ lower.25.CDF)),
                                ifelse(Upper_125_ratio==0, 0, .25 * ( abs(Upper_125_ratio)/ upper.125.CDF)),
                                ifelse(Lower_125_ratio==0, 0, .25 * ( abs(Lower_125_ratio)/ lower.125.CDF)))
                             ) / 2.25

        pop.adjustment <- ((length(control) + length(treatment) - 2) / (length(control)  + length(treatment) )) ^ 2

  #Graphs

        if(plot){
            original.par <- par(no.readonly = TRUE)

            boxplot(list(control, treatment), las = 2, names = c("Control", "Treatment"), xlab = "Means", horizontal = TRUE, main = "NNS ANOVA and Effect Size", col = c("grey", "white"), cex.axis = 0.75)

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
                "Certainty" = min(1, NNS.ANOVA.rho * pop.adjustment)))}


    if(!is.null(confidence.interval)){
        #Upper end of CDF confidence interval for control mean
        if(tails == "both"){
            CI <- confidence.interval+(1-confidence.interval)/2
        }
        if(tails == "left" | tails == "right") {
            CI <- confidence.interval
        }


            a <- UPM.VaR(CI, 1, control)
            b <- mean(control)
            if(plot){
                if(tails == "both" | tails == "right"){
                    abline(v = max(a, b), col = "green", lwd = 4, lty = 3)
                    text(max(a, b), pos = 2, 0.75, "mu+", col = "green")
                    text(max(a, b), pos = 4, 0.75, paste0((1 - CI) * 100, "% --->"), col = "green")}
            }

            #Lower end of CDF confidence interval for control mean
            c <- LPM.VaR(CI, 1, control)
            d <- mean(control)

            if(plot){
                if(tails == "both" | tails == "left"){
                    abline(v = min(c, d), col = "blue", lwd = 4, lty = 3)
                    text(min(c, d), pos = 4, 0.75, "mu-", col = "blue")
                    text(min(c, d), pos=2, 0.75, paste0( "<--- ", (1 - CI) * 100, "%"), col = 'blue')}

                par(original.par)
            }

            #Effect Size Lower Bound
            if(tails == "both"){
                Lower.Bound.Effect <- min(mean(treatment) - max(a, b), 0)
            }
            if(tails == "left"){
                Lower.Bound.Effect <- min(mean(treatment) - max(c, d), 0)
            }
            if(tails == "right"){
                Lower.Bound.Effect <- min(mean(treatment) - max(a, b), 0)
            }


  #Effect Size Upper Bound
            if(tails == "both"){
                Upper.Bound.Effect <- max(mean(treatment) - min(c, d), 0)
            }
            if(tails == "left"){
                Upper.Bound.Effect <- max(mean(treatment) - min(c, d), 0)
            }
            if(tails == "right"){
                Upper.Bound.Effect <- max(mean(treatment) - min(a, b), 0)
            }



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
