NNS.ANOVA.bin<- function(control, treatment,
                         mean.of.means = NULL,
                         upper.target = NULL,
                         lower.target = NULL,
                         confidence.interval = NULL, tails = NULL, plot = TRUE){

  if(is.null(upper.target) && is.null(lower.target)){
        mean.of.means <- mean(c(mean(control), mean(treatment)))
        upper.target <- mean(c(UPM.VaR(.75, 1, control), UPM.VaR(.75, 1, treatment)))
        lower.target <- mean(c(LPM.VaR(.75, 1, control), LPM.VaR(.75, 1, treatment)))
  }


  #Continuous CDF for each variable from Mean of Means
        LPM_ratio.1 <- LPM.ratio(1, mean.of.means, control)
        LPM_ratio.2 <- LPM.ratio(1, mean.of.means, treatment)

        Upper_ratio.1 <- UPM.ratio(1, upper.target, control)
        Upper_ratio.2 <- UPM.ratio(1, upper.target, treatment)
        Upper_ratio <- sqrt((Upper_ratio.1* Upper_ratio.2))

        Lower_ratio.1 <- LPM.ratio(1, lower.target, control)
        Lower_ratio.2 <- LPM.ratio(1, lower.target, treatment)
        Lower_ratio <- sqrt((Lower_ratio.1* Lower_ratio.2))


  #Continuous CDF Deviation from 0.5
        MAD.CDF <- mean(c(abs(LPM_ratio.1 - 0.5), abs(LPM_ratio.2 - 0.5)))
        upper.CDF <- mean(c(Upper_ratio.1, Upper_ratio.2))
        lower.CDF <- mean(c(Lower_ratio.1, Lower_ratio.2))


  #Certainty associated with samples
        NNS.ANOVA.rho <- sum(c( ((.5- MAD.CDF)^2) / .25 ,
                                ifelse(Upper_ratio==0,0,.5 * ( abs(Upper_ratio)/ upper.CDF)),
                                ifelse(Lower_ratio==0,0,.5 * ( abs(Lower_ratio)/ lower.CDF))))/2

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
