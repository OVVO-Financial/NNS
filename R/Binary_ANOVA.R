
NNS.ANOVA.bin<- function(control,treatment,confidence.interval=NULL,tails=NULL,plot=TRUE){

        mean.of.means <- mean(c(mean(control),mean(treatment)))

  #Continuous CDF for each variable from Mean of Means
        LPM_ratio.1 <- LPM.ratio(1,mean.of.means,control)

        LPM_ratio.2 <- LPM.ratio(1,mean.of.means,treatment)


  #Continuous CDF Deviation from 0.5
        MAD.CDF<- mean(c(abs(LPM_ratio.1 - 0.5),abs(LPM_ratio.2 - 0.5)))


  #Certainty associated with samples
        NNS.ANOVA.rho <- (0.5 - MAD.CDF)^2/0.25


  #Graphs
        if(plot==TRUE){
        boxplot(list(control,treatment), las=2, names=c("Control","Treatment"),
              xlab= "Means", horizontal = TRUE, main= "NNS ANOVA and Effect Size",
              col=c("grey","white"),
              cex.axis= 0.75)

        #For ANOVA Visualization
        abline(v=mean.of.means,col="red",lwd=4)
        mtext("Grand Mean", side = 3,col = "red",at=mean.of.means)}

if(is.null(confidence.interval)){
    return(list("Control Mean" = mean(control),"Treatment Mean" = mean(treatment),"Grand Mean" = mean.of.means,"Control CDF" =LPM_ratio.1,"Treatment CDF" = LPM_ratio.2, "Certainty" = NNS.ANOVA.rho))}


if(!is.null(confidence.interval)){
        #Upper end of CDF confidence interval for control mean
        if(tails=="Both"){
            CI=confidence.interval+(1-confidence.interval)/2
        }
        if(tails=="Left" | tails=="Right") {
            CI=confidence.interval
        }


            a=UPM.VaR(CI,1,control)
            b=mean(control)
            if(plot==TRUE){
            if(tails=="Both"|tails=="Right"){
              abline(v=max(a,b),col="green",lwd=4,lty=3)
              text(max(a,b),pos=2,0.75,"mu+",col="green")
              text(max(a,b),pos=4,0.75,paste0((1-CI)*100,"% --->"),col="green")}
            }

        #Lower end of CDF confidence interval for control mean
            c=LPM.VaR(CI,1,control)
            d=mean(control)
            if(plot==TRUE){
            if(tails=="Both"|tails=="Left"){
              abline(v=min(c,d),col="blue",lwd=4,lty=3)
              text(min(c,d),pos=4,0.75,"mu-",col="blue")
              text(min(c,d),pos=2,0.75,paste0( "<--- ",(1-CI)*100,"%"),col='blue')}
            }

  #Effect Size Lower Bound
        if(tails=="Both"){
          Lower.Bound.Effect=min(mean(treatment)-max(a,b),0) }
        if(tails=="Left"){
          Lower.Bound.Effect=min(mean(treatment)-max(c,d),0) }
        if(tails=="Right"){
          Lower.Bound.Effect=min(mean(treatment)-max(a,b),0) }


  #Effect Size Upper Bound
        if(tails=="Both"){
          Upper.Bound.Effect=max(mean(treatment)-min(c,d),0) }
        if(tails=="Left"){
          Upper.Bound.Effect=max(mean(treatment)-min(c,d),0) }
        if(tails=="Right"){
          Upper.Bound.Effect=max(mean(treatment)-min(a,b),0) }

  #Certainty Statistic and Effect Size Given Confidence Interval
        return(list("Control Mean" = mean(control),"Treatment Mean" = mean(treatment),"Grand Mean" = mean.of.means,"Control CDF" =LPM_ratio.1,"Treatment CDF" = LPM_ratio.2, "Certainty" = NNS.ANOVA.rho,"Lower Bound Effect"=Lower.Bound.Effect,"Upper Bound Effect"=Upper.Bound.Effect))
}
}

