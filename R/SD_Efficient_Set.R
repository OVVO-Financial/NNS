#' NNS SD Efficient Set
#'
#' Determines the set of stochastic dominant variables for various degrees.
#' @param x a numeric matrix or data frame.
#' @param degree numeric options: (1,2,3); Degree of stochastic dominance test from (1,2 or 3).
#' @return Returns set of stochastic dominant variable names.
#' @keywords stochastic dominance
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2016) "LPM Density Functions for the Computation of the SD Efficient Set." Journal of Mathematical Finance, 6, 105-126. \url{http://www.scirp.org/Journal/PaperInformation.aspx?PaperID=63817}.
#' @examples
#' set.seed(123)
#' x<-rnorm(100); y<-rnorm(100); z<-rnorm(100)
#' x<-data.frame(x,y,z)
#' NNS.SD.Efficient.Set(x,1)
#' @export



NNS.SD.Efficient.Set <- function(x,degree) {
  n <- ncol(x)
  max_target <- max(x)
  LPM_order<- numeric(0)
  Dominated_set<- numeric(0)
  current_base<- numeric(0)


  LPM_order=sapply(1:n,function(i) LPM(1,max(x),x[,i]))

  final_ranked <- x[,order(LPM_order)]

  current_base<-1

  if(degree==1){
  for (i in 1:(n-1)) {

      base<- final_ranked[,current_base[length(current_base)]]

      challenger <- final_ranked[,i+1]

      if (NNS.FSD.uni(base,challenger)==1){ current_base[i]<- current_base[length(current_base)]}


      if (NNS.FSD.uni(base,challenger)==0){

        for (j in current_base){
          base<- final_ranked[,j]
          if (NNS.FSD.uni(base,challenger)==0){ next }
          else
            {Dominated_set[i] <- i+1  }
          }

       current_base[i]<- i+1}

      else {Dominated_set[i]<- i+1 }
  }
}


 if(degree==2){
   for (i in 1:(n-1)) {

     base<- final_ranked[,current_base[length(current_base)]]

     challenger <- final_ranked[,i+1]

     if (NNS.SSD.uni(base,challenger)==1){ current_base[i]<- current_base[length(current_base)]}


     if (NNS.SSD.uni(base,challenger)==0){

       for (j in current_base){
         base<- final_ranked[,j]
         if (NNS.SSD.uni(base,challenger)==0){ next }
         else
         {Dominated_set[i] <- i+1  }
       }

       current_base[i]<- i+1}

     else {Dominated_set[i]<- i+1 }
   }
 }


 if(degree==3){
   for (i in 1:(n-1)) {

     base<- final_ranked[,current_base[length(current_base)]]

     challenger <- final_ranked[,i+1]

     if (NNS.TSD.uni(base,challenger)==1){ current_base[i]<- current_base[length(current_base)]}


     if (NNS.TSD.uni(base,challenger)==0){

       for (j in current_base){
         base<- final_ranked[,j]
         if (NNS.TSD.uni(base,challenger)==0){ next }
         else
         {Dominated_set[i] <- i+1  }
       }

       current_base[i]<- i+1}

     else {Dominated_set[i]<- i+1 }
   }
 }


  if(length(Dominated_set)>0){
    SD_x  = data.frame(final_ranked[-na.omit(Dominated_set)])
    return(SD_x)
    }

    else {print(colnames(final_ranked))}


}
