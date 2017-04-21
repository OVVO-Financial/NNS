#' NNS Regression
#'
#' Generates a nonlinear regression based on partial moment quadrant means.
#'
#' @param x a vector, matrix or data frame of variables of numeric or factor data types.
#' @param y a numeric or factor vector with compatible dimsensions to \code{x}.
#' @param order integer; Controls the number of partial moment quadrant means.  Users are encouraged to try different \code{(order=...)} integer settings with \code{(noise.reduction="off")}.  \code{(order="max")} will force a limit condition perfect fit.
#' @param stn numeric [0,1]; Signal to noise parameter, sets the threshold of \code{(NNS.dep)} which reduces \code{("order")} when \code{(order=NULL)}.  Defaults to 0.99 to ensure high dependence for higher \code{("order")} and endpoint determination.  \code{(noise.reduction="off")} sets \code{(stn=0)} to allow for maximum fit.
#' @param type \code{NULL} (default).  To perform logistic regression, set to \code{(type = "LOGIT")}.  To perform a classification, set to \code{(type = "CLASS")}.
#' @param point.est a numeric or factor vector with compatible dimsensions to \code{x}.  Returns the fitted value \code{y.hat} for any value of \code{x}.
#' @param location Sets the legend location within the plot, per the \code{x} and \code{y} co-ordinates used in base graphics \link{legend}.
#' @param return.values logical; \code{TRUE} (default), set to \code{FALSE} in order to only display a regression plot and call values as needed.
#' @param plot  logical; \code{TRUE} (default) To plot regression.
#' @param plot.regions logical; \code{FALSE} (default).  Generates 3d regions associated with each regression point for multivariate regressions.  Note, adds significant time to routine.
#' @param residual.plot logical; \code{TRUE} (default) To plot \code{y.hat} and \code{Y}.
#' @param threshold  numeric [0,1]; \code{threshold=0} (default) Sets the correlation threshold for independent variables.
#' @param n.best integer; \code{NULL} (default) Sets the number of nearest regression points to use in weighting for multivariate regression at 2*(# of regressors).  \code{(n.best="all")} will select and weight all generated regression points.  Analogous to \code{k} in \code{k Nearest Neighbors} algorithm and different values are tested using cross-validation in \link{NNS.stack}.
#' @param noise.reduction the method of determing regression points options: ("mean","median","mode","off"); In low signal:noise situations,\code{(noise.reduction="mean")}  uses means for \link{NNS.dep} restricted partitions, \code{(noise.reduction="median")}  uses medians instead of means for \link{NNS.dep} restricted partitions, while \code{(noise.reduction="mode")}  uses modes instead of means for \link{NNS.dep} restricted partitions.  \code{(noise.reduction="off")}  allows for maximum possible fit with a specific \code{order}.
#' @param norm \code{NULL} (default) the method of normalization options: ("NNS","std"); Normalizes \code{x} between 0 and 1 for multivariate regression when set to \code{(norm="std")}, or normalizes \code{x} according to \link{NNS.norm} when set to \code{(norm="NNS")}.
#' @param dist options:("L1","L2") the method of distance calculation; Selects the distance calculation used. \code{dist="L2"} (default) selects the Euclidean distance and \code{(dist="L1")} seclects the Manhattan distance.
#' @param multivariate.call Internal parameter for multivariate regressions.
#' @return UNIVARIATE REGRESSION RETURNS THE FOLLOWING VALUES:
#'
#'      \code{"R2"} provides the goodness of fit;
#'
#'      \code{"MSE"} returns the MSE between \code{y} and \code{y.hat};
#'
#'      \code{"Prediction.Accuracy"} returns the correct rounded \code{"Point.est"} used in classifications versus the categorical \code{y};
#'
#'      \code{"derivative"} for the coefficient of the \code{x} and its applicable range;
#'
#'      \code{"Point"} returns the \code{x} point(s) being evaluated;
#'
#'      \code{"Point.est"} for the predicted value generated;
#'
#'      \code{"regression.points"} provides the points used in the regression equation for the given order of partitions;
#'
#'      \code{"partition"} returns the \code{"NNS.ID"} assigned to the observation and \code{y};
#'
#'      \code{"Fitted"} returns a vector containing only the fitted values, \code{y.hat};
#'
#'      \code{"Fitted.xy"} returns a \link{data.table} of \code{x},\code{y} and \code{y.hat};
#'
#'
#' MULTIVARIATE REGRESSION RETURNS THE FOLLOWING VALUES:
#'
#' \code{"R2"} provides the goodness of fit;
#'
#' \code{"equation"} returns the synthetic X* dimension reduction equation;
#'
#' \code{"rhs.partitions"} returns the partition points for each \code{x};
#'
#' \code{"RPM"} provides the Regression Point Matrix, the points for each \code{x} used in the regression equation for the given order of partitions;
#'
#' \code{"partition"} returns the \code{"NNS.ID"} assigned to the observation and \code{y};
#'
#' \code{"Point.est"} returns the predicted value generated;
#'
#' \code{"Fitted"} returns a vector containing only the fitted values, \code{y.hat};
#'
#' \code{"Fitted.xy"} reutnrs a \link{data.table} of \code{y} and fitted values.
#'
#' @note Please ensure \code{point.est} is of compatible dimensions to \code{x}, error message will ensue if not compatible.  Also, upon visual inspection of the data, if a highly periodic variable is observed set \code{(stn=0)} or \code{(order="max")} to ensure a proper fit.
#' @keywords nonlinear regression, classifier
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' set.seed(123)
#' x<-rnorm(100); y<-rnorm(100)
#' NNS.reg(x,y)
#'
#' ## Manual {order} selection
#' NNS.reg(x,y,order=2)
#'
#' ## Maximum {order} selection
#' NNS.reg(x,y,order="max")
#'
#' ## x-only paritioning (Univariate only)
#' NNS.reg(x,y,type="XONLY")
#'
#' ## Logistic Regression (Univariate only)
#' NNS.reg(x,y,type="LOGIT")
#'
#' ## For Multiple Regression:
#' x<-cbind(rnorm(100),rnorm(100),rnorm(100)); y<-rnorm(100)
#' NNS.reg(x,y,point.est=c(.25,.5,.75))
#'
#' ## For Multiple Regression based on Synthetic X* (Dimension Reduction):
#' x<-cbind(rnorm(100),rnorm(100),rnorm(100)); y<-rnorm(100)
#' NNS.reg(x,y,point.est=c(.25,.5,.75),type="CLASS")
#'
#' ## IRIS dataset example:
#' #Dimension Reduction:
#' NNS.reg(iris[,1:4],iris[,5],type="CLASS",order=5)
#' #Multiple Regression:
#' NNS.reg(iris[,1:4],iris[,5],order=2,noise.reduction="off")
#'
#' ## To call fitted values:
#' x<-rnorm(100); y<-rnorm(100)
#' NNS.reg(x,y)$Fitted
#'
#' ## To call partial derivative (univariate regression only):
#' x<-rnorm(100); y<-rnorm(100)
#' NNS.reg(x,y)$derivative
#'
#' @export


NNS.reg = function (x,y,
                    order=NULL,
                    stn=.99,
                    type = NULL,
                    point.est = NULL,
                    location = 'top',
                    return.values = TRUE,
                    plot = TRUE, plot.regions=FALSE,residual.plot=TRUE,
                    threshold = 0,
                    n.best=NULL,
                    noise.reduction='mean',
                    norm=NULL,
                    dist="L2",multivariate.call=FALSE){

  if(plot.regions==TRUE && order=='max'){stop('Please reduce the "order" or set "plot.regions = FALSE".')}

  original.columns = ncol(x)
  original.variable = x
  np = nrow(point.est)

  if(noise.reduction=='off'){stn=0}else{stn=stn}

  y=as.numeric(y)

  if(class(x)=='factor'){
    if(is.null(ncol(x))){
      x=as.numeric(x)
      if(!is.null(point.est)){point.est=as.numeric(point.est)}
    }else{
        x=apply(x,2,as.numeric)
        if(!is.null(point.est)){point.est=apply(point.est,2,as.numeric)}
        }

  }

  if(!is.null(ncol(original.variable))){
    if(ncol(original.variable)==1){
      x=original.variable
    } else {
      if(is.null(type)){
        if(!is.null(original.columns)){
        if(is.null(n.best)){n.best=2*original.columns} else {n.best=n.best}}
        if(is.null(original.columns)){if(is.null(n.best)){n.best=2} else {n.best=n.best}}

        return(NNS.M.reg(x,y,point.est=point.est,plot=plot,residual.plot=plot,order=order,n.best=n.best,type=type,location=location,noise.reduction=noise.reduction,norm = norm,dist=dist,stn = stn,return.values=return.values,plot.regions = plot.regions))
      } # Multivariate NULL type
          else{
              if(type=="CLASS"){
                  x= data.matrix(x)
                  y= as.numeric(y)

                  x.star.matrix = matrix(nrow=length(y))

                  x.star.dep = sapply(1:original.columns, function(i) abs(NNS.dep(unlist(x[,i]),unlist(y))$Dependence))
                  x.star.coef=sapply(1:original.columns, function(i) abs(NNS.cor(unlist(x[,i]),unlist(y))))

                  x.star.coef[x.star.coef<threshold]<- 0

                  x.star.matrix=t(t(original.variable) * x.star.coef)

                  synthetic.x.equation=(paste0("Synthetic Independent Variable X* = (",paste(format(x.star.coef[1:original.columns],digits = 4),paste("X",1:original.columns,sep = ''),sep='*',collapse = "  "),")/",sum(abs(x.star.coef)>0)))

                  #In case all IVs have 0 correlation to DV
                  if(all(x.star.matrix==0)){
                      x.star.matrix=x
                      x.star.coef[x.star.coef==0]<- 1
                  }

                  #Above threshold coefficients
                  coefs=abs(x.star.coef)
                  coefs=coefs[coefs>0]

                  if(!is.null(point.est)){
                      new.point.est=numeric()
                      if(is.null(np)){

                          new.point.est= sum(point.est*x.star.coef)/sum(abs(x.star.coef)>0)
                      }
                        else{
                          new.point.est=sapply(1:np,function(i) sum(point.est[i,]*x.star.coef)/sum(abs(x.star.coef)>0))

                        }
                    point.est=new.point.est
                  }

                  x = rowSums(x.star.matrix/sum(abs(x.star.coef)>0))

              } # Multivariate "CLASS" type
          } #Multivariate Not NULL type
      }

  } #Multivariate



  if(is.null(original.columns)){
    synthetic.x.equation=NULL

    dependence = NNS.dep(x,y,print.map = F)$Dependence

  } else {
    if(type=="CLASS") dependence=mean(x.star.dep)}

  if(is.null(order)){
      dep.reduced.order=floor(NNS.part(x,y,order='max')$order*dependence)}
      else {
      dep.reduced.order=order
      }

  if(dependence>stn ){
    if(is.null(type)){
      part.map = NNS.part(x,y,noise.reduction='off',order=dep.reduced.order)
      if(length(part.map$regression.points$x)==0){
        part.map=NNS.part(x,y,noise.reduction='off',order = min(nchar(part.map$dt$quadrant)),min.obs = 1)
      } else {part.map=part.map
      }
    } # NULL type

    if(!is.null(type)){
      part.map=NNS.part(x,y,type = "XONLY",noise.reduction='off',order = dep.reduced.order)
      if(length(part.map$regression.points$x)==0){
        part.map=NNS.part(x,y,noise.reduction='off',type="XONLY",order = min(nchar(part.map$dt$quadrant)),min.obs = 1)
      } else {part.map=part.map
      }
    } # type
  } # Dependence > stn

  if(dependence<=stn){
    if(is.null(type)){
      part.map = NNS.part(x,y,noise.reduction=noise.reduction, order=dep.reduced.order,type = "XONLY")
      if(length(part.map$regression.points$x)==0){
        part.map=NNS.part(x,y,type = "XONLY",noise.reduction=noise.reduction,order = min(nchar(part.map$dt$quadrant)),min.obs = 1)
      } else {part.map=part.map
      }
    } # NULL type


    if(!is.null(type)){
      part.map = NNS.part(x,y,type = "XONLY",noise.reduction=noise.reduction, order = dep.reduced.order)
      if(length(part.map$regression.points$x)==0){
        part.map=NNS.part(x,y,type = "XONLY",noise.reduction=noise.reduction,order = min(nchar(part.map$dt$quadrant)),min.obs = 1)
      } else {part.map=part.map
      }
    } # type
  } # Dependence < stn


  Regression.Coefficients = data.frame(matrix(ncol=3))
  colnames(Regression.Coefficients) = c('Coefficient','X Lower Range','X Upper Range')

  regression.points=part.map$regression.points[,.(x,y)]

  min.range = min(regression.points$x)
  max.range = max(regression.points$x)


  Dynamic.average.min = mean(y[x<=min.range])
  Dynamic.average.max = mean(y[x>=max.range])

  ###Endpoints
  if(length(x[x<min.range])>0){
    if(dependence<stn){
      x0 = Dynamic.average.min} else {
        x0 = unique(y[x==min(x)])} }  else {x0 = unique(y[x==min(x)])}

  if(length(x[x>max.range])>0){
    if(dependence<stn){x.max = Dynamic.average.max} else {x.max = unique(y[x==max(x)])}}  else { x.max = unique(y[x==max(x)])}


  regression.points=rbindlist(list(regression.points,list(max(x),mean(x.max))))
  regression.points=rbindlist(list(regression.points,list(min(x),mean(x0))))

  setkey(regression.points,x)

  ### Regression Equation

  regression.points =regression.points[complete.cases(regression.points*0)]

  if(multivariate.call==T){
    return(regression.points$x)
  }

  ### Consolidate possible duplicated points
  regression.points=regression.points[,y := mean(y),by=x]
  regression.points=unique(regression.points)

  rise = regression.points[,'rise' := y - shift(y)]
  run = regression.points[,'run' := x - shift(x)]


  Regression.Coefficients=regression.points[,.(rise,run)]


  Regression.Coefficients=Regression.Coefficients[complete.cases(Regression.Coefficients)]


  upper.x = regression.points[(2:.N),x]

  Regression.Coefficients=Regression.Coefficients[, `:=` ('Coefficient'=(rise/run),'X.Lower.Range'=regression.points[-.N,x],'X.Upper.Range'=upper.x)]

  Regression.Coefficients=Regression.Coefficients[,.(Coefficient,X.Lower.Range,X.Upper.Range)]


  Regression.Coefficients=unique(Regression.Coefficients)
  regression.points=regression.points[,.(x,y)]
  setkey(regression.points,'x')
  regression.points=unique(regression.points)

  ### Fitted Values
  p = length(regression.points[,x])

  fitted = numeric()
  fitted.order = numeric()

  point.order=numeric()
  point.est.y=numeric()

  if(is.na(Regression.Coefficients[1,Coefficient])){Regression.Coefficients[1,Coefficient:= Regression.Coefficients[2,Coefficient] ]}
  if(is.na(Regression.Coefficients[.N,Coefficient])){Regression.Coefficients[.N,Coefficient:= Regression.Coefficients[.N-1,Coefficient] ]}

  coef.interval=findInterval(x,Regression.Coefficients[,(X.Lower.Range)],left.open = FALSE)
  reg.interval=findInterval(x,regression.points[,x],left.open = FALSE)

  estimate=((x- regression.points[reg.interval,x])*Regression.Coefficients[coef.interval,Coefficient])+regression.points[reg.interval,y]

  if(!is.null(point.est)){
  coef.point.interval=findInterval(point.est,Regression.Coefficients[,(X.Lower.Range)],left.open = FALSE)
  reg.point.interval=findInterval(point.est,regression.points[,x],left.open=FALSE)
  coef.point.interval[coef.point.interval==0]<- 1
  reg.point.interval[reg.point.interval==0]<- 1
  point.est.y=((point.est - regression.points[reg.point.interval,x])*Regression.Coefficients[coef.point.interval,Coefficient])+regression.points[reg.point.interval,y]
  }


  fitted=(data.table(x=x,y=y,y.hat=estimate))

  Values = (cbind(x,Fitted=fitted[,y.hat],Actual=fitted[,y],Difference=fitted[,y.hat]-fitted[,y], Accuracy=abs(round(fitted[,y.hat])-fitted[,y])))


  MSE = fitted[,mean(y.hat-y)^2]
  y.fitted=fitted[,y.hat]

  Prediction.Accuracy=(length(y)-sum(abs(round(y.fitted)-(y))>0))/length(y)

  R2=  (sum((fitted[,y.hat]-mean(y))*(y-mean(y)))^2)/(sum((y-mean(y))^2)*sum((fitted[,y.hat]-mean(y))^2))

  R2.adj = R2

  ###Plotting and regression equation
  if(plot==TRUE){
    r2.leg=bquote(bold(R^2 == .(format(R2,digits=4))))
    xmin= min(c(point.est,x))
    xmax= max(c(point.est,x))
    ymin= min(c(point.est.y,y))
    ymax= max(c(point.est.y,y))

    if(is.null(order)){
      plot.order = dep.reduced.order} else {plot.order=order}

    plot(x,y,xlim=c(xmin,xmax),ylim=c(ymin,ymax),col='steelblue',main=paste(paste0("NNS Order = ", plot.order),sep="\n"),
         xlab = if(!is.null(original.columns))
            {if(original.columns>1){paste0("Synthetic X* ","(Segments = ",(p-1),')')}}else{paste0("X  ","(Segments = ",(p-1),")",sep='')},
         ylab="Y",mgp=c(2.5,0.5,0),
         cex.lab=2,cex.main=2)

    ### Plot Regression points and fitted values and legend
    points(na.omit(regression.points[,.(x,y)]),col='red',pch=15)
    lines(na.omit(regression.points[,.(x,y)]),col='red',lwd=2,lty = 2)


    if(!is.null(point.est)){
      points(point.est,point.est.y,
                                    col='green',pch=18,cex=1.5)

      legend(location, bty="n", y.intersp = 0.75,legend=r2.leg)}

    if(is.null(point.est)){
      legend(location, bty="n", y.intersp = 0.75,legend=r2.leg)
    }

    if(!is.null(point.est)){
        if(sum(point.est>max(x))>0){
          segments(point.est[point.est>max(x)],point.est.y[point.est>max(x)],regression.points[.N,x],regression.points[.N,y],col="green",lty=2)}
        if(sum(point.est<min(x))>0){
          segments(point.est[point.est<min(x)],point.est.y[point.est<min(x)],regression.points[1,x],regression.points[1,y],col="green",lty=2)
       }
    }
  }# plot TRUE bracket

  ### Return Values
  if(return.values == TRUE){
    return(list("R2"=R2, "MSE"=MSE, "Prediction.Accuracy"=Prediction.Accuracy,"equation"=synthetic.x.equation, "derivative"=Regression.Coefficients[],"Point"=point.est,"Point.est"=point.est.y,"regression.points"=regression.points[],"partition"=part.map$dt[,.(y,NNS.ID=quadrant)],"Fitted"=fitted[,.(y.hat)],"Fitted.xy"=fitted))
  } else {
    invisible(list("R2"=R2, "MSE"=MSE, "Prediction.Accuracy"=Prediction.Accuracy,"equation"=synthetic.x.equation, "derivative"=Regression.Coefficients[],"Point"=point.est,"Point.est"=point.est.y,"regression.points"=regression.points[],"partition"=part.map$dt[,.(y,NNS.ID=quadrant)],"Fitted"=fitted[,.(y.hat)],"Fitted.xy"=fitted))
  }

}
