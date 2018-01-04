NNS.M.reg <- function (X_n,Y,order=NULL,stn=0.99,n.best=1,type=NULL,point.est=NULL, plot=FALSE,residual.plot=TRUE,location=NULL,noise.reduction='mean',norm=NULL,dist="L2",return.values=FALSE,plot.regions=FALSE){

  if(is.null(ncol(X_n))){X_n=t(t(X_n))}
  n=ncol(X_n)

  if(is.null(names(Y))){y.label="Y"}else{y.label=names(Y)}

  ### Mode calculation
  mode=function(x) {
    if(length(x)>1){
      d <- density(x)
      d$x[which.max(d$y)]
    }else{x}
  }

  np=nrow(point.est)
  if(is.null(np)&!is.null(point.est)){point.est=t(point.est)
  }else{point.est=point.est}


  if(!is.null(point.est)){
    if(ncol(point.est)!=n){
      stop("Please ensure 'point.est' is of compatible dimensions to 'x'")
    }
  }

  ### For Multiple regressions
  ###  Turn each column into numeric values
  original.IVs = data.matrix(X_n)
  original.DV=as.numeric(Y)

  if(!is.null(norm)){
    if(norm=='std'){
      original.IVs=apply(original.IVs,2,function(b) (b-min(b))/(max(b)-min(b)))}else{
        original.IVs=NNS.norm(original.IVs)}
    if(!is.null(point.est)){
      point.B=rbind(point.est,original.IVs)
      colnames(point.B)=colnames(point.est)
      if(norm=='std'){
        point.est=apply(point.B,2,function(c) (c-min(c))/(max(c)-min(c)))[1:np,]}else{
          point.est=NNS.norm(point.B)[1:np,]
        }
    }}

  original.matrix=data.frame(original.IVs,original.DV)


  reg.points=list()
  sections = list()

  ###  Regression Point Matrix
  if(is.numeric(order)|is.null(order)){

    reg.points=apply(original.IVs,2,function(b) NNS.reg(b,original.DV,order=order,type=type,noise.reduction=noise.reduction,plot = F,multivariate.call = TRUE)$x)

    if(all(sapply(reg.points, length) == length(reg.points[[1]]))==FALSE){
      reg.points.matrix=do.call('cbind',lapply(reg.points, `length<-`, max(lengths(reg.points))))}
    else {reg.points.matrix=reg.points}
  } else {
    reg.points.matrix=original.IVs
  }


  ### If regression points are error (not likely)...
  if(length(reg.points.matrix[,1])==0){
    for(i in 1:n){
      part.map=NNS.part(original.IVs[,i],original.DV,order=order,type=type,noise.reduction=noise.reduction)
      dep=NNS.dep(original.IVs[,i],original.DV,order=3)$Dependence
      if(dep>stn){
        reg.points[[i]] = NNS.part(original.IVs[,i],original.DV,order=round(dep*max(nchar(part.map$df$quadrant))),type=type,noise.reduction='off',max.obs = 1)$regression.points$x}
      else{reg.points[[i]] = NNS.part(original.IVs[,i],original.DV,order=round(dep*max(nchar(part.map$df$quadrant))),noise.reduction=noise.reduction,type="XONLY",max.obs = 1)$regression.points$x}
    }
    reg.points.matrix=do.call('cbind',lapply(reg.points, `length<-`, max(lengths(reg.points))))
  }

  if(is.null(colnames(original.IVs))){
    colnames.list=list()
    for(i in 1:n){
      colnames.list[i]=paste0("X",i)
    }
    colnames(reg.points.matrix)=as.character(colnames.list)}

  reg.points.matrix = unique(reg.points.matrix)
  ### Find intervals in regression points for each variable, use left.open T and F for endpoints.
  NNS.ID = list()

  for(j in 1:n){
    sorted.reg.points=sort(reg.points.matrix[,j])
    sorted.reg.points=sorted.reg.points[!is.na(sorted.reg.points)]
    NNS.ID[[j]]= findInterval(original.IVs[,j],sorted.reg.points,left.open = FALSE)
  }

  NNS.ID = matrix(unlist(NNS.ID),nrow = length(Y),ncol = n)

  ### Create unique identifier of each observation's interval
  NNS.ID = apply(NNS.ID, 1 , paste , collapse = "." )


  ### Match y to unique identifier
  obs=c(1:length(Y))

  mean.by.id.matrix = data.table(original.IVs,original.DV,NNS.ID,obs)
  setkey(mean.by.id.matrix,'NNS.ID','obs')

  if(noise.reduction=='mean'|noise.reduction=='off'){
    mean.by.id.matrix=mean.by.id.matrix[,y.hat := mean(original.DV),by='NNS.ID']}
  if(noise.reduction=='median'){
    mean.by.id.matrix=mean.by.id.matrix[,y.hat := median(original.DV),by='NNS.ID']
  }
  if(noise.reduction=='mode'){
    mean.by.id.matrix=mean.by.id.matrix[,y.hat := mode(original.DV),by='NNS.ID']
  }

  y.identifier=mean.by.id.matrix[,NNS.ID]

  ###Order y.hat to order of original Y
  resid.plot=mean.by.id.matrix[]
  setkey(resid.plot,'obs')
  y.hat=mean.by.id.matrix[,.(y.hat)]

  fitted.matrix = data.table(original.IVs,y=original.DV,y.hat,mean.by.id.matrix[,.(NNS.ID)])

  setkey(mean.by.id.matrix,'NNS.ID')
  REGRESSION.POINT.MATRIX=mean.by.id.matrix[,obs := NULL]



  if(noise.reduction=='mean'|noise.reduction=='off'){
    REGRESSION.POINT.MATRIX=REGRESSION.POINT.MATRIX[, lapply(.SD, mean), by=NNS.ID]}
  if(noise.reduction=='median'){
    REGRESSION.POINT.MATRIX=REGRESSION.POINT.MATRIX[, lapply(.SD, median), by=NNS.ID]
  }
  if(noise.reduction=='mode'){
    REGRESSION.POINT.MATRIX=REGRESSION.POINT.MATRIX[, lapply(.SD, mode), by=NNS.ID]
  }

  REGRESSION.POINT.MATRIX=REGRESSION.POINT.MATRIX[,NNS.ID := NULL]
  REGRESSION.POINT.MATRIX=REGRESSION.POINT.MATRIX[,original.DV := NULL]


  if(!is.numeric(n.best)){n.best=REGRESSION.POINT.MATRIX[,.N]} else {n.best=n.best}

  ### DISTANCES
  ### Calculate distance from each point in REGRESSION.POINT.MATRIX
  if(!is.null(point.est)){
    distance<- function(dist.est){

      for(j in 1:n){set(REGRESSION.POINT.MATRIX,j=j,value=REGRESSION.POINT.MATRIX[[j]]-as.numeric(dist.est)[j])}


      if(dist=="L1"){
        row.sums=as.numeric(rowSums(abs(REGRESSION.POINT.MATRIX[,.SD,.SDcols=c(1:n)])))}
      else{row.sums=as.numeric(rowSums(REGRESSION.POINT.MATRIX[,.SD,.SDcols=c(1:n)]^2))
      }

      for(j in 1:n){set(REGRESSION.POINT.MATRIX,j=j,value=REGRESSION.POINT.MATRIX[[j]]+as.numeric(dist.est)[j])}

      row.sums[row.sums==0]<- 1e-10
      total.row.sums = sum(1/row.sums)
      weights = (1/row.sums)/total.row.sums


      highest=rev(order(weights))[1:min(n.best,length(weights))]

      weights[-highest]<-0
      weights.sum=sum(weights)

      weights=weights/weights.sum
      single.estimate = sum(weights*REGRESSION.POINT.MATRIX$y.hat)

      return(single.estimate)
    }


    ### Point estimates

    predict.fit=numeric()
    predict.fit.iter=numeric()
    if(is.null(np)){
      predict.fit = distance(dist.est = point.est)
    }
    if(!is.null(np)){
      predict.fit.iter=apply(point.est,1,function(p) distance(dist.est = as.vector(p) ))

      predict.fit=as.vector(predict.fit.iter)
    }

  } else {predict.fit=NULL} #is.null point.est


  R2=  (sum((y.hat-mean(original.DV))*(original.DV-mean(original.DV)))^2)/(sum((original.DV-mean(original.DV))^2)*sum((y.hat-mean(original.DV))^2))


  ### 3d plot
  if(plot && n==2){
    region.1=mean.by.id.matrix[[1]]
    region.2=mean.by.id.matrix[[2]]
    region.3=mean.by.id.matrix[,y.hat]




    plot3d(x=original.IVs[,1],y=original.IVs[,2],z=original.DV,box=F,size = 3,col='steelblue',xlab=colnames(reg.points.matrix)[1], ylab=colnames(reg.points.matrix)[2], zlab=y.label )


    if(plot.regions){
      region.matrix = data.table(original.IVs,original.DV,NNS.ID)
      region.matrix[,`:=`(min.x1=min(.SD), max.x1=max(.SD)), by=NNS.ID,.SDcols=1]
      region.matrix[,`:=`(min.x2=min(.SD), max.x2=max(.SD)), by=NNS.ID,.SDcols=2]
      if(noise.reduction=='off'|noise.reduction=='mean'){
        region.matrix[,`:=`(y.hat=mean(original.DV)),by=NNS.ID]}
      if(noise.reduction=='median'){
        region.matrix[,`:=`(y.hat=median(original.DV)),by=NNS.ID]
      }
      if(noise.reduction=='mode'){
        region.matrix[,`:=`(y.hat=mode(original.DV)),by=NNS.ID]
      }
      setkey(region.matrix,NNS.ID,min.x1,max.x1,min.x2,max.x2)
      region.matrix[ ,{
        quads3d(x=.(min.x1[1],min.x1[1],max.x1[1],max.x1[1]),y=.(min.x2[1],max.x2[1],max.x2[1],min.x2[1]),z=.(y.hat[1],y.hat[1],y.hat[1],y.hat[1]),col='pink',alpha=1)
        if(identical(min.x1[1],max.x1[1])|identical(min.x2[1],max.x2[1])){
          segments3d(x=.(min.x1[1],max.x1[1]),y=.(min.x2[1],max.x2[1]),z=.(y.hat[1],y.hat[1]),col='pink',alpha=1)}
      }
      , by = NNS.ID]
    }#plot.regions = T


    points3d(x=as.numeric(unlist(REGRESSION.POINT.MATRIX[,.SD,.SDcols=1])),y=as.numeric(unlist(REGRESSION.POINT.MATRIX[,.SD,.SDcols=2])),z=as.numeric(unlist(REGRESSION.POINT.MATRIX[,.SD,.SDcols=3])),col='red',size=5)
    if(!is.null(point.est)){
      if(is.null(np)){
        points3d(x=point.est[1],y=point.est[2],z=predict.fit,col='green',size=5)
      } else {
        points3d(x=point.est[,1],y=point.est[,2],z=predict.fit,col='green',size=5)
      }
    }
  }

  ### Residual plot
  if(residual.plot){
    resids=cbind(original.DV,y.hat)
    r2.leg=bquote(bold(R^2 == .(format(R2,digits=4))))
    matplot(resids,type = 'l',xlab="Index",ylab=expression(paste("y (black)   ",hat(y), " (red)")),cex.lab=1.5,mgp=c(2,.5,0))

    title(main = paste0("NNS Order = ",order),cex.main=2)
    legend(location,legend =r2.leg,bty = 'n')
  }


  rhs.partitions=data.table(reg.points.matrix)
  ### Return Values
  if(return.values){

    return(list(R2=R2,rhs.partitions=rhs.partitions, RPM=REGRESSION.POINT.MATRIX[] ,Point.est=predict.fit,Fitted=fitted.matrix[,.(y.hat)],Fitted.xy=fitted.matrix[]))}

  else{
    invisible(list(R2=R2,rhs.partitions=rhs.partitions, RPM=REGRESSION.POINT.MATRIX[],Point.est=predict.fit,Fitted=fitted.matrix[,.(y.hat)],Fitted.xy=fitted.matrix[]))}

}
