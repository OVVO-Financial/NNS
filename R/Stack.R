#' NNS stack
#'
#' Prediction model using the predictions of the NNS base models \link{NNS.reg} as features (i.e. meta-features) for the stacked model.
#'
#' @param IVs.train a vector, matrix or data frame of variables of numeric or factor data types.
#' @param IVs.test a vector, matrix or data frame of variables of numeric or factor data types.
#' @param DV.train a numeric or factor vector with compatible dimsensions to \code{(IVs.train)}.
#' @param CV.size numeric [0,1]; Sets the cross-validation size if \code{(IVs.test=NULL)}.  Defaults to 0.2 for a 20 percent random sampling of the training set.
#' @param weight options: ("MSE","Features") method for selecting model output weight; Set \code{(weight="MSE")} for optimum parameters and weighting based on each base model's \code{"MSE"}.  \code{(weight="Feautures")} uses a weighting based on the number of features present, whereby logistic \link{NNS.reg} receives higher relative weights for more regressors.  Defaults to \code{"MSE"}.
#' @param precision options: ("LOW","HIGH"); 2 settings offered: \code{"LOW"} (Default) ,and \code{"HIGH"}.  \code{"HIGH"} is the limit condition of every observation as a regression point and uses a \code{(norm="NNS")} while \code{(precision="LOW")} uses a \code{(norm="std")} in \link{NNS.reg}.  Errors/warnings can generally be reconciled with \code{(precision="LOW")}.
#' @param method numeric options: (1,2); Select the NNS method to include in stack.  \code{(method=1)} selects \link{NNS.reg}; \code{(method=2)} selects \link{NNS.reg} dimension reduction regression.  Defaults to \code{method=c(1,2)}, including both NNS regression methods in the stack.
#' @param threshold  numeric [0,1]; Sets the correlation threshold for independent variables in \link{NNS.reg}.  Defaults to \code{(threshold=0)}.
#' @param seed numeric; 123 (default) Sets seed for CV sampling.
#' @return Returns a vector of fitted values for the dependent variable test set for all models.
#' \itemize{
#' \item{\code{"NNS.reg.n.best"}} returns the optimum \code{"n.best"} paramater for the \link{NNS.reg} multivariate regression.  \code{"MSE.reg"} returns the MSE for the \link{NNS.reg} multivariate regression.
#' \item{\code{"NNS.dim.red.order"}} returns the optimum \code{"order"} from the \link{NNS.reg} dimension reduction regression.
#' \item{\code{"MSE.dim.red"}} returns the MSE for the \link{NNS.reg} dimension reduction regression.
#' \item{\code{"reg"}} returns \link{NNS.reg} output.
#' \item{\code{"dim.red"}} returns \link{NNS.reg} dimension reduction regression output.
#' \item{\code{"stack"}} returns the output of the stacked model.
#' }
#'
#' @keywords classifier
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. (2016) "Classification Using NNS Clustering Analyis"
#' \url{https://ssrn.com/abstract=2864711}
#' @note If character variables are used, transform them first to factors using \link{as.factor}, or \link{data.matrix} to ensure overall dataset is numeric.  A multifunction \link{sapply} can also be applied to the overall dataset: \code{data<- sapply(data,function(x){as.factor(x);as.numeric(x)})}.  Then run \code{NNS.stack} with transormed variables.
#' @examples
#'  ## Using 'iris' dataset where test set [IVs.test] is 'iris' rows 141:150.
#'  \dontrun{
#'  NNS.stack(iris[1:140,1:4],iris[1:140,5],IVs.test=iris[141:150,1:4])}
#'
#'  ## Using 'iris' dataset to determine [n.best] and [logistic.order] with no test set.
#'  \dontrun{
#'  NNS.stack(iris[,1:4],iris[,5])}
#'
#'  ## Selecting NNS.reg and dimension reduction techniques.
#'  \dontrun{
#'  NNS.stack(iris[1:140,1:4],iris[1:140,5],iris[141:150,1:4],method=c(1,2))}
#' @export

NNS.stack <- function(IVs.train,DV.train,IVs.test=NULL,CV.size=.2,weight="MSE",precision="LOW",method=c(1,2),threshold=0,seed=123){

  IVs.train<- apply(IVs.train,2,as.numeric)
  DV.train<- as.numeric(DV.train)

  n<- ncol(IVs.train)

  l=length(IVs.train[,1])
  #IV test provided...
  if(!is.null(IVs.test)){
  IVs.test<- data.matrix(IVs.test)
  }
  #No IV test provided...

  set.seed(seed*l)
  test.set=sample(1:l,as.integer(CV.size*l),replace = FALSE)

  CV.IVs.train<- IVs.train[c(-test.set),]

  CV.IVs.test<- IVs.train[c(test.set),]
  CV.DV.train<- DV.train[c(-test.set)]
  CV.DV.test<- DV.train[c(test.set)]

  IVs.train<- CV.IVs.train
  DV.train<- CV.DV.train

  ### NORMALIZATION OF VARIABLES and SELECTION OF ORDER:
  np=nrow(CV.IVs.test)
  points.norm=rbind(CV.IVs.test,CV.IVs.train)
  colnames(points.norm)=colnames(CV.IVs.test)
  if(precision=="LOW"){
    order=NULL
    CV.IVs.train=apply(CV.IVs.train,2,function(b) (b-min(b))/(max(b)-min(b)))
    CV.IVs.test=apply(points.norm,2,function(b) (b-min(b))/(max(b)-min(b)))[1:np,]
    }

  if(precision=="HIGH"){
    order="max"
    CV.IVs.train=NNS.norm(CV.IVs.train)
    CV.IVs.test=NNS.norm(points.norm)[1:np,]
  }


  if(1%in%method){

    nns.cv=sapply(1:(2*n),function(i) mean((NNS.reg(CV.IVs.train, CV.DV.train,point.est = CV.IVs.test, plot=F, n.best = i,order=order,noise.reduction = 'off')$Point.est-CV.DV.test)^2))

    nns.cv=c(nns.cv,mean((NNS.reg(CV.IVs.train, CV.DV.train,point.est = CV.IVs.test, plot=F, n.best = 'all',order=order,noise.reduction = 'off')$Point.est-CV.DV.test)^2))

    k=which.min(nns.cv)
    if(k==length(nns.cv)){k='all'}else{k=k}

    nns.1=NNS.reg(IVs.train, DV.train,point.est = IVs.test, plot=F, n.best = k,order=order,noise.reduction = 'off')$Point.est
    nns.cv.mse=min(nns.cv)
  } else {
    nns.1=0
    nns.cv=1e-10
    nns.cv.mse='N/A'}

  # Logistic Regression Output
  if(2%in%method){

    nns.ord=sapply(1:log(l,2),function(i) mean((NNS.reg(CV.IVs.train, CV.DV.train,point.est = CV.IVs.test,plot=F, order = i, type = "CLASS",threshold = threshold,noise.reduction = 'off')$Point.est-CV.DV.test)^2))

    nns.2=NNS.reg(IVs.train, DV.train,point.est = IVs.test, type = "CLASS",noise.reduction = 'off',plot=F,order=which.min(nns.ord),threshold = threshold)$Point.est
    nns.ord.mse=min(nns.ord)
  } else {nns.2=0
  nns.ord=1e-10
  nns.ord.mse='N/A'}



  ### Weights for combining NNS techniques
  if(weight=="Features"){weights=c(n,n-1)}

  nns.cv[nns.cv==0]<- 1e-10
  nns.ord[nns.ord==0]<- 1e-10


  if(weight=="MSE"){weights=c(max(1e-10,1/min(na.omit(nns.cv))),max(1e-10,1/min(na.omit(nns.ord))))}


  weights=pmax(weights,c(0,0))
  weights[!(c(1,2)%in%method)]<- 0
  weights=weights/sum(weights)


  estimates = (weights[1]*nns.1+weights[2]*nns.2)

  return(list(NNS.reg.n.best=k,MSE.reg=nns.cv.mse,NNS.dim.red.order=which.min(nns.ord),MSE.dim.red=nns.ord.mse,reg=nns.1,dim.red=nns.2,stack=estimates))

}
