require(NNS)
require(Deriv)
require(xtable)
require(np)


PartialDrep=function(x,y,f){
  out=matrix(nrow=1,ncol=8)
  
  Analytic= Deriv(f)(median(x))
  PDyx= dy.dx(x,y,deriv.method = "NNS")
  diff=abs(Analytic-PDyx)
  PDyx2=dy.dx(x,y,deriv.method = "FS")
  diff2=abs(Analytic-PDyx2)
  np.d=gradients(npreg(tydat=y,txdat=x,exdat=median(x),gradients=TRUE))
  np.diff=abs(Analytic-np.d)
  out[1,1]=Analytic
  out[1,2]=PDyx
  out[1,3]=diff
  out[1,4]=PDyx2
  out[1,5]=diff2
  out[1,6]=np.d
  out[1,7]=np.diff
  out[1,8]=NNS.dep(x,y,print.map = FALSE)[2]
  
  return(out)}

## Simulation parameters
for(ii in c(160,80,40,20,10,1)){
n=300           # Sample Size
reps=1000        # Iterations for each experiment
m=ii            # Noise setting
print(m)
options(digits = 3)
myexpm=c(1:6) 
nexpm=length(myexpm)
out=matrix(nrow=reps,ncol=8)
ou1.1=matrix(ncol=8)
ou1.2=matrix(ncol=8)
ou1.3=matrix(ncol=8)
ou1.4=matrix(ncol=8)
ou1.5=matrix(ncol=8)
ou1.6=matrix(ncol=8)


for (expm in 1:nexpm){
  for(i in 1:reps){set.seed(234*i)
    
    if (expm==1){x1=1:n
    y1=3+4*x1+rnorm(n,mean=runif(1,-1,1),sd=1/m)+runif(n,-1,1)
    f= function(x) 3+4*x

    ou1.1=rbind(ou1.1,PartialDrep(x1,y1,f))}
    
    if (expm==2){x1=1:n
    y1=3+4*x1-3*x1^2+rnorm(n,mean=runif(1,-1,1),sd=1/m)+runif(n,-1,1)
    f= function(x) 3+4*x-3*x^2

    ou1.2=rbind(ou1.2,PartialDrep(x1,y1,f))}
    
    if (expm==3){x1=runif(n)
    y1=3^(1+2*x1+3*x1^2)+rnorm(n,mean=runif(1,-1,1),sd=1/m)+runif(n,-1,1)
    f= function(x) 3^(1+2*x+3*x^2)

    ou1.3=rbind(ou1.3,PartialDrep(x1,y1,f))}
    
    if (expm==4){x1=runif(n)
    y1=exp(1)^(.1+.2*x1)+rnorm(n,mean=runif(1,-1,1),sd=1/m)+runif(n,-1,1)
    f= function(x) exp(1)^(.1+.2*x)

    ou1.4=rbind(ou1.4,PartialDrep(x1,y1,f))}
    
    if (expm==5){x1=runif(n)
    y1=(1/(exp(1)^(1+2*x1)))+rnorm(n,mean=runif(1,-1,1),sd=1/m)+runif(n,-1,1)
    f= function(x) (1/(exp(1)^(1+2*x)))

    ou1.5=rbind(ou1.5,PartialDrep(x1,y1,f))}
    
    if (expm==6){x1=rnorm(n)
    y1=exp(1)^(.1-.2*x1)+rnorm(n,mean=runif(1,-1,1),sd=1/m)+runif(n,-1,1)
    f= function(x) exp(1)^(.1-.2*x)

    ou1.6=rbind(ou1.6,PartialDrep(x1,y1,f))}
  }
} 

#First row is na, needs to be eliminated...
ou2.1=apply(apply((ou1.1[-1,]),2,as.numeric),2,mean)
ou2.2=apply(apply((ou1.2[-1,]),2,as.numeric),2,mean)
ou2.3=apply(apply((ou1.3[-1,]),2,as.numeric),2,mean)
ou2.4=apply(apply((ou1.4[-1,]),2,as.numeric),2,mean)
ou2.5=apply(apply((ou1.5[-1,]),2,as.numeric),2,mean)
ou2.6=apply(apply((ou1.6[-1,]),2,as.numeric),2,mean)


out=rbind(ou2.1,ou2.2,ou2.3,ou2.4,ou2.5,ou2.6)

#MAPE columns
out[,3]=abs(out[,1]-out[,2])/abs(out[,1])
out[,5]=abs(out[,1]-out[,4])/abs(out[,1])
out[,7]=abs(out[,1]-out[,6])/abs(out[,1])

colnames(out)=c("Analytic dy/dx","NNS dy/dx","NNS MAPE","NNS (FS) dy/dx","NNS (FS) MAPE","np dy/dx","np MAPE","Dependence")
rownames(out)=c(paste("DT",1:2,sep=""),paste("ST",1:4,sep=""))

print(xtable(out,caption= paste0("Experiments from section 4.0.2 where $m=",m,"$"),
             caption.placement='top',
             table.placement="H"
             ,digits = xdigits(out)),type='latex')
}
