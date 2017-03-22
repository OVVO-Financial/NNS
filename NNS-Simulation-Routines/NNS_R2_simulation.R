require(NNS)
require(np)
require(xtable)

R2rep=function(x,y,reps=1){
  out=matrix(nrow=reps,ncol=4)
  
  for (i in 1:reps){
    R2.NNS=NNS.reg(x,y,plot=F)$R2
    R2.np= npreg(y~x)$R2
    NNS.np.diff=R2.NNS-R2.np
    out[i,1]=R2.NNS
    out[i,2]=R2.np
    out[i,3]=NNS.np.diff
    out[i,4]=NNS.dep(x,y,print.map = FALSE)[2]}
  
  return(out)}



## Simulation parameters
for(ii in c(30,50,100,200,300)){
n=ii           # Sample Size
reps=1000      # Iterations for each experiment

print(ii)

options(digits = 6)


myexpm=c(1:6) 
nexpm=length(myexpm)

out=matrix(nrow=reps,ncol=6)
ou1.1=matrix(ncol=4)
ou1.2=matrix(ncol=4)
ou1.3=matrix(ncol=4)
ou1.4=matrix(ncol=4)
ou1.5=matrix(ncol=4)
ou1.6=matrix(ncol=4)

for (expm in 1:nexpm){
  for(i in 1:reps){set.seed(234*i)

    if (expm==1){x1=1:n
    y1=3+4*x1+rnorm(n,mean=runif(1,-1,1),sd=1)+runif(n,-1,1)
    ou1.1=rbind(ou1.1,R2rep(x1,y1,reps=1))}
    if (expm==2){x1=1:n
    y1=3+4*x1-3*x1^2+rnorm(n,mean=runif(1,-1,1),sd=1)+runif(n,-1,1)
    ou1.2=rbind(ou1.2,R2rep(x1,y1,reps=1))}
    if (expm==3){x1=runif(n)
    y1=3^(1+2*x1 +3*x1^2)+rnorm(n,mean=runif(1,-1,1),sd=1)+runif(n,-1,1)
    ou1.3=rbind(ou1.3,R2rep(x1,y1,reps=1))}
    if (expm==4){x1=runif(n)
    y1=exp(.1+.2*x1)+rnorm(n,mean=runif(1,-1,1),sd=1)+runif(n,-1,1)
    ou1.4=rbind(ou1.4,R2rep(x1,y1,reps=1))}
    if (expm==5){x1=runif(n)
    y1=1 /exp(1+2*x1)+rnorm(n,mean=runif(1,-1,1),sd=1)+runif(n,-1,1)
    ou1.5=rbind(ou1.5,R2rep(x1,y1,reps=1))}
    if (expm==6){x1=rnorm(n)
    y1=1/ (1+exp(.1-.2*x1))+rnorm(n,mean=runif(1,-1,1),sd=1)+runif(n,-1,1)
    ou1.6=rbind(ou1.6,R2rep(x1,y1,reps=1))}
  }
} 

ou2.1=apply(apply((ou1.1[-1,]),2,as.numeric),2,mean)
ou2.2=apply(apply((ou1.2[-1,]),2,as.numeric),2,mean)
ou2.3=apply(apply((ou1.3[-1,]),2,as.numeric),2,mean)
ou2.4=apply(apply((ou1.4[-1,]),2,as.numeric),2,mean)
ou2.5=apply(apply((ou1.5[-1,]),2,as.numeric),2,mean)
ou2.6=apply(apply((ou1.6[-1,]),2,as.numeric),2,mean)

out=rbind(ou2.1,ou2.2,ou2.3,ou2.4,ou2.5,ou2.6)

colnames(out)=c("NNS R2","np R2", "NNS np diff","Dependence")
rownames(out)=c(paste("DT",1:2,sep=""),paste("ST",1:4,sep=""))
print(c(reps=reps,n=n))
print(xtable(out,caption= n,
             caption.placement='top',
             table.placement="H"
             ,digits = xdigits(out)),type='latex')
}
