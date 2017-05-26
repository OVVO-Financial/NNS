
NNS.dep.matrix = function(x, order = NULL,
                   degree= NULL){

n= ncol(x)
if(is.null(n)){stop("supply both 'x' and 'y' or a matrix-like 'x'")}

rhos = list()
deps = list()

for(i in 1:n){
        rhos[[i]]=sapply(1:n, function(b) NNS.dep(x[,i],x[,b],print.map = F,order=order,degree = degree)$Correlation)
        deps[[i]]=sapply(1:n, function(b) NNS.dep(x[,i],x[,b],print.map = F,order=order,degree = degree)$Dependence)
        }

rhos=do.call(rbind,rhos)
deps=do.call(rbind,deps)

colnames(rhos) = colnames(x);colnames(deps) = colnames(x)
rownames(rhos) = colnames(x);rownames(deps) = colnames(x)

return(list("Correlation"=rhos,"Dependence"=deps))
}



