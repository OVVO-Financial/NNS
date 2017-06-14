NNS.dep.matrix = function(x, order = NULL,
                   degree= NULL){

n= ncol(x)
if(is.null(n)){stop("supply both 'x' and 'y' or a matrix-like 'x'")}

raw.rhos = list()
raw.deps = list()

for(i in 1:(n-1)){
        raw.rhos[[i]]=sapply((i+1):n, function(b) NNS.dep(x[,i],x[,b],print.map = F,order=order,degree = degree)$Correlation)
        raw.deps[[i]]=sapply((i+1):n, function(b) NNS.dep(x[,i],x[,b],print.map = F,order=order,degree = degree)$Dependence)
        }


rhos <- matrix(, n, n)
rhos[lower.tri(rhos, diag=FALSE)] <- unlist(raw.rhos)
diag(rhos) <- 1
rhos = pmax(rhos, t(rhos), na.rm=TRUE)

deps <- matrix(, n, n)
deps[lower.tri(deps, diag=FALSE)] <- unlist(raw.deps)
diag(deps) <- 1
deps = pmax(deps, t(deps), na.rm=TRUE)

colnames(rhos) = colnames(x);colnames(deps) = colnames(x)
rownames(rhos) = colnames(x);rownames(deps) = colnames(x)

return(list("Correlation"=rhos,"Dependence"=deps))
}



