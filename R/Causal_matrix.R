NNS.caus.matrix<- function(x, tau = tau){

n <- ncol(x)
    if(is.null(n)){
        stop("supply both 'x' and 'y' or a matrix-like 'x'")
    }

    indiv.causes <- list()

    for(i in 1 : (n - 1)){
        indiv.causes[[i]] <- sapply((i + 1) : n, function(b) NNS.caus(x[ , i], x[ , b], plot = FALSE, tau = tau))
    }


    for(i in 1 : length(indiv.causes)){
        indiv.causes[[i]] <- indiv.causes[[i]][2, ]-indiv.causes[[i]][1, ]
    }

    causes <- matrix(NA, n, n)
    causes[lower.tri(causes, diag = FALSE)] <- unlist(indiv.causes)


    diag(causes) <- 1
    causes <- pmax(causes, t(-causes), na.rm = TRUE)
    colnames(causes) <- colnames(x)
    rownames(causes) <- colnames(x)

    return(causes)

}
