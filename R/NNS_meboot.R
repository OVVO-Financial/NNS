#' NNS meboot
#'
#' Adapted maximum entropy bootstrap routine from \code{meboot} \url{https://cran.r-project.org/package=meboot}.
#'
#' @param x vector of data, \code{ts} object or \code{pdata.frame} object.
#' @param reps numeric; number of replicates to generate.
#' @param setSpearman numeric [0,1]; The default setting \code{setSpearman = NULL} assumes that
#' the user does not want to generate replicates that are perfectly dependent on original time series, \code{setSpearman=1} recovers the original \code{meboot(...)} settings.
#' \code{setSpearman < 1} admits less perfect (more realistic for some purposes) dependence.
#' @param drift logical; \code{TRUE} default preserves the drift of the original series.
#' @param trim numeric [0,1]; The mean trimming proportion, defaults to \code{trim=0.1}.
#' @param xmin numeric; the lower limit for the left tail.
#' @param xmax numeric; the upper limit for the right tail.
#' @param reachbnd logical; If \code{TRUE} potentially reached bounds (xmin = smallest value - trimmed mean and
#' xmax = largest value + trimmed mean) are given when the random draw happens to be equal to 0 and 1, respectively.
#' @param expand.sd logical; If \code{TRUE} the standard deviation in the ensemble is expanded. See \code{expand.sd} in \code{meboot::meboot}.
#' @param force.clt logical; If \code{TRUE} the ensemble is forced to satisfy the central limit theorem. See \code{force.clt} in \code{meboot::meboot}.
#' @param scl.adjustment logical; If \code{TRUE} scale adjustment is performed to ensure that the population variance of the transformed series equals the variance of the data.
#' @param sym logical; If \code{TRUE} an adjustment is peformed to ensure that the ME density is symmetric.
#' @param elaps logical; If \code{TRUE} elapsed time during computations is displayed.
#' @param colsubj numeric; the column in \code{x} that contains the individual index. It is ignored if the input data \code{x} is not a \code{pdata.frame} object.
#' @param coldata numeric; the column in \code{x} that contains the data of the variable to create the ensemble. It is ignored if the input data \code{x} is not a \code{pdata.frame} object.
#' @param coltimes numeric; an optional argument indicating the column that contains the times at which the observations for each individual are observed. It is ignored if the input data \code{x}
#' is not a \code{pdata.frame} object.
#' @param ... possible argument \code{fiv} to be passed to \code{expand.sd}.
#'
#' @return
#' \itemize{
#'   \item{x} original data provided as input.
#' \item{replicates} maximum entropy bootstrap replicates.
#' \item{ensemble} average observation over all replicates.
#' \item{xx} sorted order stats (xx[1] is minimum value).
#' \item{z} class intervals limits.
#' \item{dv} deviations of consecutive data values.
#' \item{dvtrim} trimmed mean of dv.
#' \item{xmin} data minimum for ensemble=xx[1]-dvtrim.
#' \item{xmax} data x maximum for ensemble=xx[n]+dvtrim.
#' \item{desintxb} desired interval means.
#' \item{ordxx} ordered x values.
#' \item{kappa} scale adjustment to the variance of ME density.
#' \item{elaps} elapsed time.
#' }
#'
#' @references
#' \itemize{
#' \item Vinod, H.D. (2013), Maximum Entropy Bootstrap Algorithm Enhancements.
#' \href{http://ssrn.com/abstract=2285041}{http://ssrn.com/abstract=2285041}.
#'
#' \item Vinod, H.D. (2006), Maximum Entropy Ensembles for Time Series Inference in Economics,
#' \emph{Journal of Asian Economics}, \bold{17}(6), pp. 955-978.
#'
#' \item Vinod, H.D. (2004), Ranking mutual funds using unconventional utility theory and stochastic dominance, \emph{Journal of Empirical Finance}, \bold{11}(3), pp. 353-377.
#' }
#'
#' @examples
#' \dontrun{
#' # To generate an orthogonal rank correlated time-series to AirPassengers
#' boots <- NNS.meboot(AirPassengers, reps=100, setSpearman = 0, xmin = 0)
#'
#' # Verify correlation of replicates ensemble to original
#' cor(boots$ensemble, AirPassengers, method = "spearman")
#'
#' # Plot all replicates
#' matplot(boots$replicates, type = 'l')
#'
#' # Plot ensemble
#' lines(boots$ensemble, lwd = 3)
#'}
#'
#' @export

 NNS.meboot <- function(x,
                        reps=999,
                        setSpearman=1,
                        drift=TRUE,
                        trim=0.10,
                        xmin=NULL,
                        xmax=NULL,
                        reachbnd=TRUE,
                        expand.sd=TRUE, force.clt=TRUE,
                        scl.adjustment = FALSE, sym = FALSE, elaps=FALSE,
                        colsubj, coldata, coltimes,...)
  {
    if ("pdata.frame" %in% class(x))
    {
      res <- meboot::meboot.pdata.frame (x, reps, trim$trim, reachbnd,
                                 expand.sd, force.clt, scl.adjustment, sym, elaps,
                                 colsubj, coldata, coltimes, ...)
      return(res)
    }


    if(is.null(setSpearman)) setSpearman <- -99

    trim <- list(trim=trim, xmin=xmin, xmax=xmax)

    trimval <- if (is.null(trim$trim)) 0.1 else trim$trim


    ptm1 <- proc.time()

    n <- length(x)

    # Sort the original data in increasing order and
    # store the ordering index vector.

    xx <- sort(x)
    ordxx <- order(x)


    ### Fred Viole SUGGESTION PART 1 of 2

    if(setSpearman <1){
      ordxx_2 <- order(ordxx)
    }

    #ordxx <- sort.int(x, index.return=TRUE)

    # symmetry

    if (sym)
    {
      xxr <- rev(xx) #reordered values
      xx.sym <- mean(xx) + 0.5*(xx - xxr) #symmetrized order stats
      xx <- xx.sym #replace order stats by symmetrized ones
    }

    # Compute intermediate points on the sorted series.

    z <- rowMeans(embed(xx, 2))

    # Compute lower limit for left tail ('xmin') and
    # upper limit for right tail ('xmax').
    # This is done by computing the 'trim' (e.g. 10%) trimmed mean
    # of deviations among all consecutive observations ('dv').
    # Thus the tails are uniform distributed.

    dv <- abs(diff(as.numeric(x)))
    dvtrim <- mean(dv, trim=trimval)

    if (is.list(trim))
    {
      if (is.null(trim$xmin))
      {
        xmin <- xx[1] - dvtrim
      } else
        xmin <- trim$xmin

      if (is.null(trim$xmax))
      {
        xmax <- xx[n] + dvtrim
      } else
        xmax <- trim$xmax

      if (!is.null(trim$xmin) || !is.null(trim$xmax))
      {
        if (isTRUE(force.clt))
        {
          expand.sd <- FALSE
          force.clt <- FALSE
          warning("expand.sd and force.clt were set to FALSE in order to ",
                  "enforce the limits xmin/xmax.")
        }
      }
    } else {
      xmin <- xx[1] - dvtrim
      xmax <- xx[n] + dvtrim
    }

    # do this here so that this warnings are printed after
    # the above warnings (if necessary)

    if (is.list(trim))
    {
      if (!is.null(trim$xmin) && trim$xmin > min(x))
        warning("the lower limit trim$xmin may not be satisfied in the replicates ",
                "since it is higher than the minimum value observed ",
                "in the input series x")
      if (!is.null(trim$xmax) && trim$xmax < max(x))
        warning("the upper limit trim$xmax may not be satisfied in the replicates ",
                "since it is lower than the maximum value observed ",
                "in the input series x")
    }


    # Compute the mean of the maximum entropy density within each
    # interval in such a way that the 'mean preserving constraint'
    # is satisfied. (Denoted as m_t in the reference paper.)
    # The first and last interval means have distinct formulas.
    # See Theil and Laitinen (1980) for details.

    aux <- colSums( t(embed(xx, 3))*c(0.25,0.5,0.25) )
    desintxb <- c(0.75*xx[1]+0.25*xx[2], aux, 0.25*xx[n-1]+0.75*xx[n])

    # Generate random numbers from the [0,1] uniform interval and
    # compute sample quantiles at those points.

    # Generate random numbers from the [0,1] uniform interval.

    ensemble <- matrix(x, nrow=n, ncol=reps)
    ensemble <- apply(ensemble, 2, NNS.meboot.part,
                      n, z, xmin, xmax, desintxb, reachbnd)

    # So far the object 'ensemble' contains the quantiles.
    # Now give them time series dependence and heterogeneity.

    qseq <- apply(ensemble, 2, sort)


    # 'qseq' has monotonic series, the correct series is obtained
    # after applying the order according to 'ordxx' defined above.

    ensemble[ordxx,] <- qseq


    ### Pilot Spearman
    if(setSpearman==-99){
      y <- NNS.meboot(x, reps = 30, setSpearman = 1)$ensemble
      pilot <- cbind(x,y)
      setSpearman <-  fivenum(apply(pilot, 2, function(z) (cor(pilot[,1],z)))[-1])[2]
    }



    ### Fred Viole SUGGESTION  PART 2 of 2
    ### Average two ordxx ensemble matrices

    if(setSpearman<1){
      matrix2 = matrix(, nrow=length(x), ncol = reps)
      matrix2[ordxx_2,] = qseq

      # Intial search

      e <- c(ensemble)
      m <- c(matrix2)
      l <- length(e)

      func <- function(ab, d=drift){
        a <- ab[1]
        b <- ab[2]

        ifelse(d,
              (abs(cor((a*m + b*e)/(a + b), e, method = "spearman") - setSpearman) +
                  abs(mean((a*m + b*e))/mean(e) - 1) +
                    abs( cor((a*m + b*e)/(a + b), 1:l) - cor(e, 1:l))
              ),
              abs(cor((a*m + b*e)/(a + b), e, method = "spearman") - setSpearman) +
                abs(mean((a*m + b*e))/mean(e) - 1)
              )

      }

      res <- optim(c(.01,.01), func, control=list(abstol = .01))

      ensemble <- (res$par[1]*matrix2 +
                    res$par[2]*ensemble) / (sum(abs(res$par)))

      if(identical(ordxx_2, ordxx)){
        if(reps>1) ensemble <- t(apply(ensemble, 1, function(x) sample(x, size = reps, replace = TRUE)))
      }

    }


    if(expand.sd)
      ensemble <- NNS.meboot.expand.sd(x=x, ensemble=ensemble, ...)

    if(force.clt && reps > 1)
      ensemble <- meboot::force.clt(x=x, ensemble=ensemble)

    # scale adjustment

    if (scl.adjustment)
    {
      zz <- c(xmin,z,xmax) #extended list of z values
      v <- diff(zz^2) / 12
      xb <- mean(x)
      s1 <- sum((desintxb - xb)^2)
      uv <- (s1 + sum(v)) / n
      desired.sd <- sd(x)
      actualME.sd <- sqrt(uv)
      if (actualME.sd <= 0)
        stop("actualME.sd<=0 Error")
      out <- desired.sd / actualME.sd
      kappa <- out - 1

      ensemble <- ensemble + kappa * (ensemble - xb)
    } else
      kappa <- NULL


    # Force min / max values
    if(!is.null(trim[[2]])) ensemble <- apply(ensemble, 2, function(x) pmax(trim[[2]], x))
    if(!is.null(trim[[3]])) ensemble <- apply(ensemble, 2, function(x) pmin(trim[[3]], x))

    if(is.ts(x)){
      ensemble <- ts(ensemble, frequency=frequency(x), start=start(x))
      if(reps>1) dimnames(ensemble)[[2]] <- paste("Series", 1:reps)
    } else {
      if(reps>1) dimnames(ensemble)[[2]] <- paste("Replicate", 1:reps)
    }



    # Computation time
    ptm2 <- proc.time(); elapsr <- meboot::elapsedtime(ptm1, ptm2)
    if(elaps)
      cat("\n  Elapsed time:", elapsr$elaps,
          paste(elapsr$units, ".", sep=""), "\n")

    list(x=x, replicates=ensemble, ensemble=rowMeans(ensemble), xx=xx, z=z, dv=dv, dvtrim=dvtrim, xmin=xmin,
         xmax=xmax, desintxb=desintxb, ordxx=ordxx, kappa = kappa, elaps=elapsr)
  }
