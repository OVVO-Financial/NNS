## ----setup, include=FALSE, message=FALSE--------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(NNS)
options(mc.cores = 1)
RcppParallel::setThreadOptions(numThreads = 1)
Sys.setenv("OMP_THREAD_LIMIT" = 1)

## ----setup2, message=FALSE, warning = FALSE-----------------------------------
library(NNS)
require(knitr)
require(rgl)

## -----------------------------------------------------------------------------
set.seed(123); x = rnorm(100)
ecdf(x)
P = ecdf(x)
P(0); P(1)

## ----message=FALSE------------------------------------------------------------
LPM.ratio(degree = 0, target = 0, variable = x); LPM.ratio(degree = 0, target = 1, variable = x)

## ----fig.align='center', fig.width=6, fig.height=6, echo = FALSE--------------
LPM.CDF = LPM.ratio(degree = 0, target = sort(x), variable = x)

plot(ecdf(x))
points(sort(x), LPM.CDF, col='red')
legend('left', legend = c('ecdf', 'LPM.ratio'), fill=c('black','red'), border=NA, bty='n')

## ----fig.align='center', fig.height=8, fig.width=8, echo=FALSE, warning=FALSE, message = FALSE, eval=FALSE----
# zzz = rnorm(length(x), mean = 0, sd = 1)
# norm_approx = pnorm(sort(zzz), mean=0, sd=1) #pnorm(sort(x),mean=-mean(x),sd=sd(x))
# 
# plot(ecdf(x), main = "eCDF via LPM.ratio()", lwd = 4)
# 
# 
# # Altering shape of distribution with LPM degree
# for(i in c(0, 0.25, .5, 1, 2)){
#   idx <- which(i == c(0, 0.25, .5, 1, 2))
#   lines(sort(x), LPM.ratio(i, sort(x),x), col = rainbow(5, alpha = 1)[idx], lty = 1, lwd = 3)
# }
# 
#  lines(sort(zzz), norm_approx ,col='black', lty = 3, lwd = 2)
# 
# 
# legend("topleft",c("LPM.ratio(degree = 0)","LPM.ratio(degree = 0.25)","LPM.ratio(degree = 0.5)","LPM.ratio(degree = 1)","LPM.ratio(degree = 2)", "N(0,1) approximation"),
#        col = c(rainbow(5)[1:5], "black"), lwd = 3, lty = c(rep(1, 5), 3))

## ----fig.align='center', echo=FALSE, fig.width=10, fig.height=8, message=FALSE, warning=FALSE, eval=FALSE----
# layout(matrix(c(1, 1, 1,1,1,
#                 2, 3, 4,5,6,
#                 2, 3, 4,5,6), nrow=5, byrow=FALSE),widths = c(2,rep(1,5)))
# 
# 
# plot(ecdf(x), main = "eCDF via LPM.ratio()", lwd = 4)
# 
# 
# # Altering shape of distribution with LPM degree
# for(i in c(0, 0.25, .5, 1, 2)){
#   idx <- which(i == c(0, 0.25, .5, 1, 2))
#   lines(sort(x), LPM.ratio(i, sort(x),x), col = rainbow(5, alpha = 1)[idx], lty = 1, lwd = 3)
# }
# 
#  lines(sort(zzz), norm_approx ,col='black', lty = 3, lwd = 2)
# 
# 
# legend("topleft",c("LPM.ratio(degree = 0)","LPM.ratio(degree = 0.25)","LPM.ratio(degree = 0.5)","LPM.ratio(degree = 1)","LPM.ratio(degree = 2)", "N(0,1) approximation"),
#        col = c(rainbow(5)[1:5], "black"), lwd = 3, lty = c(rep(1, 5), 3))
# 
# 
# 
# 
# y = hist(LPM.VaR(seq(0,1,length.out = 100), 0, x), plot = FALSE, breaks = 15)
# 
# plot(y$breaks,
#      c(y$counts,0), type = "s",
#     col="black",lwd = 3, ylim = c(0,50), main = "Inverse CDF via LPM.VaR(degree 0)", breaks = 15, xlab = "x", ylab = "freq")
# hist(LPM.VaR(seq(0,1,length.out = 100), 0, x), add = TRUE, col =  rainbow(5, alpha = .5)[1], breaks = 15)
# 
# y = hist(LPM.VaR(seq(0,1,length.out = 100), 0, x), border = NA, plot = FALSE, breaks = 15)
# plot(y$breaks,
#      c(y$counts,0)
#      ,type="s",col="black",lwd = 3, ylim = c(0,50), main = "Inverse CDF via LPM.VaR(degree 0.25)", breaks = 15, xlab = "x", ylab = "freq")
# hist(LPM.VaR(seq(0,1,length.out = 100), .25, x), border = rainbow(5)[2], add = TRUE, col =  rainbow(5, alpha = .5)[2], breaks = 15)
# 
# y = hist(LPM.VaR(seq(0,1,length.out = 100), 0, x), plot = FALSE, breaks = 15)
# plot(y$breaks,
#      c(y$counts,0)
#      ,type="s",col="black",lwd = 3, ylim = c(0,50), main = "Inverse CDF via LPM.VaR(degree 0.5)", breaks = 15, xlab = "x", ylab = "freq")
# hist(LPM.VaR(seq(0,1,length.out = 100), .5, x), border = rainbow(5)[3], add = TRUE, col =  rainbow(5, alpha = .5)[3], breaks = 15)
# 
# y = hist(LPM.VaR(seq(0,1,length.out = 100), 0, x), plot = FALSE, breaks = 15)
# plot(y$breaks,
#      c(y$counts,0)
#      ,type="s",col="black",lwd = 3, ylim = c(0,50), main = "Inverse CDF via LPM.VaR(degree 1)", breaks = 15, xlab = "x", ylab = "freq")
# hist(LPM.VaR(seq(0,1,length.out = 100), 1, x), border = rainbow(5)[4], add = TRUE, col =  rainbow(5, alpha = .5)[4], breaks = 15)
# 
# y = hist(LPM.VaR(seq(0,1,length.out = 100), 0, x), plot = FALSE, breaks = 15)
# plot(y$breaks,
#      c(y$counts,0)
#      ,type="s",col="black",lwd = 3, ylim = c(0,50), main = "Inverse CDF via LPM.VaR(degree 2)", breaks = 15, xlab = "x", ylab = "freq")
# hist(LPM.VaR(seq(0,1,length.out = 100), 2, x), border = rainbow(5)[5], add = TRUE, col =  rainbow(5, alpha = .5)[5], breaks = 15)

## ----eval=FALSE---------------------------------------------------------------
# degree.0.samples = LPM.VaR(percentile = seq(0, 1, length.out = 100), degree = 0, x = x)
# degree.0.25.samples = LPM.VaR(percentile = seq(0, 1, length.out = 100), degree = 0.25, x = x)
# degree.0.5.samples = LPM.VaR(percentile = seq(0, 1, length.out = 100), degree = 0.5, x = x)
# degree.1.samples = LPM.VaR(percentile = seq(0, 1, length.out = 100), degree = 1, x = x)
# degree.2.samples = LPM.VaR(percentile = seq(0, 1, length.out = 100), degree = 2, x = x)
# 
# head(data.frame(check.names = FALSE, cbind("original x" = sort(x), degree.0.samples,
#                                                           degree.0.25.samples,
#                                                           degree.0.5.samples,
#                                                           degree.1.samples,
#                                                           degree.2.samples)), 10)
# 
# #>           original x degree.0.samples degree.0.25.samples degree.0.5.samples degree.1.samples
# #> 0.000000%  -2.309169        -2.309169           -2.309169         -2.3091689       -2.3091689
# #> 1.010101%  -1.966617        -1.966617           -1.941187         -1.6935534       -1.4744855
# #> 2.020202%  -1.686693        -1.686693           -1.599486         -1.4541184       -1.2160087
# #> 3.030303%  -1.548753        -1.548753           -1.382558         -1.2462537       -1.0823025
# #> 4.040404%  -1.265396        -1.265396           -1.250804         -1.1454066       -0.9968077
# #> 5.050505%  -1.265061        -1.265061           -1.176437         -1.0745549       -0.9290647
# #> 6.060606%  -1.220718        -1.220718           -1.119649         -1.0252637       -0.8666875
# #> 7.070707%  -1.138137        -1.138137           -1.067817         -0.9868700       -0.8090628
# #> 8.080808%  -1.123109        -1.123109           -1.026420         -0.9322146       -0.7556675
# #> 9.090909%  -1.071791        -1.071791           -1.014290         -0.8710985       -0.7070109
# #>           degree.2.samples
# #> 0.000000%       -2.3091689
# #> 1.010101%       -1.1614823
# #> 2.020202%       -0.9710028
# #> 3.030303%       -0.8610022
# #> 4.040404%       -0.7810056
# #> 5.050505%       -0.7169391
# #> 6.060606%       -0.6632142
# #> 7.070707%       -0.6170757
# #> 8.080808%       -0.5765594
# #> 9.090909%       -0.5403627

## ----fig.align='center', fig.width=8, fig.height=8, eval=FALSE----------------
# set.seed(123)
# boots = NNS.MC(x, reps = 1, lower_rho = -1, upper_rho = 1, by = .5)$replicates
# reps = do.call(cbind, boots)
# 
# 
# matplot(reps, type = "l", col = rainbow(length(boots)))
# lines(x, type = "l", lwd = 3, ylim = c(min(reps), max(reps)))

## ----eval = FALSE-------------------------------------------------------------
# sapply(boots, function(r) cor(r, x, method = "spearman"))
# 
# #>     rho = 1   rho = 0.5     rho = 0  rho = -0.5    rho = -1
# #>  0.99840384  0.50357036  0.01478548 -0.46829883 -0.97912991

## ----tgt_drift, fig.align='center', fig.width=8, fig.height=8, eval=FALSE-----
# boots = NNS.MC(x, reps = 1, lower_rho = -1, upper_rho = 1, by = .5, target_drift = 0.05)$replicates
# reps = do.call(cbind, boots)
# 
# plot(x, type = "l", lwd = 3, ylim = c(min(c(x, reps)), max(c(x, reps))))
# matplot(reps, type = "l", col = rainbow(length(boots)), add = TRUE)

## ----multisim, eval=FALSE-----------------------------------------------------
# set.seed(123)
# x <- rnorm(1000); y <- rnorm(1000); z <- rnorm(1000)
# original.data <- cbind(x, y, z, x_copy = x)
# n <- nrow(original.data)
# dep.structure <- apply(original.data, 2, function(v) {
#   LPM.ratio(degree = 0, target = v, variable = v) - .5 / n
# })
# new.data <- replicate(ncol(original.data), rnorm(2000, mean = 10, sd = 20))
# new.dep.data <- vapply(seq_len(ncol(original.data)), function(j) {
#   LPM.VaR(percentile = dep.structure[, j], degree = 0, x = new.data[, j])
# }, numeric(n))
# 
# cor(original.data, method = "spearman")
# cor(new.dep.data, method = "spearman")

## ----eval=FALSE---------------------------------------------------------------
# # Apply bootstrap to each variable
# new.boot.dep.data = apply(original.data, 2, function(r) NNS.meboot(r, reps = 100, rho = .95))
# 
# # Reformat into vectors
# boot.ensemble.vectors = lapply(new.boot.dep.data, function(z) unlist(z["ensemble",]))
# 
# # Create matrix from vectors
# new.boot.dep.matrix = do.call(cbind, boot.ensemble.vectors)

## ----eval=FALSE---------------------------------------------------------------
# for(i in 1:4) print(cor(new.boot.dep.matrix[,i], original.data[,i], method = "spearman"))
# 
# #> [1] 0.9885851
# #> [1] 0.9965075
# #> [1] 0.9906302
# #> [1] 0.9901535

## ----eval=FALSE---------------------------------------------------------------
# NNS.copula(original.data)
# NNS.copula(new.boot.dep.matrix)
# 
# #> [1] 0.4452914
# #> [1] 0.4423774

## ----eval=FALSE---------------------------------------------------------------
# head(original.data)
# head(new.boot.dep.matrix)
# 
# #>                x           y          z      x_copy
# #> [1,] -0.56047565 -0.99579872 -0.5116037 -0.56047565
# #> [2,] -0.23017749 -1.03995504  0.2369379 -0.23017749
# #> [3,]  1.55870831 -0.01798024 -0.5415892  1.55870831
# #> [4,]  0.07050839 -0.13217513  1.2192276  0.07050839
# #> [5,]  0.12928774 -2.54934277  0.1741359  0.12928774
# #> [6,]  1.71506499  1.04057346 -0.6152683  1.71506499
# #>                     x           y          z      x_copy
# #> ensemble1 -0.50534918 -0.99114837 -0.5599005 -0.53374834
# #> ensemble2 -0.25336404 -1.10495674  0.1277390 -0.20297478
# #> ensemble3  1.62035927  0.07462431 -0.5952796  1.61734635
# #> ensemble4 -0.02471890  0.01657279  1.2773063  0.05897732
# #> ensemble5  0.05586323 -1.48844305  0.1384859  0.12920777
# #> ensemble6  1.80768833  0.82770106 -0.6056424  1.83212766

