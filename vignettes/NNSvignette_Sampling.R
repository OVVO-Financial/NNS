## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----setup2, message=FALSE, warning = FALSE-----------------------------------
library(NNS)
library(data.table)
require(knitr)
require(rgl)
require(meboot)

## -----------------------------------------------------------------------------
set.seed(123); x = rnorm(100)
ecdf(x)
P = ecdf(x)
P(0); P(1)

## ---- message=FALSE-----------------------------------------------------------
LPM.ratio(degree = 0, target = 0, variable = x); LPM.ratio(degree = 0, target = 1, variable = x)

## ---- fig.align='center', fig.width=6, fig.height=6, echo = FALSE-------------
LPM.CDF = LPM.ratio(degree = 0, target = sort(x), variable = x)

plot(ecdf(x))
points(sort(x), LPM.CDF, col='red')
legend('left', legend = c('ecdf', 'LPM.ratio'), fill=c('black','red'), border=NA, bty='n')

## ---- fig.align='center', fig.height=8, fig.width=8, echo=FALSE, warning=FALSE, message = FALSE----
zzz= rnorm(length(x), mean = 0, sd = 1)
norm_approx = pnorm(sort(zzz), mean=0, sd=1) #pnorm(sort(x),mean=-mean(x),sd=sd(x))

plot(ecdf(x), main = "eCDF via LPM.ratio()", lwd = 4)


# Altering shape of distribution with LPM degree
for(i in c(0, 0.25, .5, 1, 2)){
  idx <- which(i == c(0, 0.25, .5, 1, 2))
  lines(sort(x), LPM.ratio(i, sort(x),x), col = rainbow(5, alpha = 1)[idx], lty = 1, lwd = 3)
}

 lines(sort(zzz), norm_approx ,col='black', lty = 3, lwd = 2)


legend("topleft",c("LPM.ratio(degree = 0)","LPM.ratio(degree = 0.25)","LPM.ratio(degree = 0.5)","LPM.ratio(degree = 1)","LPM.ratio(degree = 2)", "N(0,1) approximation"),
       col = c(rainbow(5)[1:5], "black"), lwd = 3, lty = c(rep(1, 5), 3))

## ---- fig.align='center', echo=FALSE, fig.width=10, fig.height=8, message=FALSE, warning=FALSE----
layout(matrix(c(1, 1, 1,1,1,
                2, 3, 4,5,6,
                2, 3, 4,5,6), nrow=5, byrow=FALSE),widths = c(2,rep(1,5))) 
 
 
plot(ecdf(x), main = "eCDF via LPM.ratio()", lwd = 4)


# Altering shape of distribution with LPM degree
for(i in c(0, 0.25, .5, 1, 2)){
  idx <- which(i == c(0, 0.25, .5, 1, 2))
  lines(sort(x), LPM.ratio(i, sort(x),x), col = rainbow(5, alpha = 1)[idx], lty = 1, lwd = 3)
}

 lines(sort(zzz), norm_approx ,col='black', lty = 3, lwd = 2)


legend("topleft",c("LPM.ratio(degree = 0)","LPM.ratio(degree = 0.25)","LPM.ratio(degree = 0.5)","LPM.ratio(degree = 1)","LPM.ratio(degree = 2)", "N(0,1) approximation"),
       col = c(rainbow(5)[1:5], "black"), lwd = 3, lty = c(rep(1, 5), 3))




y = hist(LPM.VaR(seq(0,1,length.out = 100), 0, x), plot = FALSE, breaks = 15)

plot(y$breaks,
     c(y$counts,0), type = "s",
    col="black",lwd = 3, ylim = c(0,50), main = "Inverse CDF via LPM.VaR(degree 0)", breaks = 15, xlab = "x", ylab = "freq")
hist(LPM.VaR(seq(0,1,length.out = 100), 0, x), add = TRUE, col =  rainbow(5, alpha = .5)[1], breaks = 15)

y = hist(LPM.VaR(seq(0,1,length.out = 100), 0, x), border = NA, plot = FALSE, breaks = 15)
plot(y$breaks,
     c(y$counts,0)
     ,type="s",col="black",lwd = 3, ylim = c(0,50), main = "Inverse CDF via LPM.VaR(degree 0.25)", breaks = 15, xlab = "x", ylab = "freq")
hist(LPM.VaR(seq(0,1,length.out = 100), .25, x), border = rainbow(5)[2], add = TRUE, col =  rainbow(5, alpha = .5)[2], breaks = 15)

y = hist(LPM.VaR(seq(0,1,length.out = 100), 0, x), plot = FALSE, breaks = 15)
plot(y$breaks,
     c(y$counts,0)
     ,type="s",col="black",lwd = 3, ylim = c(0,50), main = "Inverse CDF via LPM.VaR(degree 0.5)", breaks = 15, xlab = "x", ylab = "freq")
hist(LPM.VaR(seq(0,1,length.out = 100), .5, x), border = rainbow(5)[3], add = TRUE, col =  rainbow(5, alpha = .5)[3], breaks = 15)

y = hist(LPM.VaR(seq(0,1,length.out = 100), 0, x), plot = FALSE, breaks = 15)
plot(y$breaks,
     c(y$counts,0)
     ,type="s",col="black",lwd = 3, ylim = c(0,50), main = "Inverse CDF via LPM.VaR(degree 1)", breaks = 15, xlab = "x", ylab = "freq")
hist(LPM.VaR(seq(0,1,length.out = 100), 1, x), border = rainbow(5)[4], add = TRUE, col =  rainbow(5, alpha = .5)[4], breaks = 15)

y = hist(LPM.VaR(seq(0,1,length.out = 100), 0, x), plot = FALSE, breaks = 15)
plot(y$breaks,
     c(y$counts,0)
     ,type="s",col="black",lwd = 3, ylim = c(0,50), main = "Inverse CDF via LPM.VaR(degree 2)", breaks = 15, xlab = "x", ylab = "freq")
hist(LPM.VaR(seq(0,1,length.out = 100), 2, x), border = rainbow(5)[5], add = TRUE, col =  rainbow(5, alpha = .5)[5], breaks = 15)

## -----------------------------------------------------------------------------
degree.0.samples = LPM.VaR(percentile = seq(0, 1, length.out = 100), degree = 0, x = x)
degree.0.25.samples = LPM.VaR(percentile = seq(0, 1, length.out = 100), degree = 0.25, x = x)
degree.0.5.samples = LPM.VaR(percentile = seq(0, 1, length.out = 100), degree = 0.5, x = x)
degree.1.samples = LPM.VaR(percentile = seq(0, 1, length.out = 100), degree = 1, x = x)
degree.2.samples = LPM.VaR(percentile = seq(0, 1, length.out = 100), degree = 2, x = x)

head(data.table::data.table(cbind("original x" = sort(x), degree.0.samples, degree.0.25.samples, degree.0.5.samples, degree.1.samples, degree.2.samples)), 10)

## ---- fig.align='center', fig.width=8, fig.height=8---------------------------
boots = NNS.MC(x, reps = 1, lower_rho = -1, upper_rho = 1, by = .25)$replicates
reps = do.call(cbind, boots)

plot(x, type = "l", lwd = 3, ylim = c(min(reps), max(reps)))
matplot(reps, type = "l", col = rainbow(length(boots)), add = TRUE)

## -----------------------------------------------------------------------------
sapply(boots, function(r) cor(r, x, method = "spearman"))

## ----multisim-----------------------------------------------------------------
set.seed(123)
x <- rnorm(1000); y <- rnorm(1000); z <- rnorm(1000)

# Add variable x to original data to avoid total independence (example only)
original.data <- cbind(x, y, z, x)

# Determine dependence structure
dep.structure <- apply(original.data, 2, function(x) LPM.ratio(degree = 1, target = x, variable = x))
  
# Generate new data with different mean, sd and length (or distribution type)
new.data <- sapply(1:ncol(original.data), function(x) rnorm(nrow(original.data)*2, mean = 10, sd = 20))

# Apply dependence structure to new data
new.dep.data <- sapply(1:ncol(original.data), function(x) LPM.VaR(percentile = dep.structure[,x], degree = 1, x = new.data[,x]))

## ----comparison, warning=FALSE------------------------------------------------
NNS.copula(original.data)
NNS.copula(new.dep.data)

## -----------------------------------------------------------------------------
head(original.data)
head(new.dep.data)

## -----------------------------------------------------------------------------
# Apply bootstrap to each variable
new.boot.dep.data = apply(original.data, 2, function(r) NNS.meboot(r, reps = 10, rho = .95))

# Reformat into vectors
boot.ensemble.vectors = lapply(new.boot.dep.data, function(z) unlist(z["ensemble",]))

# Create matrix from vectors
new.boot.dep.matrix = do.call(cbind, boot.ensemble.vectors)

## -----------------------------------------------------------------------------
for(i in 1:4) print(cor(new.boot.dep.matrix[,i], original.data[,i], method = "spearman"))

## -----------------------------------------------------------------------------
NNS.copula(original.data)
NNS.copula(new.boot.dep.matrix)

## -----------------------------------------------------------------------------
head(original.data)
head(new.boot.dep.matrix)

