## ----setup, message=FALSE-----------------------------------------------------
# Prereqs (uncomment if needed):
# install.packages("NNS")

library(NNS)

## ----include=FALSE, message=FALSE---------------------------------------------
options(mc.cores = 1)
RcppParallel::setThreadOptions(numThreads = 1)
Sys.setenv("OMP_THREAD_LIMIT" = 1)

## -----------------------------------------------------------------------------
set.seed(42)

# Normal sample
y <- rnorm(3000)
mu <- mean(y)
L2 <- LPM(2, mu, y); U2 <- UPM(2, mu, y)
cat(sprintf("LPM2 + UPM2 = %.6f vs var(y)=%.6f\n", (L2+U2)*(length(y) / (length(y) - 1)), var(y)))

# Empirical CDF via LPM.ratio(0, t, x)
for (t in c(-1,0,1)) {
  cdf_lpm <- LPM.ratio(0, t, y)
  cat(sprintf("CDF at t=%+.1f : LPM.ratio=%.4f | empirical=%.4f\n", t, cdf_lpm, mean(y<=t)))
}

# Asymmetry on a skewed distribution
z <- rexp(3000)-1; mu_z <- mean(z)
cat(sprintf("Skewed z: LPM2=%.4f, UPM2=%.4f (expect imbalance)\n", LPM(2,mu_z,z), UPM(2,mu_z,z)))

## -----------------------------------------------------------------------------
M <- NNS.moments(y)
M

## -----------------------------------------------------------------------------
set.seed(23)
multimodal <- c(rnorm(1500,-2,.5), rnorm(1500,2,.5))
NNS.mode(multimodal,multi = TRUE)

## -----------------------------------------------------------------------------
qgrid <- LPM.VaR(seq(0.05,0.95,.1),0,z) # equivalent to quantile(z,probs = seq(0.05,0.95,by=0.1))
CDF_tbl <- data.frame(threshold = as.numeric(qgrid), CDF = LPM.ratio(0,qgrid,z))
CDF_tbl

## -----------------------------------------------------------------------------
set.seed(1)
x <- runif(2000,-1,1)
y <- x^2 + rnorm(2000, sd=.05)
cat(sprintf("Pearson r = %.4f\n", cor(x,y)))
cat(sprintf("NNS.dep  = %.4f\n", NNS.dep(x,y)$Dependence))

X <- data.frame(a=x, b=y, c=x*y + rnorm(2000, sd=.05))
pm <- PM.matrix(1, 1, target = "means", variable=X, pop_adj=TRUE)
pm

cop <- NNS.copula(X, continuous=TRUE, plot=FALSE)
cop

## ----eval=FALSE---------------------------------------------------------------
# # Data
# set.seed(123); x = rnorm(100); y = rnorm(100); z = expand.grid(x, y)
# 
# # Plot
# rgl::plot3d(z[,1], z[,2], Co.LPM(0, z[,1], z[,2], z[,1], z[,2]), col = "red")
# 
# # Uniform values
# u_x = LPM.ratio(0, x, x); u_y = LPM.ratio(0, y, y); z = expand.grid(u_x, u_y)
# 
# # Plot
# rgl::plot3d(z[,1], z[,2], Co.LPM(0, z[,1], z[,2], z[,1], z[,2]), col = "blue")

## -----------------------------------------------------------------------------
A <- rnorm(100, mean = 0, sd = 1)
B <- rnorm(100, mean = 0, sd = 5)
C <- rnorm(100, mean = 10, sd = 1)
D <- rnorm(100, mean = 10, sd = 10)

X <- data.frame(A, B, C, D)

# Linear scaling
lin_norm <- NNS.norm(X, linear = TRUE, chart.type=NULL, location=NULL)

## -----------------------------------------------------------------------------
px <- 100 + cumsum(rnorm(260, sd = 1))
rn <- NNS.rescale(px, a=100, b=0.03, method="riskneutral", T=1, type="Terminal")
c( target = 100*exp(0.03*1), mean_rn = mean(rn) )

## -----------------------------------------------------------------------------
ctrl <- rnorm(200, 0, 1)
trt  <- rnorm(180, 0.35, 1.2)
NNS.ANOVA(control=ctrl, treatment=trt, means.only=FALSE, plot=FALSE)

A <- list(g1=rnorm(150,0.0,1.1), g2=rnorm(150,0.2,1.0), g3=rnorm(150,-0.1,0.9))
NNS.ANOVA(control=A, means.only=TRUE, plot=FALSE)

## ----stochsuperiority, echo=TRUE----------------------------------------------
set.seed(123)
x = rnorm(1000, mean = 0, sd = 1)
y = rnorm(1000, mean = 1, sd = 1)

NNS.SS(x, y)

## ----stochsuperiorityci, echo=TRUE, eval=FALSE--------------------------------
# NNS.SS(x, y, confidence.interval = TRUE, reps = 999, ci = 0.95)[1:5]
# 
# $p_gt
# [1] 0.233915
# 
# $p_tie
# [1] 0
# 
# $p_star
# [1] 0.233915
# 
# $lower
# [1] 0.2105631
# 
# $upper
# [1] 0.2537789

## ----stochsuperioritydiscrete, echo=TRUE--------------------------------------
set.seed(123)
x = sample(1:5, 100, replace = TRUE)
y = sample(1:5, 100, replace = TRUE)

NNS.SS(x, y)

## ----fig.width=7, fig.height=5, fig.align='center'----------------------------
# Example 1: Nonlinear regression
set.seed(123)
x_train <- runif(1000, -2, 2)
y_train <- sin(pi * x_train) + rnorm(1000, sd = 0.2)

x_test <- seq(-2, 2, length.out = 100)

NNS.reg(x = x_train, y = y_train, order = NULL, point.est = x_test)

## ----eval = FALSE-------------------------------------------------------------
# # Simple train/test for boosting & stacking
# test.set = 141:150
# 
# boost <- NNS.boost(IVs.train = iris[-test.set, 1:4],
#               DV.train = iris[-test.set, 5],
#               IVs.test = iris[test.set, 1:4],
#               epochs = 10, learner.trials = 10,
#               status = FALSE, balance = TRUE,
#               type = "CLASS", folds = 5)
# 
# 
# mean(boost$results == as.numeric(iris[test.set,5]))
# # [1] 1
# 
# 
# boost$feature.weights; boost$feature.frequency
# 
# stacked <- NNS.stack(IVs.train = iris[-test.set, 1:4],
#                      DV.train = iris[-test.set, 5],
#                      IVs.test = iris[test.set, 1:4],
#                      type = "CLASS", balance = TRUE,
#                      ncores = 1, folds = 1)
# mean(stacked$stack == as.numeric(iris[test.set,5]))
# # [1] 1

## -----------------------------------------------------------------------------
NNS.caus(mtcars$hp,  mtcars$mpg)  # hp -> mpg
NNS.caus(mtcars$mpg, mtcars$hp)   # hp -> mpg

## ----fig.width=7, fig.align='center'------------------------------------------
# Univariate nonlinear ARMA
set.seed(42)
z <- as.numeric(scale(sin(1:480/8) + rnorm(480, sd=.35)))

# Seasonality detection (prints a summary)
seasonal_period <- NNS.seas(z, plot = FALSE)
head(seasonal_period$all.periods)

# Validate seasonal periods
NNS.ARMA.optim(z, h = 48, seasonal.factor = seasonal_period$periods, plot = TRUE, ncores = 1)

## -----------------------------------------------------------------------------
x_ts <- cumsum(rnorm(350, sd=.7))
mb <- NNS.meboot(x_ts, reps=5, rho = 1)
dim(mb["replicates", ]$replicates)

## -----------------------------------------------------------------------------
mc <- NNS.MC(x_ts, reps=5, lower_rho=-1, upper_rho=1, by=.5, exp=1)
length(mc$ensemble); names(mc$replicates)

head(mc$replicates$`rho = 0`)

## -----------------------------------------------------------------------------
set.seed(42)
RA <- rnorm(240, 0.005, 0.03)
RB <- rnorm(240, 0.003, 0.02)
RC <- rnorm(240, 0.006, 0.04)

NNS.FSD.uni(RA, RB)
NNS.SSD.uni(RA, RB)
NNS.TSD.uni(RA, RB)

Rmat <- cbind(A=RA, B=RB, C=RC)
try(NNS.SD.cluster(Rmat, degree = 1))
try(NNS.SD.efficient.set(Rmat, degree = 1))

