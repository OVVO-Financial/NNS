## ----setup, include=FALSE, message=FALSE--------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(NNS)
options(mc.cores = 1)
RcppParallel::setThreadOptions(numThreads = 1)
Sys.setenv("OMP_THREAD_LIMIT" = 1)


## ----setup2, message=FALSE, warning=FALSE-------------------------------------
library(NNS)
require(knitr)


## ----linear-------------------------------------------------------------------
x = seq(-5, 5, .05); y = x ^ 3

for(i in 1 : 4){NNS.part(x, y, order = i, Voronoi = TRUE, obs.req = 0)}


## ----x part,results='hide'----------------------------------------------------
for(i in 1 : 4){NNS.part(x, y, order = i, type = "XONLY", Voronoi = TRUE)}


## ----res2, echo=FALSE---------------------------------------------------------
NNS.part(x,y,order = 4, type = "XONLY")


## ----depreg,results='hide'----------------------------------------------------
for(i in 1 : 3){NNS.part(x, y, order = i, obs.req = 0, Voronoi = TRUE, type = "XONLY") ; NNS.reg(x, y, order = i, ncores = 1)}


## ----nonlinear,fig.width=5,fig.height=3,fig.align = "center"------------------
NNS.reg(x, y, ncores = 1)


## ----nonlinear multi,fig.width=5,fig.height=3,fig.align = "center"------------
f = function(x, y) x ^ 3 + 3 * y - y ^ 3 - 3 * x
y = x ; z <- expand.grid(x, y)
g = f(z[ , 1], z[ , 2])
NNS.reg(z, g, order = "max", plot = FALSE, ncores = 1)


## ----dimension_reduction,fig.width=5,fig.height=3,fig.align="center"----------
regression_x <- iris[, 2:4]
regression_y <- iris$Sepal.Length
reduced <- NNS.reg(regression_x, regression_y, dim.red.method = "cor",
                   plot = FALSE, ncores = 1)
reduced$equation


## ----dimension_reduction_threshold,fig.width=5,fig.height=3,fig.align="center"----
reduced_threshold <- NNS.reg(regression_x, regression_y,
                            dim.red.method = "cor", threshold = .75,
                            point.est = regression_x[1:10, ],
                            plot = FALSE, ncores = 1)
reduced_threshold$equation
reduced_threshold$Point.est


## ----class,fig.width=5,fig.height=3,fig.align="center",message=FALSE----------
classification <- NNS.reg(iris[, 1:4], iris$Species, type = "CLASS",
                          point.est = iris[1:10, 1:4], plot = FALSE, ncores = 1)
classification$class.levels[classification$Point.est]
classification$Point.probabilities


## ----classification_indices---------------------------------------------------
class_reduced <- NNS.reg(iris[, 1:4], iris$Species, dim.red.method = "cor",
                         plot = FALSE, ncores = 1)
class_reduced$equation


## ----stack,fig.width=5,fig.height=3,fig.align="center",message=FALSE----------
set.seed(123)
stack_fit <- NNS.stack(IVs.train = iris[, 1:4],
                       DV.train = iris$Species,
                       IVs.test = iris[1:10, 1:4],
                       dim.red.method = "cor",
                       obj.fn = expression(mean(predicted == actual)),
                       objective = "max", type = "CLASS",
                       folds = 1, ncores = 1, status = FALSE)
stack_fit[c("OBJfn.reg", "OBJfn.dim.red", "NNS.reg.n.best",
            "NNS.dim.red.count", "weights", "class.levels")]
stack_fit$class.levels[stack_fit$stack]


## ----stack2, message = FALSE,fig.width=5,fig.height=3,fig.align = "center",results='hide'----
set.seed(123)
x = rnorm(100); y = rnorm(100)

nns.params = NNS.stack(IVs.train = cbind(x, x),
                        DV.train = y,
                        method = 1, ncores = 1)


## ----stack2res, fig.width=5,fig.height=3,fig.align = "center",results='hide'----
duplicated_fit <- NNS.reg(cbind(x, x), y,
                          n.best = nns.params$NNS.reg.n.best,
                          point.est = cbind(x, x), plot = FALSE,
                          residual.plot = FALSE, ncores = 1,
                          confidence.interval = .95)
plot(x, y, pch = 19, col = "grey60", main = "Tuned multivariate prediction path")
o <- order(x)
lines(x[o], duplicated_fit$Point.est[o], col = "steelblue", lwd = 2)


## ----smooth, fig.width=5,fig.height=3,fig.align = "center",results='hide'-----
NNS.reg(x, y, smooth = TRUE)


## ----uniimpute, fig.width=7, fig.height=5-------------------------------------
set.seed(123)

# Univariate predictor with nonlinear signal
n <- 400
x <- sort(runif(n, -3, 3))
y <- sin(x) + 0.2 * x^2 + rnorm(n, 0, 0.25)

# Induce ~25% MCAR missingness in y
miss <- rbinom(n, 1, 0.25) == 1
y_mis <- y
y_mis[miss] <- NA

# ---- Increasing dimensions trick ----
# Duplicate x so the distance operates in a 2D space: cbind(x, x).
# This chooses the multivariate neighbor path; duplication adds no information.
x2_train <- cbind(x[!miss], x[!miss])
x2_miss  <- cbind(x[miss],  x[miss])

# 1-NN donor imputation with NNS.reg
y_hat_uni <- NNS::NNS.reg(
  x         = x2_train,             # predictors (duplicated x)
  y         = y[!miss],             # observed responses
  point.est = x2_miss,              # rows to impute
  order     = "max",                # maximum requested partition order
  n.best    = 1,                    # 1-NN donor
  plot      = FALSE
)$Point.est

# Fill back
y_completed_uni <- y_mis
y_completed_uni[miss] <- y_hat_uni

# Plot observed vs imputed (NNS 1-NN)
plot(x, y, pch = 1, col = "steelblue", cex = 1.5, lwd = 2,
     xlab = "x", ylab = "y", main = "NNS 1-NN Imputation")
points(x[miss], y_hat_uni, col = "red", pch = 15, cex = 1.3)

legend("topleft",
       legend = c("Observed", "Imputed (NNS 1-NN)"),
       col    = c("steelblue", "red"),
       pch    = c(1, 15),
       pt.lwd = c(2, NA),
       bty    = "n")


## ----multiimpute, fig.width=7, fig.height=5-----------------------------------
set.seed(123)

# Multivariate predictors with nonlinear & interaction structure
n <- 800
X <- cbind(
  x1 = rnorm(n),
  x2 = runif(n, -2, 2),
  x3 = rnorm(n, 0, 1)
)

f <- function(x1, x2, x3) 1.1*x1 - 0.8*x2 + 0.5*x3 + 0.6*x1*x2 - 0.4*x2*x3 + 0.3*sin(1.3*x1)
y <- f(X[,1], X[,2], X[,3]) + rnorm(n, 0, 0.4)

# Induce ~30% MCAR missingness in y
miss <- rbinom(n, 1, 0.30) == 1
y_mis <- y
y_mis[miss] <- NA

# Training (observed) vs rows to impute
X_obs <- X[!miss, , drop = FALSE]
y_obs <- y[!miss]
X_mis <- X[ miss, , drop = FALSE]

# 1-NN donor imputation with NNS.reg
y_hat_mv <- NNS::NNS.reg(
  x         = X_obs,     # all observed predictors
  y         = y_obs,     # observed responses
  point.est = X_mis,     # rows to impute
  order     = "max",     # maximum requested partition order
  n.best    = 1,         # 1-NN donor
  plot      = FALSE
)$Point.est

# Completed vector
y_completed_mv <- y_mis
y_completed_mv[miss] <- y_hat_mv

# Plot observed vs imputed (multivariate, NNS 1-NN)
plot(seq_along(y), y,
     pch = 1, col = "steelblue", cex = 1.5, lwd = 2,
     xlab = "Observation index", ylab = "y",
     main = "NNS 1-NN Multivariate Imputation")

# Overlay imputed values
points(which(miss), y_hat_mv, pch = 15, col = "red", cex = 1.2)

# Legend
legend("topleft",
       legend = c("Observed", "Imputed (NNS 1-NN)"),
       col    = c("steelblue", "red"),
       pch    = c(1, 15),
       pt.lwd = c(2, NA),
       bty    = "n")

