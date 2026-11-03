## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7, fig.height = 5)
suppressPackageStartupMessages(library(NNS))
options(mc.cores = 1)
RcppParallel::setThreadOptions(numThreads = 1)
Sys.setenv("OMP_THREAD_LIMIT" = 1)

## ----install,message=FALSE,warning = FALSE------------------------------------
library(NNS)
require(knitr)
require(rgl)

## ----basic-example, eval=FALSE------------------------------------------------
# set.seed(123)
# 
# A <- rnorm(100, mean = 5, sd = 1)
# B <- rnorm(100, mean = 20, sd = 5)
# C <- rnorm(100, mean = 10, sd = 1)
# D <- rnorm(100, mean = 50, sd = 10)
# 
# X <- data.frame(A, B, C, D)
# 
# # Linear scaling
# lin_norm <- NNS.norm(X, linear = TRUE, chart.type = NULL)
# head(lin_norm)
#      A Normalized B Normalized C Normalized D Normalized
# [1,]     18.38264     17.81323     25.40626     18.19440
# [2,]     19.75029     23.05113     23.56018     18.03539
# [3,]     27.15749     20.32425     20.27462     17.24622
# [4,]     20.99534     19.77814     21.95814     16.76224
# [5,]     21.23872     16.50706     19.96390     19.37523
# [6,]     27.80492     21.41627     19.83496     22.63783
# 
# # Verify means are equal
# apply(lin_norm, 2, function(x) c(mean = mean(x), sd = sd(x)))
# 
#      A Normalized B Normalized C Normalized D Normalized
# mean    21.077727    21.077727    21.077727    21.077727
# sd       3.779676     5.236256     1.978297     4.410985

## ----nonlinear-example, eval=FALSE--------------------------------------------
# nonlin_norm <- NNS.norm(X, linear = FALSE, chart.type = NULL)
# head(nonlin_norm)
#      A Normalized B Normalized C Normalized D Normalized
# [1,]     2.082162     4.690344     4.098285     11.04240
# [2,]     2.237074     6.069519     3.800493     10.94589
# [3,]     3.076072     5.351514     3.270500     10.46693
# [4,]     2.378098     5.207720     3.542068     10.17321
# [5,]     2.405665     4.346421     3.220377     11.75906
# [6,]     3.149404     5.639050     3.199579     13.73917
# 
# apply(nonlin_norm, 2, function(x) c(mean = mean(x), sd = sd(x)))
# 
#      A Normalized B Normalized C Normalized D Normalized
# mean    2.3874296     5.549909    3.4000486     12.79232
# sd      0.4281159     1.378742    0.3191192      2.67708

## ----unequal, eval = FALSE----------------------------------------------------
# set.seed(123)
# vec1 <- rnorm(n = 10, mean = 0, sd = 1)
# vec2 <- rnorm(n = 5, mean = 5, sd = 5)
# vec3 <- rnorm(n = 8, mean = 10, sd = 10)
# 
# vec_list <- list(vec1, vec2, vec3)
# 
# NNS.norm(vec_list)
# 
# $`x_1 Normalized`
#  [1]  13.074058  -3.004912 -11.745878  25.406891  -4.647966  -5.481229   6.225165   5.920719   6.113733   9.640242
# 
# $`x_2 Normalized`
# [1]  2.875960212  0.008876158  1.230826150  5.855582361 10.779166523
# 
# $`x_3 Normalized`
# [1]  4.0749062  2.2395840  0.4067264  0.7457562 15.6445780  5.1941416  2.3326665  2.5622994

## ----rescale-minmax-----------------------------------------------------------
raw_vals <- c(-2.5, 0.2, 1.1, 3.7, 5.0)

scaled_minmax <- NNS.rescale(
  x = raw_vals,
  a = 5,
  b = 10,
  method = "minmax",
  T = NULL,
  type = "Terminal"
)

cbind(raw_vals, scaled_minmax)
range(scaled_minmax)

## ----rescale-riskneutral, eval=FALSE------------------------------------------
# set.seed(123)
# S0 <- 100
# r <- 0.05
# T <- 1
# 
# # Simulate a price path
# prices <- S0 * exp(cumsum(rnorm(250, 0.0005, 0.02)))
# 
# rn_terminal <- NNS.rescale(
#   x = prices,
#   a = S0,
#   b = r,
#   method = "riskneutral",
#   T = T,
#   type = "Terminal"
# )
# 
# c(
#   mean_original = mean(prices),
#   mean_rescaled = mean(rn_terminal),
#   target = S0 * exp(r * T)
# )
# 
# mean_original mean_rescaled        target
#      109.7019      105.1271      105.1271

## ----rescale-discounted, eval=FALSE-------------------------------------------
# rn_discounted <- NNS.rescale(
#   x = prices,
#   a = S0,
#   b = r,
#   method = "riskneutral",
#   T = T,
#   type = "Discounted"
# )
# 
# c(
#   mean_rescaled = mean(rn_discounted),
#   target_discounted_mean = S0
# )
# 
#          mean_rescaled target_discounted_mean
#                    100                    100

## ----image--------------------------------------------------------------------
set.seed(123)

x <- rnorm(1000, 5, 2)
y <- rgamma(1000, 3, 1)

# Combine variables
X <- cbind(x, y)

# NNS normalization
X_norm_lin <- NNS.norm(X, linear = TRUE)
X_norm_nonlin <- NNS.norm(X, linear = FALSE)

# Standard min-max normalization
minmax <- function(v) (v - min(v)) / (max(v) - min(v))
X_minmax <- apply(X, 2, minmax)

## ----plotting, echo=FALSE-----------------------------------------------------
par(mfrow = c(2,2))

steelblue_alpha <- rgb(1,0,0,0.4)
red_alpha <- rgb(0,0,1,0.4)

# Breaks for original data
br_orig <- pretty(range(c(x, y)), n = 15)

# Original variables
hist(x,
     col = steelblue_alpha,
     breaks = br_orig,
     main = "Original Variables",
     xlab = "")

hist(y,
     col = red_alpha,
     breaks = br_orig,
     add = TRUE)


# Breaks for NNS normalized variables
br_norm <- pretty(range(c(X_norm_lin[,1], X_norm_lin[,2])), n = 15)

# NNS normalized
hist(X_norm_lin[,1],
     col = steelblue_alpha,
     breaks = br_norm,
     main = "NNS.norm(..., Linear=TRUE)",
     xlab = "")

hist(X_norm_lin[,2],
     col = red_alpha,
     breaks = br_norm,
     add = TRUE)

# Breaks for NNS normalized variables
br_norm <- pretty(range(c(X_norm_nonlin[,1], X_norm_nonlin[,2])), n = 15)

# NNS normalized
hist(X_norm_nonlin[,1],
     col = steelblue_alpha,
     breaks = br_norm,
     main = "NNS.norm(..., Linear=FALSE)",
     xlab = "")

hist(X_norm_nonlin[,2],
     col = red_alpha,
     breaks = br_norm,
     add = TRUE)

# Breaks for min-max normalized variables
br_minmax <- pretty(range(c(X_minmax[,1], X_minmax[,2])), n = 15)

# Standard min-max normalization
hist(X_minmax[,1],
     col = steelblue_alpha,
     breaks = br_minmax,
     main = "Standard Min-Max",
     xlab = "")

hist(X_minmax[,2],
     col = red_alpha,
     breaks = br_minmax,
     add = TRUE)

