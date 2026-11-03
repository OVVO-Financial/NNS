# NNS.diff: derivative diagnostics and a reproducible scope

`NNS.diff` uses projected secant geometry to infer a perturbation scale and
returns several numerical derivative estimates with convergence diagnostics.
Its outputs answer different numerical questions. Complex-step differentiation
requires a function with an appropriate analytic complex extension. A real-only,
thresholded, or nonsmooth function does not satisfy that requirement simply
because it accepts a complex argument.

## An exact-function demonstration

```r
library(NNS)
f <- function(z) exp(z)
x0 <- 0.7
estimate <- NNS.diff(f, point = x0, h = 0.1, tol = 1e-10,
                     max.iter = 1000, print.trace = FALSE, plot = FALSE)
estimate
exp(x0)                       # known derivative for comparison
```

Inspect both termination diagnostics and error against known derivatives across
scales. A returned inferred step is a numerical choice, not a proof that the
search found a statistically optimal bandwidth.

## Why frozen real noise cannot establish complex-step robustness

```r
set.seed(20261103)
epsilon <- rnorm(1)
g <- function(z) exp(z) + epsilon
h <- 1e-12
c(clean = Im(f(x0 + 1i*h))/h,
  frozen_real_noise = Im(g(x0 + 1i*h))/h)
```

The real additive perturbation contributes zero imaginary part. The two outputs
therefore agree independently of any noise-adaptive property of NNS. A frozen
real noise surface evaluated only through `Re(z)` has the same limitation.
For a stochastic simulator, distinguish fresh independent noise per evaluation,
common random numbers, and a fixed perturbed function. Supply the complete noise
generator, seeds, target derivative, evaluation budget, and complex-input rule.
Compare repeated errors and failures on the same evaluation protocol.

The earlier numerical tables and universal ranking language on this page are
withdrawn because the noise generator needed to interpret them was absent.
They remain in Git history at commit
`1075a338d932db13a8210e8ed23fb94f70a3f21f`. The directional method remains available;
claims about noisy derivative accuracy require the stated reproducible evidence.
