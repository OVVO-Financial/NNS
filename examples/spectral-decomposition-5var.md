# 5-Variable NNS Spectral Decomposition

This example extends the directional spectral decomposition of Chapter 11 from the
bivariate case to five variables. The mean split — the same benchmark used throughout
the book — partitions each variable at its component mean, placing every observation
into one of 2⁵ = 32 orthants. All 32 orthants are occupied. The between-orthant and
within-orthant covariance matrices sum exactly to the population covariance matrix,
and the recovered eigensystem matches classical PCA to floating-point precision.

The second part of this example addresses a practical question the five-variable case
immediately raises: how does this scale to high-dimensional multivariate data? The
answer introduces the `DPM_nD` aggregation strategy and states its tradeoff precisely.

---

## Setup

Generate `n = 10,000` observations from a five-dimensional distribution with
uniform pairwise correlation `ρ ≈ 0.5`.

```r
set.seed(2024)
n <- 10000
d <- 5

R <- matrix(0.5, d, d) + diag(0.5, d)
L <- chol(R)
Z <- matrix(rnorm(n * d), n, d) %*% L
```

---

## Classical PCA

Compute the population-denominator covariance matrix and its eigensystem.

```r
mu    <- colMeans(Z)
Zc    <- sweep(Z, 2, mu)
Sigma <- crossprod(Zc) / n
pca   <- eigen(Sigma)

cat("Classical eigenvalues:\n")
print(round(pca$values, 6))
cat("\nClassical eigenvectors (first 2 cols):\n")
print(round(pca$vectors[, 1:2], 6))
```

```
Classical eigenvalues:
[1] 2.965193 0.510783 0.499329 0.493351 0.488601

Classical eigenvectors (first 2 cols):
          [,1]      [,2]
[1,] -0.451062  0.522912
[2,] -0.449908 -0.756151
[3,] -0.448047 -0.191485
[4,] -0.446032  0.329460
[5,] -0.440949  0.097919
```

The leading eigenvalue (2.965) is nearly six times the next largest (0.511),
consistent with a strong common factor driven by the uniform pairwise correlation.
The four residual eigenvalues are nearly equal (0.51, 0.50, 0.49, 0.49), indicating
approximately exchangeable residual structure across the minor axes.

---

## Mean-Split Orthant Partition

Partition each of the five variables at its component mean. Each observation falls
into one of up to 2⁵ = 32 orthants, encoded as a binary integer over which variables
are above their respective means.

```r
above_mean    <- sweep(Z, 2, mu, ">")
orthant_label <- apply(above_mean, 1, function(row) sum(row * 2^(0:(d - 1))))
n_orthants    <- length(unique(orthant_label))
cat("\nNumber of occupied orthants:", n_orthants, "out of", 2^d, "\n")
```

```
Number of occupied orthants: 32 out of 32
```

All 32 orthants are occupied. This is the direct multivariate generalization of the
four-quadrant CUPM/CLPM/DLPM/DUPM partition: each of the five sign combinations of
deviations from the component means receives positive empirical mass.

The mean split is the natural benchmark here because the between-orthant covariance
`Sigma_Q` is defined relative to the global mean `mu`. Using the mean as the split
point ensures that the displacement vectors `u_r = m_r - mu` measure deviations from
the same reference as the covariance matrix itself.

---

## Between-Within Decomposition

For each orthant `r`, compute the orthant probability `p_r`, the conditional mean
displacement `u_r = m_r - mu`, and the within-orthant covariance `\mathrm{Cov}_r`.
Accumulate the between-orthant and within-orthant covariance matrices.

```r
Sigma_Q <- matrix(0, d, d)
Sigma_W <- matrix(0, d, d)

for (lab in unique(orthant_label)) {
  mask  <- orthant_label == lab
  n_r   <- sum(mask)
  p_r   <- n_r / n

  Zr    <- Z[mask, , drop = FALSE]
  m_r   <- colMeans(Zr)
  u_r   <- m_r - mu
  Cov_r <- crossprod(sweep(Zr, 2, m_r)) / n_r

  Sigma_Q <- Sigma_Q + p_r * tcrossprod(u_r)
  Sigma_W <- Sigma_W + p_r * Cov_r
}

Sigma_rec <- Sigma_Q + Sigma_W
```

The decomposition identity is


```math
\Sigma
=
\underbrace{\sum_r p_r u_r u_r^\top}_{\Sigma_Q}
+
\underbrace{\sum_r p_r \mathrm{Cov}(Z \mid r)}_{\Sigma_W}.
```


Each of the 32 orthants contributes a rank-one spectral primitive
`B_r = p_r u_r u_r^T` to the between-orthant covariance, and a weighted
within-orthant scatter matrix to `Sigma_W`. The identity is the law of total
covariance applied to the mean-split orthant partition.

---

## Verification

```r
cat("\nMax absolute difference between original Σ and Σ_Q + Σ_W:",
    max(abs(Sigma - Sigma_rec)), "\n")
cat("Are they identical?", all.equal(Sigma, Sigma_rec, tolerance = 1e-12), "\n")
```

```
Max absolute difference between original Σ and Σ_Q + Σ_W: 3.330669e-15
Are they identical? TRUE
```

The recovery error is below `4e-15`, well within double-precision floating-point
tolerance. The identity holds exactly in finite arithmetic up to rounding.

---

## Eigensystem Recovery

```r
pca_rec <- eigen(Sigma_rec)

cat("\nRecovered eigenvalues:\n")
print(round(pca_rec$values, 6))
cat("Difference (classical - recovered):\n")
print(round(pca$values - pca_rec$values, 10))

for (j in 1:d) {
  if (sum(pca$vectors[, j] * pca_rec$vectors[, j]) < 0)
    pca_rec$vectors[, j] <- -pca_rec$vectors[, j]
}
cat("\nAbsolute difference in aligned eigenvectors (first 2):\n")
print(round(abs(pca$vectors[, 1:2] - pca_rec$vectors[, 1:2]), 10))
```

```
Recovered eigenvalues:
[1] 2.965193 0.510783 0.499329 0.493351 0.488601

Difference (classical - recovered):
[1] 0 0 0 0 0

Absolute difference in aligned eigenvectors (first 2):
     [,1] [,2]
[1,]    0    0
[2,]    0    0
[3,]    0    0
[4,]    0    0
[5,]    0    0
```

Eigenvalues and eigenvectors match to 10 decimal places after sign alignment.
The sign alignment step is cosmetic — eigenvectors are defined up to sign and
`eigen()` makes an arbitrary choice. After aligning, the absolute differences
are zero at every displayed precision.

---

## Eigenvalue Attribution

Each classical eigenvalue decomposes into a between-orthant contribution
(conditional mean displacement along the eigenvector) and a within-orthant
contribution (residual scatter projected onto the eigenvector):


```math
\lambda_i
=
\underbrace{v_i^\top \Sigma_Q v_i}_{\lambda_{i,Q}}
+
\underbrace{v_i^\top \Sigma_W v_i}_{\lambda_{i,W}}.
```


```r
attrib <- data.frame(eigenvalue = pca$values, between = NA, within = NA)
for (i in seq_along(pca$values)) {
  v                 <- pca$vectors[, i]
  attrib$between[i] <- drop(t(v) %*% Sigma_Q %*% v)
  attrib$within[i]  <- drop(t(v) %*% Sigma_W %*% v)
}
attrib$total       <- attrib$between + attrib$within
attrib$between_pct <- 100 * attrib$between / attrib$total

cat("\nEigenvalue attribution:\n")
print(round(attrib, 6))
```

```
Eigenvalue attribution:
  eigenvalue  between   within    total between_pct
1   2.965193 2.439089 0.526104 2.965193    82.25735
2   0.510783 0.243905 0.266878 0.510783    47.75114
3   0.499329 0.241130 0.258199 0.499329    48.29081
4   0.493351 0.237770 0.255581 0.493351    48.19483
5   0.488601 0.235968 0.252633 0.488601    48.29467
```

**PC1 is 82.3% between-orthant.** The leading eigenvalue is dominated by
separation among the 32 orthant conditional means. The mean-split orthant
partition captures the relevant second-moment geometry of the common factor
almost entirely through conditional mean displacement; within-orthant residual
scatter contributes only 17.7%.

**PC2–PC5 split approximately evenly.** Each of the four residual eigenvectors
receives roughly 48% from between-orthant and 52% from within-orthant scatter.
This is consistent with the near-equal residual eigenvalues: with no dominant
directional structure in the residuals, between and within contributions balance.

The diagnostic ratio for PC1 is


```math
D_{\mathrm{spectral}}
=
\frac{v_1^\top \Sigma_Q v_1}{\lambda_1}
=
\frac{2.4391}{2.9652}
\approx 0.823.
```


A value this high confirms that the 32-orthant mean-split partition is capturing
the dominant covariance geometry almost entirely through conditional mean structure.
The within-orthant residual is small because the common factor aligns tightly with
the direction of maximum conditional mean separation across orthants.

---

## Orthant-Level Attribution of PC1

The between-orthant contribution to `lambda_1` decomposes further into per-orthant
rank-one terms:


```math
\lambda_{1,Q}
=
\sum_r p_r (v_1^\top u_r)^2.
```


```r
v1              <- pca$vectors[, 1]
orthant_contrib <- data.frame(orthant      = unique(orthant_label),
                              contribution = NA_real_)
for (i in seq_along(orthant_contrib$orthant)) {
  lab  <- orthant_contrib$orthant[i]
  mask <- orthant_label == lab
  p_r  <- sum(mask) / n
  u_r  <- colMeans(Z[mask, , drop = FALSE]) - mu
  orthant_contrib$contribution[i] <- p_r * (sum(v1 * u_r))^2
}

cat("\nSum of orthant-level between contributions for PC1:",
    sum(orthant_contrib$contribution), "\n")
cat("Direct Σ_Q between contribution for PC1:",
    attrib$between[1], "\n")
```

```
Sum of orthant-level between contributions for PC1: 2.439089
Direct Σ_Q between contribution for PC1: 2.439089
```

The per-orthant terms sum to exactly the Rayleigh quotient `v1^T Sigma_Q v1`.
This confirms that the between-orthant attribution is not an approximation — it is
an exact partition of `lambda_{1,Q}` into 32 identifiable orthant contributions,
each equal to `p_r (v1^T u_r)^2`.

The largest contributions come from the two fully concordant orthants: all five
variables above their means (binary label 31) and all five below (binary label 0).
These two orthants have the largest displacement magnitudes `||u_r||` and align most
strongly with the common-factor eigenvector `v_1`, which loads approximately equally
on all five variables.

---

## Converse Failure

```r
cat("\nGiven only the eigenvalues/eigenvectors, can we recover the orthant means?\n")
cat("No. Example: conditional mean of most populated orthant is:\n")
lab_max <- names(which.max(table(orthant_label)))
print(round(colMeans(Z[orthant_label == lab_max, , drop = FALSE]), 6))
cat("This cannot be deduced from PCA output.\n")
```

```
Given only the eigenvalues/eigenvectors, can we recover the orthant means?
No. Example: conditional mean of most populated orthant is:
[1] 1.092094 1.088237 1.049830 1.084412 1.071189
This cannot be deduced from PCA output.
```

The most populated orthant — all five variables above their component means — has a
conditional mean approximately 1.07 standard deviations above `mu` in every
dimension. This is exactly the all-concordant-upper orthant, the five-dimensional
analogue of CUPM. Classical PCA reports only `(lambda_i, v_i)`. Neither the orthant
assignment, the orthant probabilities, nor the 32 conditional means appear anywhere
in that output. The directional decomposition runs in one direction:


```math
\{p_r, m_r, \mathrm{Cov}(Z \mid r)\}_r
\;\longrightarrow\;
\Sigma
\;\longrightarrow\;
(\lambda_i, v_i).
```


The map does not reverse.

---


## Scaling to Higher Dimensions: The Curse of Dimensionality

The five-variable example is a useful proof-of-concept. It confirms that the full
mean-split orthant decomposition works cleanly for:

```math
2^5 = 32
```

orthants and recovers the PCA eigensystem to floating-point precision.

A natural statistical question follows:

> How does this scale to 50 variables?

The full orthant decomposition scales exponentially. For a 50-dimensional dataset,
the number of possible orthants is:

```math
2^{50}
\approx
1.13 \times 10^{15}.
```

That is not a practical partition to estimate directly. Even if the computation were
possible, most orthants would be empty or too sparsely populated to support stable
conditional mean and covariance estimates. This is the curse of dimensionality in
its most direct form.

### State Aggregation with DPM_nD

The n-dimensional partial moment framework provides a practical aggregation strategy.
Rather than retaining all orthants separately, the high-dimensional state space can
be collapsed into three macroscopic, observable directional states:

```math
\mathrm{CLPM}_{nD}
=
\mathrm{all\ variables\ simultaneously\ below\ target},
```

```math
\mathrm{CUPM}_{nD}
=
\mathrm{all\ variables\ simultaneously\ above\ target},
```

```math
\mathrm{DPM}_{nD}
=
\mathrm{all\ mixed\ sign\ configurations}.
```

The state count is therefore reduced from:

```math
2^d
```

to:

```math
3.
```

This is a dimensionality reduction by aggregation.

### General Statistical Interpretation

The three aggregate states have direct statistical meanings:

- `CLPM_nD` measures joint lower-tail concordance.
- `CUPM_nD` measures joint upper-tail concordance.
- `DPM_nD` measures mixed-sign divergence or dispersion.

The normalized quantities are shares of total directional mass:

```math
\mathrm{CLPM}^{\,\mathrm{norm}}_{nD}
=
\frac{\mathrm{CLPM}_{nD}}
{\mathrm{CLPM}_{nD} + \mathrm{CUPM}_{nD} + \mathrm{DPM}_{nD}},
```

```math
\mathrm{CUPM}^{\,\mathrm{norm}}_{nD}
=
\frac{\mathrm{CUPM}_{nD}}
{\mathrm{CLPM}_{nD} + \mathrm{CUPM}_{nD} + \mathrm{DPM}_{nD}},
```

```math
\mathrm{DPM}^{\,\mathrm{norm}}_{nD}
=
\frac{\mathrm{DPM}_{nD}}
{\mathrm{CLPM}_{nD} + \mathrm{CUPM}_{nD} + \mathrm{DPM}_{nD}}.
```

At degree zero, these are probability shares. At higher degrees, they are
severity-weighted directional mass shares.

This gives a compact high-dimensional summary of whether observations tend to
cluster in all-lower, all-upper, or mixed-sign regions relative to a target vector.

### The Tradeoff

The aggregation is useful, but it comes with a clear tradeoff.

`DPM_nD` preserves computational tractability and captures the two fully concordant
tail states: all variables below target and all variables above target. These are
often important summary states in multivariate dependence analysis.

However, `DPM_nD` sacrifices granular geometric resolution inside the mixed-sign
region. It does not identify exactly which variables diverged from which others. For
example, in a five-variable setting, all of the following are mixed-sign states:

```text
Variable1+ Variable2+ Variable3- Variable4- Variable5-
Variable1+ Variable2- Variable3+ Variable4- Variable5+
Variable1- Variable2+ Variable3- Variable4+ Variable5+
```

The full orthant-level decomposition keeps these states separate. `DPM_nD` aggregates
them.

Therefore, the full spectral genealogy requires:

```math
\{p_r,\;u_r,\;\Sigma_r\}_{r=1}^{2^d},
```

while the scalable directional summary uses:

```math
\mathrm{CLPM}_{nD},
\qquad
\mathrm{CUPM}_{nD},
\qquad
\mathrm{DPM}_{nD}.
```

The full orthant-level between covariance is:

```math
\Sigma_Q
=
\sum_{r=1}^{2^d} p_r u_r u_r^\top.
```

The `DPM_nD` aggregation does not preserve each individual mixed-orthant displacement
vector `u_r`. It is therefore best understood as a scalable descriptive statistic,
not as a replacement for exact orthant-level spectral attribution.

### General Workflow

The two approaches are complementary.

A useful high-dimensional workflow is:

1. Use `CLPM_nD`, `CUPM_nD`, and `DPM_nD` as compact directional summaries.
2. Use the full orthant decomposition when `d` is small enough for reliable
   estimation.
3. Use grouped orthants, selected orthants, or lower-dimensional factor partitions
   when intermediate resolution is needed.
4. Reserve exact per-orthant spectral attribution for settings where the number of
   occupied orthants is statistically manageable.

In short:

```math
\mathrm{full\ orthant\ decomposition}
=
\mathrm{exact\ but\ exponential}.
```

```math
\mathrm{DPM}_{nD}\ \mathrm{aggregation}
=
\mathrm{scalable\ but\ coarser}.
```

This section bridges the gap between theoretical exactness and high-dimensional
statistical practice. The five-variable decomposition proves the identity. The
`DPM_nD` aggregation explains how related directional information can still be
summarized when full orthant enumeration is not feasible.

---

## Full R Code

```r
# =============================================================================
# 5-Variable NNS Spectral Decomposition (mean split)
# =============================================================================
set.seed(2024)
n <- 10000
d <- 5

# Generate correlated data
R <- matrix(0.5, d, d) + diag(0.5, d)
L <- chol(R)
Z <- matrix(rnorm(n * d), n, d) %*% L

# ---------------------------------------------------------------------------
# Classical PCA (population-denominator covariance)
# ---------------------------------------------------------------------------
mu    <- colMeans(Z)
Zc    <- sweep(Z, 2, mu)
Sigma <- crossprod(Zc) / n
pca   <- eigen(Sigma)

cat("Classical eigenvalues:\n")
print(round(pca$values, 6))
cat("\nClassical eigenvectors (first 2 cols):\n")
print(round(pca$vectors[, 1:2], 6))

# ---------------------------------------------------------------------------
# Mean-split orthant partition
# ---------------------------------------------------------------------------
above_mean    <- sweep(Z, 2, mu, ">")
orthant_label <- apply(above_mean, 1, function(row) sum(row * 2^(0:(d - 1))))
n_orthants    <- length(unique(orthant_label))
cat("\nNumber of occupied orthants:", n_orthants, "out of", 2^d, "\n")

# ---------------------------------------------------------------------------
# Between-orthant Sigma_Q and within-orthant Sigma_W
# ---------------------------------------------------------------------------
Sigma_Q <- matrix(0, d, d)
Sigma_W <- matrix(0, d, d)

for (lab in unique(orthant_label)) {
  mask  <- orthant_label == lab
  n_r   <- sum(mask)
  p_r   <- n_r / n

  Zr    <- Z[mask, , drop = FALSE]
  m_r   <- colMeans(Zr)
  u_r   <- m_r - mu
  Cov_r <- crossprod(sweep(Zr, 2, m_r)) / n_r

  Sigma_Q <- Sigma_Q + p_r * tcrossprod(u_r)
  Sigma_W <- Sigma_W + p_r * Cov_r
}

Sigma_rec <- Sigma_Q + Sigma_W

# ---------------------------------------------------------------------------
# Verification
# ---------------------------------------------------------------------------
cat("\nMax absolute difference between original Sigma and Sigma_Q + Sigma_W:",
    max(abs(Sigma - Sigma_rec)), "\n")
cat("Are they identical?", all.equal(Sigma, Sigma_rec, tolerance = 1e-12), "\n")

# ---------------------------------------------------------------------------
# Recover eigensystem
# ---------------------------------------------------------------------------
pca_rec <- eigen(Sigma_rec)
cat("\nRecovered eigenvalues:\n")
print(round(pca_rec$values, 6))
cat("Difference (classical - recovered):\n")
print(round(pca$values - pca_rec$values, 10))

for (j in 1:d) {
  if (sum(pca$vectors[, j] * pca_rec$vectors[, j]) < 0)
    pca_rec$vectors[, j] <- -pca_rec$vectors[, j]
}
cat("\nAbsolute difference in aligned eigenvectors (first 2):\n")
print(round(abs(pca$vectors[, 1:2] - pca_rec$vectors[, 1:2]), 10))

# ---------------------------------------------------------------------------
# Eigenvalue attribution
# ---------------------------------------------------------------------------
attrib <- data.frame(eigenvalue = pca$values, between = NA, within = NA)
for (i in seq_along(pca$values)) {
  v                 <- pca$vectors[, i]
  attrib$between[i] <- drop(t(v) %*% Sigma_Q %*% v)
  attrib$within[i]  <- drop(t(v) %*% Sigma_W %*% v)
}
attrib$total       <- attrib$between + attrib$within
attrib$between_pct <- 100 * attrib$between / attrib$total

cat("\nEigenvalue attribution:\n")
print(round(attrib, 6))

# ---------------------------------------------------------------------------
# Orthant-level contribution to PC1 between part
# ---------------------------------------------------------------------------
v1              <- pca$vectors[, 1]
orthant_contrib <- data.frame(orthant      = unique(orthant_label),
                              contribution = NA_real_)
for (i in seq_along(orthant_contrib$orthant)) {
  lab  <- orthant_contrib$orthant[i]
  mask <- orthant_label == lab
  p_r  <- sum(mask) / n
  u_r  <- colMeans(Z[mask, , drop = FALSE]) - mu
  orthant_contrib$contribution[i] <- p_r * (sum(v1 * u_r))^2
}

cat("\nSum of orthant-level between contributions for PC1:",
    sum(orthant_contrib$contribution), "\n")
cat("Direct Sigma_Q between contribution for PC1:",
    attrib$between[1], "\n")

# ---------------------------------------------------------------------------
# Converse failure illustration
# ---------------------------------------------------------------------------
cat("\nGiven only the eigenvalues/eigenvectors, can we recover the orthant means?\n")
cat("No. Example: conditional mean of most populated orthant is:\n")
lab_max <- names(which.max(table(orthant_label)))
print(round(colMeans(Z[orthant_label == lab_max, , drop = FALSE]), 6))
cat("This cannot be deduced from PCA output.\n")
```
