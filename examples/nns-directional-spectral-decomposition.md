# NNS Directional Spectral Decomposition: Full Orthants, Pairwise Matrices, and Scalable PCA Attribution

This note extends the Chapter 11 directional spectral decomposition workflow in three layers:

1. **Full orthant decomposition**  
   The complete mean-split partition gives an exact spectral genealogy of covariance and PCA, but it scales as $2^d$.

2. **Pairwise partial-moment matrices**  
   Pairwise directional matrices recover covariance and PCA exactly while scaling as $O(d^2)$. This is the practical high-dimensional PCA attribution layer.

3. **`DPM_nD` aggregation**  
   The n-dimensional NNS functions collapse the full orthant partition into three global states: all-lower, all-upper, and mixed-sign. This is scalable and useful, but it is a coarse directional summary, not a full spectral genealogy.

The main conclusion is:

> Full orthants explain PCA completely.  
> Pairwise partial-moment matrices scale PCA recovery and directional attribution to high dimensions.  
> `DPM_nD` gives a compact global diagnostic, not a replacement for spectral decomposition.

---

## Executive Summary

Classical PCA starts with a covariance matrix:

```math
\Sigma,
```

and diagonalizes it:

```math
\Sigma v_i = \lambda_i v_i.
```

PCA reports eigenvalues and eigenvectors, but it does not explain which directional regions of the data generated them.

The NNS directional framework begins with observable directional structure. In the full orthant version, every observation is assigned to a mean-split orthant. For $d$ variables, this gives:

```math
2^d
```

possible states.

Each occupied orthant $r$ has:

- empirical probability $p_r$;
- conditional mean $m_r$;
- centered displacement $u_r = m_r-\mu$;
- within-orthant covariance $\Sigma_r$.

The covariance matrix decomposes exactly as:

```math
\Sigma
=
\Sigma_Q+\Sigma_W,
```

where:

```math
\Sigma_Q
=
\sum_r p_r u_r u_r^\top,
```

and:

```math
\Sigma_W
=
\sum_r p_r\Sigma_r.
```

Thus:

```math
\Sigma
=
\sum_r p_r u_r u_r^\top
+
\sum_r p_r\Sigma_r.
```

Diagonalizing the recovered matrix gives the same classical PCA eigensystem.

The full orthant decomposition is exact, but it becomes infeasible in high dimensions. For $d=50$,

```math
2^{50}
\approx
1.1259\times 10^{15}.
```

The practical high-dimensional alternative is to use pairwise partial-moment matrices. These recover covariance through:

```math
\Sigma
=
\mathrm{CLPM}
+
\mathrm{CUPM}
-
\mathrm{DLPM}
-
\mathrm{DUPM}.
```

Then PCA is recovered from:

```math
(\mathrm{CLPM},\mathrm{CUPM},\mathrm{DLPM},\mathrm{DUPM})
\longrightarrow
\Sigma
\longrightarrow
(\lambda_i,v_i).
```

This preserves exact covariance and PCA recovery while scaling with four $d\times d$ matrices rather than $2^d$ orthant states.

---

# Part I — Full Orthant Decomposition

## 1. Setup

Generate $n=10{,}000$ observations from a five-dimensional distribution with approximately uniform pairwise correlation $0.5$.

```r
set.seed(2024)

n <- 10000
d <- 5

R <- matrix(0.5, d, d) + diag(0.5, d)
L <- chol(R)

Z <- matrix(rnorm(n * d), n, d) %*% L
colnames(Z) <- paste0("X", seq_len(d))

target <- colMeans(Z)
```

The component mean vector is used as the target.

In the reported run, the target vector was:

| Variable | Target |
|---|---:|
| X1 | -0.006349 |
| X2 | -0.004824 |
| X3 | -0.002800 |
| X4 | 0.004986 |
| X5 | 0.001181 |

---

## 2. Classical PCA

Compute the population-denominator covariance matrix:

```r
mu    <- colMeans(Z)
Zc    <- sweep(Z, 2, mu)
Sigma <- crossprod(Zc) / n
pca   <- eigen(Sigma, symmetric = TRUE)
```

The classical eigenvalues were:

| Component | Eigenvalue |
|---:|---:|
| PC1 | 2.965193 |
| PC2 | 0.510783 |
| PC3 | 0.499329 |
| PC4 | 0.493351 |
| PC5 | 0.488601 |

The first eigenvalue is much larger than the remaining four, consistent with a strong common direction induced by the uniform correlation structure.

---

## 3. Mean-Split Orthant Partition

Partition each variable at its component mean. For five variables, the full mean-split partition contains:

```math
2^5 = 32
```

possible orthants.

Each observation is encoded according to whether each variable is above or below its component mean:

```r
above_mean <- sweep(Z, 2, mu, ">")

orthant_label <- apply(above_mean, 1, function(row) {
  sum(row * 2^(0:(d - 1)))
})
```

The run occupied all orthants:

```text
32 out of 32
```

This is the five-dimensional analogue of the bivariate quadrant partition.

---

## 4. Between-Within Orthant Decomposition

For each orthant $r$, compute:

- empirical probability $p_r$;
- conditional mean $m_r$;
- centered displacement $u_r = m_r-\mu$;
- within-orthant covariance $\Sigma_r$.

Then accumulate:

```math
\Sigma_Q
=
\sum_{r=1}^{2^d}p_r u_r u_r^\top,
```

and:

```math
\Sigma_W
=
\sum_{r=1}^{2^d}p_r\Sigma_r.
```

The full covariance identity is:

```math
\Sigma
=
\Sigma_Q+\Sigma_W.
```

In code:

```r
Sigma_Q <- matrix(0, d, d)
Sigma_W <- matrix(0, d, d)

for (lab in unique(orthant_label)) {
  mask <- orthant_label == lab
  n_r  <- sum(mask)
  p_r  <- n_r / n

  Zr      <- Z[mask, , drop = FALSE]
  m_r     <- colMeans(Zr)
  u_r     <- m_r - mu
  Sigma_r <- crossprod(sweep(Zr, 2, m_r)) / n_r

  Sigma_Q <- Sigma_Q + p_r * tcrossprod(u_r)
  Sigma_W <- Sigma_W + p_r * Sigma_r
}

Sigma_rec <- Sigma_Q + Sigma_W
```

Each orthant contributes a rank-one between-orthant primitive:

```math
B_r
=
p_r u_r u_r^\top.
```

If $u_r\neq 0$, then:

```math
B_r u_r
=
p_r\|u_r\|^2u_r.
```

So each centered orthant conditional mean is the nonzero eigenvector of its own rank-one between-orthant contribution.

---

## 5. Covariance and PCA Recovery

The full orthant decomposition recovered the covariance matrix to floating-point precision:

```text
Max absolute covariance recovery error:
3.330669e-15
```

The recovered eigensystem also matched the classical eigensystem:

```text
Full orthant eigenvalue recovery error:
3.663736e-15

Full orthant eigenvector alignments:
1 1 1 1 1
```

Thus:

```math
\{p_r,m_r,\Sigma_r\}_{r=1}^{2^d}
\longrightarrow
\Sigma_Q+\Sigma_W
\longrightarrow
\Sigma
\longrightarrow
(\lambda_i,v_i).
```

The full orthant decomposition recovers the covariance matrix and therefore recovers the classical PCA eigensystem.

---

## 6. Eigenvalue Attribution from Full Orthants

Each classical eigenvalue decomposes into a between-orthant contribution and a within-orthant contribution.

For unit eigenvector $v_i$:

```math
\lambda_i
=
v_i^\top\Sigma_Qv_i
+
v_i^\top\Sigma_Wv_i.
```

Using the orthant decomposition:

```math
\lambda_i
=
\sum_{r=1}^{2^d}p_r(v_i^\top u_r)^2
+
\sum_{r=1}^{2^d}p_rv_i^\top\Sigma_rv_i.
```

The attribution table was:

| Component | Eigenvalue | Between | Within | Total | Between Percent |
|---:|---:|---:|---:|---:|---:|
| PC1 | 2.965193 | 2.439089 | 0.526104 | 2.965193 | 82.25735 |
| PC2 | 0.510783 | 0.243905 | 0.266878 | 0.510783 | 47.75114 |
| PC3 | 0.499329 | 0.241130 | 0.258199 | 0.499329 | 48.29081 |
| PC4 | 0.493351 | 0.237770 | 0.255581 | 0.493351 | 48.19483 |
| PC5 | 0.488601 | 0.235968 | 0.252633 | 0.488601 | 48.29467 |

For PC1:

```math
\frac{2.439089}{2.965193}
=
0.8225735.
```

So about 82.3 percent of PC1 came from between-orthant conditional mean separation.

A useful global summary is:

```math
D_{\mathrm{spectral}}
=
\frac{\mathrm{tr}(\Sigma_Q)}{\mathrm{tr}(\Sigma)}.
```

In this run:

```text
D_spectral = 0.685432
```

So about 68.5 percent of total variance came from between-orthant conditional mean displacement.

---

## 7. Orthant-Level Attribution of PC1

The between-orthant part of PC1 decomposes further into individual orthant contributions:

```math
\lambda_{1,Q}
=
\sum_{r=1}^{2^d}p_r(v_1^\top u_r)^2.
```

The top ten orthant contributions to PC1 were:

| Rank | Orthant | Pattern | Probability | Contribution | Percent of PC1 |
|---:|---:|---|---:|---:|---:|
| 1 | 31 | X1+ X2+ X3+ X4+ X5+ | 0.1663 | 0.9676085 | 32.63223 |
| 2 | 0 | X1- X2- X3- X4- X5- | 0.1646 | 0.9344073 | 31.51253 |
| 3 | 16 | X1- X2- X3- X4- X5+ | 0.0328 | 0.0539840 | 1.82059 |
| 4 | 4 | X1- X2- X3+ X4- X5- | 0.0342 | 0.0517459 | 1.74511 |
| 5 | 8 | X1- X2- X3- X4+ X5- | 0.0334 | 0.0513473 | 1.73167 |
| 6 | 2 | X1- X2+ X3- X4- X5- | 0.0321 | 0.0504620 | 1.70181 |
| 7 | 23 | X1+ X2+ X3+ X4- X5+ | 0.0332 | 0.0498662 | 1.68172 |
| 8 | 30 | X1- X2+ X3+ X4+ X5+ | 0.0340 | 0.0480060 | 1.61899 |
| 9 | 27 | X1+ X2+ X3- X4+ X5+ | 0.0330 | 0.0475473 | 1.60351 |
| 10 | 29 | X1+ X2- X3+ X4+ X5+ | 0.0311 | 0.0475126 | 1.60234 |

The top two orthants were the all-upper and all-lower states:

```text
X1+ X2+ X3+ X4+ X5+
X1- X2- X3- X4- X5-
```

Together they contributed:

```math
32.63223\% + 31.51253\%
=
64.14476\%
```

of PC1.

The orthant-level PC1 contributions summed exactly to the direct Rayleigh quotient:

```text
Sum of orthant-level between contributions for PC1:
2.439089

Direct Sigma_Q between contribution for PC1:
2.439089
```

This is the full orthant-level spectral genealogy.

---

## 8. Converse Failure

The decomposition runs in one direction.

From the full orthant decomposition, one recovers the covariance matrix:

```math
\{p_r,m_r,\Sigma_r\}_{r=1}^{2^d}
\longrightarrow
\Sigma.
```

From the covariance matrix, one recovers the PCA eigensystem:

```math
\Sigma
\longrightarrow
(\lambda_i,v_i).
```

But PCA output alone does not recover the orthant probabilities, orthant assignments, or orthant conditional means.

Therefore:

```math
\{p_r,m_r,\Sigma_r\}_{r=1}^{2^d}
\Rightarrow
\Sigma
\Rightarrow
(\lambda_i,v_i),
```

but generally:

```math
(\lambda_i,v_i)
\not\Rightarrow
\{p_r,m_r,\Sigma_r\}_{r=1}^{2^d}.
```

PCA is a downstream summary. The orthant decomposition contains strictly more directional information.

---

# Part II — `DPM_nD` as a Three-State Global Diagnostic

## 9. Why `DPM_nD` Is Useful but Coarse

The full orthant decomposition is exact, but it scales as:

```math
2^d.
```

For five variables:

```math
2^5 = 32.
```

For 50 variables:

```math
2^{50}
\approx
1.1259\times10^{15}.
```

Most high-dimensional orthants would be empty or too sparsely populated to support stable conditional mean and covariance estimates.

The NNS n-dimensional partial moment functions collapse the full orthant partition into three observable aggregate states:

```math
\mathrm{CLPM}_{nD}
=
\mathrm{all\ variables\ below\ target},
```

```math
\mathrm{CUPM}_{nD}
=
\mathrm{all\ variables\ above\ target},
```

```math
\mathrm{DPM}_{nD}
=
\mathrm{all\ mixed\ sign\ configurations}.
```

So the state count is reduced from:

```math
2^d
```

to:

```math
3.
```

For $d=5$, the full partition has 32 orthants. The three-state aggregation is:

```math
G_L=\{0\},
```

```math
G_U=\{31\},
```

```math
G_D=\{1,2,\ldots,30\}.
```

Thus:

```math
\mathrm{CLPM}_{nD}
=
\sum_{r\in G_L}p_r,
```

```math
\mathrm{CUPM}_{nD}
=
\sum_{r\in G_U}p_r,
```

```math
\mathrm{DPM}_{nD}
=
\sum_{r\in G_D}p_r.
```

At degree zero, `DPM_nD` is not separate from the full orthant decomposition. It is the coarse three-state projection of it.

---

## 10. `DPM_nD` R Example Using the NNS Package

The exported NNS package functions are:

```r
NNS::Co.LPM_nD(data, target, degree = 0, norm = TRUE)
NNS::Co.UPM_nD(data, target, degree = 0, norm = TRUE)
NNS::DPM_nD(data, target, degree = 0, norm = TRUE)
```

The direct full-orthant aggregation gave:

| State | Full Orthant Aggregation |
|---|---:|
| `CLPM_nD` | 0.1646 |
| `CUPM_nD` | 0.1663 |
| `DPM_nD` | 0.6691 |

The NNS degree-zero values were identical:

| State | NNS Degree-Zero Value |
|---|---:|
| `CLPM_nD` | 0.1646 |
| `CUPM_nD` | 0.1663 |
| `DPM_nD` | 0.6691 |

The difference between the NNS output and the full orthant aggregation was:

```text
0 0 0
```

This proves that, at degree zero, the NNS nD functions exactly recover the three-state aggregation of the full orthant partition.

---

## 11. Degree-One `DPM_nD` Directional Mass

The raw degree-one directional masses were:

| State | Raw Degree-One Mass |
|---|---:|
| `CLPM_nD` | 0.605363 |
| `CUPM_nD` | 0.575106 |
| `DPM_nD` | 0.086607 |

The normalized degree-one directional mass shares were:

| State | Normalized Share |
|---|---:|
| `CLPM_nD` | 0.477763 |
| `CUPM_nD` | 0.453885 |
| `DPM_nD` | 0.068352 |

Using `norm = TRUE` returned the same normalized shares:

| State | `norm = TRUE` |
|---|---:|
| `CLPM_nD` | 0.47776331 |
| `CUPM_nD` | 0.45388455 |
| `DPM_nD` | 0.06835214 |

The interpretation is:

> Mixed-sign observations are frequent, but they contribute relatively little first-degree directional mass in this correlated five-variable example.

This is consistent with the PCA attribution result: the leading covariance direction is dominated by concordant all-lower and all-upper conditional mean separation.

---

## 12. Limitation of `DPM_nD` for Spectral Genealogy

The full spectral decomposition requires more than probabilities. It requires each orthant displacement vector and each within-orthant covariance:

```math
\{p_r,u_r,\Sigma_r\}_{r=1}^{2^d}.
```

The full between-orthant covariance is:

```math
\Sigma_Q
=
\sum_{r=1}^{2^d}p_r u_r u_r^\top.
```

When the mixed orthants are collapsed into a single `DPM_nD` value, their individual displacement vectors are no longer separately retained.

Therefore:

```math
\mathrm{full\ orthant\ decomposition}
\Rightarrow
(\mathrm{CLPM}_{nD},\mathrm{CUPM}_{nD},\mathrm{DPM}_{nD}),
```

but:

```math
(\mathrm{CLPM}_{nD},\mathrm{CUPM}_{nD},\mathrm{DPM}_{nD})
\not\Rightarrow
\mathrm{full\ orthant\ decomposition}.
```

So the correct role of `DPM_nD` is:

```math
\mathrm{DPM}_{nD}
=
\mathrm{scalable\ three\ state\ global\ directional\ diagnostic}.
```

It is not the main high-dimensional PCA recovery method.

---

# Part III — Pairwise Partial-Moment Matrices as the Scalable PCA Method

## 13. Pairwise Directional Matrix Identity

The scalable alternative is to use pairwise directional partial-moment matrices.

For a $d$-variable matrix $Z$, compute the four $d\times d$ matrices:

```math
\mathrm{CLPM},
\qquad
\mathrm{CUPM},
\qquad
\mathrm{DLPM},
\qquad
\mathrm{DUPM}.
```

These reassemble covariance as:

```math
\Sigma
=
\mathrm{CLPM}
+
\mathrm{CUPM}
-
\mathrm{DLPM}
-
\mathrm{DUPM}.
```

In R:

```r
pm <- NNS::PM.matrix(
  LPM_degree = 1,
  UPM_degree = 1,
  target     = "mean",
  variable   = Z,
  pop_adj    = TRUE,
  norm       = FALSE
)

CLPM <- pm$clpm
CUPM <- pm$cupm
DLPM <- pm$dlpm
DUPM <- pm$dupm

Sigma_nns <- CLPM + CUPM - DLPM - DUPM
```

Then PCA is recovered by:

```r
pca_nns <- eigen(Sigma_nns, symmetric = TRUE)
```

The recovery chain is:

```math
(\mathrm{CLPM},\mathrm{CUPM},\mathrm{DLPM},\mathrm{DUPM})
\longrightarrow
\Sigma
\longrightarrow
(\lambda_i,v_i).
```

This uses four $d\times d$ matrices rather than $2^d$ orthants.

---

## 14. Pairwise Eigenvalue Attribution

For any classical PCA direction $v_i$,

```math
\lambda_i
=
v_i^\top \Sigma v_i.
```

Substituting the pairwise directional decomposition gives:

```math
\lambda_i
=
v_i^\top\mathrm{CLPM}v_i
+
v_i^\top\mathrm{CUPM}v_i
-
v_i^\top\mathrm{DLPM}v_i
-
v_i^\top\mathrm{DUPM}v_i.
```

This identifies whether a principal component is generated by:

- lower concordance;
- upper concordance;
- lower-divergent movement;
- upper-divergent movement.

This is not a full orthant genealogy. It does not recover complete simultaneous states such as:

```text
X1+ X2- X3+ ... X50-
```

But it does recover exact covariance and PCA while preserving interpretable pairwise directional sources.

---

## 15. Five-Variable Pairwise Recovery

Using the same five-variable dataset, pairwise `PM.matrix` recovered covariance and PCA to numerical precision.

```text
Pairwise PM.matrix covariance recovery error:
1.776357e-15

PM.matrix $cov.matrix vs manual reassembly error:
1.110223e-16

Directional aggregate asymmetry before symmetrization:
1.110223e-16
```

The classical sample PCA eigenvalues were:

| Component | Eigenvalue |
|---:|---:|
| PC1 | 2.965489 |
| PC2 | 0.510834 |
| PC3 | 0.499379 |
| PC4 | 0.493400 |
| PC5 | 0.488650 |

The pairwise NNS recovered eigenvalues were identical at displayed precision:

| Component | Recovered Eigenvalue |
|---:|---:|
| PC1 | 2.965489 |
| PC2 | 0.510834 |
| PC3 | 0.499379 |
| PC4 | 0.493400 |
| PC5 | 0.488650 |

The eigenvalue recovery error was:

```text
2.720046e-15
```

The eigenvector alignments were:

```text
1 1 1 1 1
```

---

## 16. Five-Variable Pairwise Directional Attribution

The pairwise attribution table was:

| PC | Eigenvalue | CLPM | CUPM | DLPM | DUPM | Recovered |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 2.965489 | 1.693504 | 1.708038 | 0.218026 | 0.218026 | 2.965489 |
| 2 | 0.510834 | 0.200757 | 0.199162 | -0.055458 | -0.055458 | 0.510834 |
| 3 | 0.499379 | 0.198809 | 0.190119 | -0.055225 | -0.055225 | 0.499379 |
| 4 | 0.493400 | 0.196324 | 0.188517 | -0.054280 | -0.054280 | 0.493400 |
| 5 | 0.488650 | 0.190190 | 0.192332 | -0.053064 | -0.053064 | 0.488650 |

The signed percentage attribution was:

| PC | CLPM % | CUPM % | DLPM signed % | DUPM signed % |
|---:|---:|---:|---:|---:|
| 1 | 57.11 | 57.60 | -7.35 | -7.35 |
| 2 | 39.30 | 38.99 | 10.86 | 10.86 |
| 3 | 39.81 | 38.07 | 11.06 | 11.06 |
| 4 | 39.79 | 38.21 | 11.00 | 11.00 |
| 5 | 38.92 | 39.36 | 10.86 | 10.86 |

For PC1:

```math
2.965489
=
1.693504
+
1.708038
-
0.218026
-
0.218026.
```

So PC1 is dominated by pairwise lower and upper concordance, with divergent pairwise movement subtracting from the common direction.

At the trace level:

| Component | Trace Contribution | Trace Percent |
|---|---:|---:|
| CLPM | 2.479584 | 50.01 |
| CUPM | 2.478168 | 49.99 |
| DLPM signed | 0.000000 | 0.00 |
| DUPM signed | 0.000000 | 0.00 |
| Total | 4.957752 | 100.00 |

Total variance is evenly split between lower and upper concordant diagonal structure in this symmetric correlated example.

---

# Part IV — 50-Variable Pairwise Scalability

## 17. Why Full Orthants Fail at 50 Variables

For $d=50$, the full orthant state count is:

```math
2^{50}
=
1.1259\times 10^{15}.
```

With $n=5000$ observations, the observation-to-orthant ratio is:

```math
\frac{5000}{2^{50}}
=
4.440892\times 10^{-12}.
```

This is not estimable as a full orthant partition. Almost every possible orthant is unobserved, and any occupied orthant would be too sparse for stable covariance estimation.

Pairwise `PM.matrix` avoids this problem by estimating four $50\times 50$ matrices.

---

## 18. 50-Variable Common/Sector Simulation

The 50-variable simulation used:

- one common market factor;
- two sector factors;
- idiosyncratic noise.

Variables 1–25 loaded on sector A. Variables 26–50 loaded on sector B.

```r
set.seed(2026)

n <- 5000
d <- 50

market  <- rnorm(n)
sector_A <- rnorm(n)
sector_B <- rnorm(n)
E <- matrix(rnorm(n * d), nrow = n, ncol = d)

Z50 <- matrix(0, nrow = n, ncol = d)

for (j in seq_len(d)) {
  if (j <= 25) {
    Z50[, j] <- 0.60 * market + 0.90 * sector_A + 0.40 * E[, j]
  } else {
    Z50[, j] <- 0.60 * market + 0.90 * sector_B + 0.40 * E[, j]
  }
}
```

The pairwise PM.matrix experiment ran quickly:

```text
Runtime:
elapsed = 0.12 seconds
```

The first five classical PCA eigenvalues were:

| Component | Eigenvalue |
|---:|---:|
| PC1 | 38.576663 |
| PC2 | 20.499903 |
| PC3 | 0.194593 |
| PC4 | 0.188941 |
| PC5 | 0.186759 |

The pairwise NNS recovered eigenvalues were identical at displayed precision.

Recovery diagnostics:

```text
Pairwise covariance recovery error:
3.774758e-15

Pairwise eigenvalue recovery error:
1.776357e-14

Pairwise eigenvector alignments, first 5 PCs:
1 1 1 1 1
```

Thus pairwise partial-moment matrices recover covariance and PCA exactly for this 50-variable example.

---

## 19. 50-Variable Directional Attribution

The directional eigenvalue attribution for the first five PCs was:

| PC | Eigenvalue | CLPM | CUPM | DLPM | DUPM | Recovered |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 38.576663 | 22.634337 | 22.669687 | 3.363681 | 3.363681 | 38.576663 |
| 2 | 20.499903 | 7.000469 | 7.585803 | -2.956815 | -2.956815 | 20.499903 |
| 3 | 0.194593 | 0.084791 | 0.082119 | -0.013842 | -0.013842 | 0.194593 |
| 4 | 0.188941 | 0.082029 | 0.082308 | -0.012302 | -0.012302 | 0.188941 |
| 5 | 0.186759 | 0.082447 | 0.080548 | -0.011882 | -0.011882 | 0.186759 |

The signed percentage attribution was:

| PC | CLPM % | CUPM % | DLPM signed % | DUPM signed % |
|---:|---:|---:|---:|---:|
| 1 | 58.67 | 58.77 | -8.72 | -8.72 |
| 2 | 34.15 | 37.00 | 14.42 | 14.42 |
| 3 | 43.57 | 42.20 | 7.11 | 7.11 |
| 4 | 43.42 | 43.56 | 6.51 | 6.51 |
| 5 | 44.15 | 43.13 | 6.36 | 6.36 |

Interpretation:

- PC1 is dominated by lower and upper concordance.
- Divergent pairwise structure subtracts from PC1.
- PC2 has substantial positive signed divergent contribution, consistent with a sector/spread direction separating the two blocks.

The trace-level decomposition was:

| Component | Trace Contribution | Trace Percent |
|---|---:|---:|
| CLPM | 33.08720 | 49.57 |
| CUPM | 33.66671 | 50.43 |
| DLPM signed | 0.00000 | 0.00 |
| DUPM signed | 0.00000 | 0.00 |
| Total | 66.75391 | 100.00 |

The trace result says that total variance is evenly split between lower and upper concordant diagonal structure, while divergence affects off-diagonal covariance and eigenvector directions rather than total variance.

---

# Part V — 50-Variable PC1-Not-Market Stress Test

## 20. Motivation

A common PCA interpretation mistake is to assume:

> PC1 is the market factor.

But PC1 is only the largest variance direction. It is not necessarily the market direction.

To test the pairwise directional framework, construct a 50-variable simulation where the spread factor is deliberately stronger than the market factor.

The market vector is:

```math
m
=
\frac{1}{\sqrt{50}}(1,1,\ldots,1).
```

The spread vector is:

```math
s
=
\frac{1}{\sqrt{50}}(\underbrace{1,\ldots,1}_{25},
\underbrace{-1,\ldots,-1}_{25}).
```

The two vectors are orthogonal:

```math
m^\top s = 0.
```

The data are generated as:

```math
Z
=
\sigma_m M m^\top
+
\sigma_s S s^\top
+
\sigma_e E,
```

with:

```math
\sigma_m = 0.70,
\qquad
\sigma_s = 1.80,
\qquad
\sigma_e = 0.35.
```

Because the spread factor is stronger than the market factor, PC1 should be the spread factor, not the market factor.

---

## 21. PCA Result: PC1 Is Not Market

The first five PCA eigenvalues were:

| Component | Eigenvalue |
|---:|---:|
| PC1 | 3.492035 |
| PC2 | 0.618253 |
| PC3 | 0.144454 |
| PC4 | 0.141775 |
| PC5 | 0.140824 |

Alignment diagnostics were:

| Quantity | Value |
|---|---:|
| $\left\lvert\langle PC1,\mathrm{market}\rangle\right\rvert$ | 0.00000742 |
| $\left\lvert\langle PC1,\mathrm{spread}\rangle\right\rvert$ | 0.99982464 |
| $\left\lvert\langle PC2,\mathrm{market}\rangle\right\rvert$ | 0.99827428 |
| $\left\lvert\langle PC2,\mathrm{spread}\rangle\right\rvert$ | 0.00016201 |

Thus PC1 is essentially the spread factor, while PC2 is essentially the market factor.

This confirms:

> PC1 is not always the market factor.

---

## 22. Pairwise Recovery in the Stress Test

The pairwise PM.matrix decomposition again recovered covariance and PCA to numerical precision:

```text
Pairwise covariance recovery error:
4.440892e-16

Pairwise eigenvalue recovery error:
2.664535e-15

Pairwise eigenvector alignments, first 5 PCs:
1 1 1 1 1
```

So even when PC1 is not market, the pairwise directional matrices recover the full PCA eigensystem exactly.

---

## 23. Directional Attribution in the Stress Test

The directional eigenvalue attribution was:

| PC | Eigenvalue | CLPM | CUPM | DLPM | DUPM | Recovered |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 3.492035 | 0.924240 | 0.908033 | -0.829881 | -0.829881 | 3.492035 |
| 2 | 0.618253 | 1.851522 | 1.846249 | 1.539759 | 1.539759 | 0.618253 |
| 3 | 0.144454 | 0.052954 | 0.050886 | -0.020307 | -0.020307 | 0.144454 |
| 4 | 0.141775 | 0.050863 | 0.050608 | -0.020152 | -0.020152 | 0.141775 |
| 5 | 0.140824 | 0.051385 | 0.052016 | -0.018711 | -0.018711 | 0.140824 |

The signed percentage attribution was:

| PC | CLPM % | CUPM % | DLPM signed % | DUPM signed % |
|---:|---:|---:|---:|---:|
| 1 | 26.47 | 26.00 | 23.76 | 23.76 |
| 2 | 299.48 | 298.62 | -249.05 | -249.05 |
| 3 | 36.66 | 35.23 | 14.06 | 14.06 |
| 4 | 35.88 | 35.70 | 14.21 | 14.21 |
| 5 | 36.49 | 36.94 | 13.29 | 13.29 |

For PC1:

```math
3.492035
=
0.924240
+
0.908033
-
(-0.829881)
-
(-0.829881).
```

So PC1 is not only concordant. It is substantially generated by divergent pairwise structure.

Interpretation:

> When PC1 is a market/common factor, concordant matrices dominate and divergence subtracts.  
> When PC1 is a spread factor, divergent pairwise structure becomes part of the positive explanatory mechanism.

This is the high-dimensional pairwise analogue of the full orthant result: PCA finds the dominant axis, while NNS identifies the directional source of that axis.

---

# Part VI — Runtime and Scaling

## 24. Runtime Scaling Grid

The final experiment compared pairwise PM.matrix recovery across dimensions:

| d | Full Orthant Count | Runtime Seconds | Covariance Error | Eigenvalue Error | PC1 Alignment |
|---:|---:|---:|---:|---:|---:|
| 5 | 32 | 0.00 | 8.8817842e-16 | 1.7763568e-15 | 1 |
| 10 | 1,024 | 0.00 | 1.1102230e-15 | 3.5527137e-15 | 1 |
| 25 | 33,554,432 | 0.02 | 1.4432899e-15 | 3.5527137e-15 | 1 |
| 50 | 1.1258999e15 | 0.07 | 1.5543122e-15 | 2.1316282e-14 | 1 |

This table shows the scaling contrast directly:

```math
\mathrm{full\ orthants}
=
O(2^d),
```

while:

```math
\mathrm{pairwise\ PM\ matrices}
=
O(d^2).
```

The pairwise method recovered covariance and PC1 exactly at displayed precision across all tested dimensions.

---

# Part VII — Final Hierarchy

The updated hierarchy is:

```math
\mathrm{full\ orthant\ decomposition}
\Rightarrow
\mathrm{complete\ spectral\ genealogy,\ exact\ but\ exponential}.
```

```math
\mathrm{pairwise\ PM\ matrices}
\Rightarrow
\mathrm{exact\ covariance/PCA\ recovery\ and\ directional\ attribution,\ scalable}.
```

```math
\mathrm{DPM}_{nD}
\Rightarrow
\mathrm{three\ state\ global\ directional\ diagnostic,\ scalable\ but\ coarse}.
```

The three layers answer different questions.

| Layer | State count / storage | What it gives | What it cannot give |
|---|---:|---|---|
| Full orthants | $2^d$ states | Complete orthant-level PCA genealogy | Practical high-dimensional estimation |
| Pairwise PM matrices | Four $d\times d$ matrices | Exact covariance/PCA recovery; directional pairwise attribution | Full simultaneous $d$-orthant labels |
| `DPM_nD` | Three values | All-lower / all-upper / mixed global directional summary | PCA recovery or full spectral genealogy |

The strongest practical conclusion is:

> Use full orthants when $d$ is small and complete regime genealogy is desired.  
> Use pairwise `PM.matrix` for scalable high-dimensional PCA recovery and directional attribution.  
> Use `DPM_nD` as a compact global diagnostic of all-lower, all-upper, and mixed directional mass.

---

# Appendix: Full R Code

```r
# =============================================================================
# Experiments for NNS Directional Spectral Decomposition
# Full Orthants, DPM_nD Aggregation, and Pairwise PM.matrix Scalability
# =============================================================================

library(NNS)

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------

cat_section <- function(title) {
  cat("\n\n============================================================\n")
  cat(title, "\n")
  cat("============================================================\n")
}

pop_cov <- function(Z) {
  Z <- as.matrix(Z)
  Zc <- sweep(Z, 2, colMeans(Z))
  crossprod(Zc) / nrow(Z)
}

align_eigenvectors <- function(V_ref, V_test, k = ncol(V_ref)) {
  sapply(seq_len(k), function(j) {
    abs(sum(V_ref[, j] * V_test[, j]))
  })
}

decode_orthant <- function(label, d, names_vec) {
  bits <- as.integer(intToBits(as.integer(label)))[1:d]
  signs <- ifelse(bits == 1, "+", "-")
  paste(paste0(names_vec, signs), collapse = " ")
}

full_orthant_decomposition <- function(Z) {
  Z <- as.matrix(Z)
  n <- nrow(Z)
  d <- ncol(Z)
  mu <- colMeans(Z)

  above_mean <- sweep(Z, 2, mu, ">")
  orthant_label <- apply(above_mean, 1, function(row) {
    sum(row * 2^(0:(d - 1)))
  })

  Sigma <- pop_cov(Z)
  Sigma_Q <- matrix(0, d, d)
  Sigma_W <- matrix(0, d, d)

  orthant_info <- data.frame(
    orthant = integer(),
    pattern = character(),
    n = integer(),
    p = numeric(),
    u_norm = numeric(),
    stringsAsFactors = FALSE
  )

  for (lab in sort(unique(orthant_label))) {
    mask <- orthant_label == lab
    n_r <- sum(mask)
    p_r <- n_r / n

    Zr <- Z[mask, , drop = FALSE]
    m_r <- colMeans(Zr)
    u_r <- m_r - mu

    Sigma_r <- pop_cov(Zr)

    Sigma_Q <- Sigma_Q + p_r * tcrossprod(u_r)
    Sigma_W <- Sigma_W + p_r * Sigma_r

    orthant_info <- rbind(
      orthant_info,
      data.frame(
        orthant = lab,
        pattern = decode_orthant(lab, d, colnames(Z)),
        n = n_r,
        p = p_r,
        u_norm = sqrt(sum(u_r^2)),
        stringsAsFactors = FALSE
      )
    )
  }

  Sigma_rec <- Sigma_Q + Sigma_W

  pca <- eigen(Sigma, symmetric = TRUE)
  pca_rec <- eigen(Sigma_rec, symmetric = TRUE)

  for (j in seq_len(d)) {
    if (sum(pca$vectors[, j] * pca_rec$vectors[, j]) < 0) {
      pca_rec$vectors[, j] <- -pca_rec$vectors[, j]
    }
  }

  attrib <- data.frame(
    PC = seq_len(d),
    eigenvalue = pca$values,
    between = NA_real_,
    within = NA_real_,
    total = NA_real_,
    between_pct = NA_real_,
    error = NA_real_
  )

  for (j in seq_len(d)) {
    v <- pca$vectors[, j, drop = FALSE]
    between_j <- drop(t(v) %*% Sigma_Q %*% v)
    within_j  <- drop(t(v) %*% Sigma_W %*% v)

    attrib$between[j] <- between_j
    attrib$within[j] <- within_j
    attrib$total[j] <- between_j + within_j
    attrib$between_pct[j] <- 100 * between_j / attrib$total[j]
    attrib$error[j] <- abs(attrib$eigenvalue[j] - attrib$total[j])
  }

  v1 <- pca$vectors[, 1]

  orthant_pc1 <- data.frame(
    orthant = integer(),
    pattern = character(),
    p = numeric(),
    contribution = numeric(),
    pct_lambda1 = numeric(),
    stringsAsFactors = FALSE
  )

  for (lab in sort(unique(orthant_label))) {
    mask <- orthant_label == lab
    p_r <- sum(mask) / n
    u_r <- colMeans(Z[mask, , drop = FALSE]) - mu
    contrib <- p_r * sum(v1 * u_r)^2

    orthant_pc1 <- rbind(
      orthant_pc1,
      data.frame(
        orthant = lab,
        pattern = decode_orthant(lab, d, colnames(Z)),
        p = p_r,
        contribution = contrib,
        pct_lambda1 = 100 * contrib / pca$values[1],
        stringsAsFactors = FALSE
      )
    )
  }

  orthant_pc1 <- orthant_pc1[order(-orthant_pc1$contribution), ]

  list(
    Sigma = Sigma,
    Sigma_Q = Sigma_Q,
    Sigma_W = Sigma_W,
    Sigma_rec = Sigma_rec,
    pca = pca,
    pca_rec = pca_rec,
    orthant_label = orthant_label,
    orthant_info = orthant_info,
    attrib = attrib,
    orthant_pc1 = orthant_pc1,
    recovery_error = max(abs(Sigma - Sigma_rec)),
    eigenvalue_error = max(abs(pca$values - pca_rec$values)),
    eigenvector_alignment = align_eigenvectors(pca$vectors, pca_rec$vectors, d),
    D_spectral = sum(diag(Sigma_Q)) / sum(diag(Sigma))
  )
}

pairwise_pm_pca <- function(Z, pcs = 5) {
  Z <- as.matrix(Z)
  d <- ncol(Z)

  Sigma_classic <- cov(Z)
  pca_classic <- eigen(Sigma_classic, symmetric = TRUE)

  pm <- NNS::PM.matrix(
    LPM_degree = 1,
    UPM_degree = 1,
    target     = "mean",
    variable   = Z,
    pop_adj    = TRUE,
    norm       = FALSE
  )

  CLPM <- pm$clpm
  CUPM <- pm$cupm
  DLPM <- pm$dlpm
  DUPM <- pm$dupm

  Sigma_nns_raw <- CLPM + CUPM - DLPM - DUPM
  Sigma_nns <- (Sigma_nns_raw + t(Sigma_nns_raw)) / 2

  pca_nns <- eigen(Sigma_nns, symmetric = TRUE)

  for (j in seq_len(d)) {
    if (sum(pca_classic$vectors[, j] * pca_nns$vectors[, j]) < 0) {
      pca_nns$vectors[, j] <- -pca_nns$vectors[, j]
    }
  }

  pcs <- min(pcs, d)

  attrib <- data.frame(
    PC = seq_len(pcs),
    eigenvalue = pca_classic$values[1:pcs],
    CLPM = NA_real_,
    CUPM = NA_real_,
    DLPM = NA_real_,
    DUPM = NA_real_,
    recovered = NA_real_,
    error = NA_real_
  )

  for (j in seq_len(pcs)) {
    v <- pca_classic$vectors[, j, drop = FALSE]

    clpm_j <- drop(t(v) %*% CLPM %*% v)
    cupm_j <- drop(t(v) %*% CUPM %*% v)
    dlpm_j <- drop(t(v) %*% DLPM %*% v)
    dupm_j <- drop(t(v) %*% DUPM %*% v)

    recovered_j <- clpm_j + cupm_j - dlpm_j - dupm_j

    attrib$CLPM[j] <- clpm_j
    attrib$CUPM[j] <- cupm_j
    attrib$DLPM[j] <- dlpm_j
    attrib$DUPM[j] <- dupm_j
    attrib$recovered[j] <- recovered_j
    attrib$error[j] <- abs(attrib$eigenvalue[j] - recovered_j)
  }

  attrib_pct <- data.frame(
    PC = attrib$PC,
    CLPM_pct = 100 * attrib$CLPM / attrib$eigenvalue,
    CUPM_pct = 100 * attrib$CUPM / attrib$eigenvalue,
    DLPM_signed_pct = -100 * attrib$DLPM / attrib$eigenvalue,
    DUPM_signed_pct = -100 * attrib$DUPM / attrib$eigenvalue
  )

  trace_fun <- function(M) sum(diag(M))

  trace_decomp <- c(
    CLPM = trace_fun(CLPM),
    CUPM = trace_fun(CUPM),
    DLPM_signed = -trace_fun(DLPM),
    DUPM_signed = -trace_fun(DUPM),
    Total = trace_fun(Sigma_nns)
  )

  list(
    pm = pm,
    CLPM = CLPM,
    CUPM = CUPM,
    DLPM = DLPM,
    DUPM = DUPM,
    Sigma_classic = Sigma_classic,
    Sigma_nns_raw = Sigma_nns_raw,
    Sigma_nns = Sigma_nns,
    pca_classic = pca_classic,
    pca_nns = pca_nns,
    covariance_error = max(abs(Sigma_classic - Sigma_nns_raw)),
    pm_cov_error = max(abs(pm$cov.matrix - Sigma_nns_raw)),
    asymmetry = max(abs(Sigma_nns_raw - t(Sigma_nns_raw))),
    eigenvalue_error = max(abs(pca_classic$values - pca_nns$values)),
    eigenvector_alignment = align_eigenvectors(
      pca_classic$vectors,
      pca_nns$vectors,
      pcs
    ),
    attrib = attrib,
    attrib_pct = attrib_pct,
    trace_decomp = trace_decomp,
    trace_pct = 100 * trace_decomp / trace_decomp["Total"]
  )
}

dpm_nd_summary <- function(Z) {
  Z <- as.matrix(Z)
  n <- nrow(Z)
  d <- ncol(Z)
  target <- colMeans(Z)

  above_target <- sweep(Z, 2, target, ">")
  orthant_label <- apply(above_target, 1, function(row) {
    sum(row * 2^(0:(d - 1)))
  })

  orthant_tab <- table(orthant_label)

  p_all_lower <- as.numeric(orthant_tab[as.character(0)]) / n
  p_all_upper <- as.numeric(orthant_tab[as.character(2^d - 1)]) / n

  if (is.na(p_all_lower)) p_all_lower <- 0
  if (is.na(p_all_upper)) p_all_upper <- 0

  p_mixed <- 1 - p_all_lower - p_all_upper

  orthant_aggregated <- c(
    CLPM_nD = p_all_lower,
    CUPM_nD = p_all_upper,
    DPM_nD  = p_mixed
  )

  clpm0 <- NNS::Co.LPM_nD(Z, target, degree = 0, norm = FALSE)
  cupm0 <- NNS::Co.UPM_nD(Z, target, degree = 0, norm = FALSE)
  dpm0  <- NNS::DPM_nD( Z, target, degree = 0, norm = FALSE)

  nns_degree0 <- c(
    CLPM_nD = clpm0,
    CUPM_nD = cupm0,
    DPM_nD  = dpm0
  )

  clpm1 <- NNS::Co.LPM_nD(Z, target, degree = 1, norm = FALSE)
  cupm1 <- NNS::Co.UPM_nD(Z, target, degree = 1, norm = FALSE)
  dpm1  <- NNS::DPM_nD( Z, target, degree = 1, norm = FALSE)

  nns_degree1_raw <- c(
    CLPM_nD = clpm1,
    CUPM_nD = cupm1,
    DPM_nD  = dpm1
  )

  nns_degree1_norm_manual <- nns_degree1_raw / sum(nns_degree1_raw)

  nns_degree1_norm_true <- c(
    CLPM_nD = NNS::Co.LPM_nD(Z, target, degree = 1, norm = TRUE),
    CUPM_nD = NNS::Co.UPM_nD(Z, target, degree = 1, norm = TRUE),
    DPM_nD  = NNS::DPM_nD( Z, target, degree = 1, norm = TRUE)
  )

  list(
    orthant_aggregated = orthant_aggregated,
    nns_degree0 = nns_degree0,
    degree0_difference = nns_degree0 - orthant_aggregated,
    nns_degree1_raw = nns_degree1_raw,
    nns_degree1_norm_manual = nns_degree1_norm_manual,
    nns_degree1_norm_true = nns_degree1_norm_true
  )
}

# =============================================================================
# Experiment 1: 5-variable full orthant decomposition
# =============================================================================

cat_section("EXPERIMENT 1: 5-VARIABLE FULL ORTHANT DECOMPOSITION")

set.seed(2024)

n <- 10000
d <- 5

R <- matrix(0.5, d, d) + diag(0.5, d)
L <- chol(R)

Z5 <- matrix(rnorm(n * d), n, d) %*% L
colnames(Z5) <- paste0("X", seq_len(d))

res5_orth <- full_orthant_decomposition(Z5)

cat("\nTarget vector:\n")
print(round(colMeans(Z5), 6))

cat("\nOccupied orthants:\n")
print(length(unique(res5_orth$orthant_label)))
cat("Out of:\n")
print(2^d)

cat("\nClassical population PCA eigenvalues:\n")
print(round(res5_orth$pca$values, 6))

cat("\nFull orthant covariance recovery error:\n")
print(res5_orth$recovery_error)

cat("\nFull orthant eigenvalue recovery error:\n")
print(res5_orth$eigenvalue_error)

cat("\nFull orthant eigenvector alignments:\n")
print(round(res5_orth$eigenvector_alignment, 12))

cat("\nEigenvalue attribution, full orthants:\n")
print(round(res5_orth$attrib, 6))

cat("\nD_spectral = trace(Sigma_Q) / trace(Sigma):\n")
print(round(res5_orth$D_spectral, 6))

cat("\nTop 10 orthant contributions to PC1:\n")
print(head(res5_orth$orthant_pc1, 10), digits = 6)

cat("\nCheck PC1 orthant contribution sum vs between attribution:\n")
print(c(
  sum_orthant_PC1 = sum(res5_orth$orthant_pc1$contribution),
  direct_between_PC1 = res5_orth$attrib$between[1]
))

# =============================================================================
# Experiment 2: 5-variable DPM_nD aggregation
# =============================================================================

cat_section("EXPERIMENT 2: 5-VARIABLE DPM_nD AGGREGATION")

res5_dpm <- dpm_nd_summary(Z5)

cat("\nThree-state aggregation from full orthant labels:\n")
print(round(res5_dpm$orthant_aggregated, 6))

cat("\nNNS degree-zero values:\n")
print(round(res5_dpm$nns_degree0, 6))

cat("\nDegree-zero difference, NNS minus full orthant aggregation:\n")
print(round(res5_dpm$degree0_difference, 12))

cat("\nNNS degree-one raw directional masses:\n")
print(round(res5_dpm$nns_degree1_raw, 6))

cat("\nNNS degree-one normalized shares, manual:\n")
print(round(res5_dpm$nns_degree1_norm_manual, 6))

cat("\nNNS degree-one normalized shares, norm = TRUE:\n")
print(round(res5_dpm$nns_degree1_norm_true, 8))

# =============================================================================
# Experiment 3: 5-variable pairwise PM.matrix recovery
# =============================================================================

cat_section("EXPERIMENT 3: 5-VARIABLE PAIRWISE PM.matrix PCA RECOVERY")

res5_pair <- pairwise_pm_pca(Z5, pcs = 5)

cat("\nPairwise PM.matrix covariance recovery error:\n")
print(res5_pair$covariance_error)

cat("\nPM.matrix $cov.matrix vs manual reassembly error:\n")
print(res5_pair$pm_cov_error)

cat("\nDirectional aggregate asymmetry before symmetrization:\n")
print(res5_pair$asymmetry)

cat("\nClassical sample PCA eigenvalues:\n")
print(round(res5_pair$pca_classic$values, 6))

cat("\nPairwise NNS recovered PCA eigenvalues:\n")
print(round(res5_pair$pca_nns$values, 6))

cat("\nPairwise eigenvalue recovery error:\n")
print(res5_pair$eigenvalue_error)

cat("\nPairwise eigenvector alignments:\n")
print(round(res5_pair$eigenvector_alignment, 12))

cat("\nDirectional eigenvalue attribution from pairwise matrices:\n")
print(round(res5_pair$attrib, 6))

cat("\nSigned percentage attribution from pairwise matrices:\n")
print(round(res5_pair$attrib_pct, 2))

cat("\nTrace-level pairwise directional decomposition:\n")
print(round(res5_pair$trace_decomp, 6))

cat("\nTrace-level pairwise directional percentages:\n")
print(round(res5_pair$trace_pct, 2))

# =============================================================================
# Experiment 4: 50-variable pairwise recovery
# =============================================================================

cat_section("EXPERIMENT 4: 50-VARIABLE PAIRWISE PM.matrix SCALABILITY")

set.seed(2026)

n <- 5000
d <- 50

asset_names <- paste0("X", seq_len(d))

market <- rnorm(n)
sector_A <- rnorm(n)
sector_B <- rnorm(n)
E <- matrix(rnorm(n * d), nrow = n, ncol = d)

Z50 <- matrix(0, nrow = n, ncol = d)

for (j in seq_len(d)) {
  if (j <= 25) {
    Z50[, j] <- 0.60 * market + 0.90 * sector_A + 0.40 * E[, j]
  } else {
    Z50[, j] <- 0.60 * market + 0.90 * sector_B + 0.40 * E[, j]
  }
}

colnames(Z50) <- asset_names

cat("\nFull orthant state count for d = 50:\n")
print(2^50)

cat("\nNumber of observations:\n")
print(n)

cat("\nObservation-to-orthant ratio n / 2^50:\n")
print(n / 2^50)

time_50 <- system.time({
  res50_pair <- pairwise_pm_pca(Z50, pcs = 5)
})

cat("\nRuntime for 50-variable pairwise PM.matrix experiment:\n")
print(time_50)

cat("\nClassical PCA eigenvalues, first 5:\n")
print(round(res50_pair$pca_classic$values[1:5], 6))

cat("\nPairwise NNS recovered PCA eigenvalues, first 5:\n")
print(round(res50_pair$pca_nns$values[1:5], 6))

cat("\nPairwise covariance recovery error, d = 50:\n")
print(res50_pair$covariance_error)

cat("\nPairwise eigenvalue recovery error, d = 50:\n")
print(res50_pair$eigenvalue_error)

cat("\nPairwise eigenvector alignments, first 5 PCs, d = 50:\n")
print(round(res50_pair$eigenvector_alignment, 12))

cat("\nDirectional eigenvalue attribution, first 5 PCs, d = 50:\n")
print(round(res50_pair$attrib, 6))

cat("\nSigned percentage attribution, first 5 PCs, d = 50:\n")
print(round(res50_pair$attrib_pct, 2))

cat("\nTrace-level pairwise directional decomposition, d = 50:\n")
print(round(res50_pair$trace_decomp, 6))

cat("\nTrace-level pairwise directional percentages, d = 50:\n")
print(round(res50_pair$trace_pct, 2))

# =============================================================================
# Experiment 5: 50-variable PC1 not market / spread-factor stress test
# =============================================================================

cat_section("EXPERIMENT 5: 50-VARIABLE PC1 NOT MARKET STRESS TEST")

set.seed(2027)

n <- 6000
d <- 50

market_vec <- rep(1, d)
market_vec <- market_vec / sqrt(sum(market_vec^2))

spread_raw <- c(rep(1, 25), rep(-1, 25))
spread_vec <- spread_raw / sqrt(sum(spread_raw^2))

cat("\nMarket-spread inner product:\n")
print(sum(market_vec * spread_vec))

sd_market <- 0.70
sd_spread <- 1.80
sd_noise <- 0.35

M <- rnorm(n)
S <- rnorm(n)
E <- matrix(rnorm(n * d), nrow = n, ncol = d)

Z50_spread <- sd_market * M %*% t(market_vec) +
  sd_spread * S %*% t(spread_vec) +
  sd_noise * E

colnames(Z50_spread) <- paste0("X", seq_len(d))

res50_spread <- pairwise_pm_pca(Z50_spread, pcs = 5)

pca_spread <- res50_spread$pca_classic
v1 <- pca_spread$vectors[, 1]
v2 <- pca_spread$vectors[, 2]

if (sum(v1 * spread_vec) < 0) v1 <- -v1
if (sum(v2 * market_vec) < 0) v2 <- -v2

cat("\nClassical PCA eigenvalues, first 5:\n")
print(round(pca_spread$values[1:5], 6))

cat("\nAlignment diagnostics:\n")
alignment_diag <- c(
  abs_PC1_market = abs(sum(v1 * market_vec)),
  abs_PC1_spread = abs(sum(v1 * spread_vec)),
  abs_PC2_market = abs(sum(v2 * market_vec)),
  abs_PC2_spread = abs(sum(v2 * spread_vec))
)
print(round(alignment_diag, 8))

cat("\nPairwise covariance recovery error, 50-variable spread test:\n")
print(res50_spread$covariance_error)

cat("\nPairwise eigenvalue recovery error, 50-variable spread test:\n")
print(res50_spread$eigenvalue_error)

cat("\nPairwise eigenvector alignments, first 5 PCs, 50-variable spread test:\n")
print(round(res50_spread$eigenvector_alignment, 12))

cat("\nDirectional eigenvalue attribution, first 5 PCs, 50-variable spread test:\n")
print(round(res50_spread$attrib, 6))

cat("\nSigned percentage attribution, first 5 PCs, 50-variable spread test:\n")
print(round(res50_spread$attrib_pct, 2))

# =============================================================================
# Experiment 6: Runtime and scaling grid
# =============================================================================

cat_section("EXPERIMENT 6: PAIRWISE PM.matrix RUNTIME SCALING GRID")

set.seed(2028)

dims <- c(5, 10, 25, 50)
n <- 3000

scaling_results <- data.frame(
  d = integer(),
  impossible_orthants = numeric(),
  runtime_elapsed_sec = numeric(),
  covariance_error = numeric(),
  eigenvalue_error = numeric(),
  pc1_alignment = numeric()
)

for (dd in dims) {
  common <- rnorm(n)
  E <- matrix(rnorm(n * dd), nrow = n, ncol = dd)

  Z <- matrix(0, nrow = n, ncol = dd)
  for (j in seq_len(dd)) {
    Z[, j] <- 0.75 * common + 0.50 * E[, j]
  }
  colnames(Z) <- paste0("X", seq_len(dd))

  timing <- system.time({
    res <- pairwise_pm_pca(Z, pcs = 1)
  })

  scaling_results <- rbind(
    scaling_results,
    data.frame(
      d = dd,
      impossible_orthants = 2^dd,
      runtime_elapsed_sec = timing[["elapsed"]],
      covariance_error = res$covariance_error,
      eigenvalue_error = res$eigenvalue_error,
      pc1_alignment = res$eigenvector_alignment[1]
    )
  )
}

cat("\nScaling results:\n")
print(scaling_results, digits = 8)

cat("\nDONE.\n")
```
