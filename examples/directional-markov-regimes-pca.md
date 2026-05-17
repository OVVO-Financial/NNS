# Directional Markov Regimes and PCA Recovery from NNS Quadrants

This note extends NNS directional spectral decomposition from static quadrants to time-indexed directional regimes.

---

## Executive Summary

Classical PCA starts with a covariance matrix and extracts abstract eigenvectors.

NNS starts with observable directional regions:

- `CUPM`: upper concordant quadrant
- `CLPM`: lower concordant quadrant
- `DLPM`: lower divergent quadrant
- `DUPM`: upper divergent quadrant

Each quadrant has an empirical probability, a conditional mean, and a conditional covariance. These observable pieces reconstruct the covariance matrix exactly.

The static covariance decomposition is:

```math
\Sigma
=
\underbrace{\sum_q p_q u_q u_q^\top}_{\Sigma^B}
+
\underbrace{\sum_q p_q \Sigma^{(q)}}_{\Sigma^W}.
```

The dynamic transition-path decomposition is:

```math
\Sigma_{\mathrm{lead}}
=
\underbrace{
\sum_{q,q'} p_{qq'} u_{q\to q'}u_{q\to q'}^\top
}_{\Sigma^B_{\mathrm{dyn}}}
+
\underbrace{
\sum_{q,q'} p_{qq'} \Sigma^{(q\to q')}
}_{\Sigma^W_{\mathrm{dyn}}}.
```

The main result is:

> PCA identifies the dominant axis. NNS identifies the regimes and transition paths that created it.

---

## 1. Hidden Markov Models Versus Directional Markov Regimes

A classical Hidden Markov Model has:

- latent states `S_t in {1, ..., K}`;
- transition probabilities `P_ij = P(S_{t+1}=j | S_t=i)`;
- emission distributions `f(Z_t | S_t=k)`.

The analyst observes `Z_t`, but not `S_t`. The states must be inferred through filtering, smoothing, or EM-type procedures.

In the NNS directional framework, the state is not latent. It is directly defined by directional geometry.

Let:

```math
Z_t =
\begin{pmatrix}
X_t\\
Y_t
\end{pmatrix},
\qquad
\mu =
\begin{pmatrix}
\bar X\\
\bar Y
\end{pmatrix}.
```

Define the observable quadrant state:

```math
Q_t \in \{\mathrm{CUPM}, \mathrm{CLPM}, \mathrm{DLPM}, \mathrm{DUPM}\}.
```

Using a mean split:

```math
\mathrm{CUPM}: X_t > \bar X,\; Y_t > \bar Y,
```

```math
\mathrm{CLPM}: X_t \leq \bar X,\; Y_t \leq \bar Y,
```

```math
\mathrm{DLPM}: X_t > \bar X,\; Y_t \leq \bar Y,
```

```math
\mathrm{DUPM}: X_t \leq \bar X,\; Y_t > \bar Y.
```

The state probabilities are empirical frequencies:

```math
p_q = P(Q_t=q).
```

The state-conditional means are:

```math
m_q = E[Z_t \mid Q_t=q].
```

The state-conditional covariances are:

```math
\Sigma^{(q)}
=
\mathrm{Cov}(Z_t \mid Q_t=q).
```

The key difference from an HMM is that `Q_t` is observed once the benchmark is specified. No latent-state inference is required.

---

## 2. Static Directional Spectral Decomposition

The covariance matrix decomposes as:

```math
\Sigma
=
\Sigma^B+\Sigma^W,
```

where:

```math
\Sigma^B
=
\sum_q p_q u_q u_q^\top,
```

and:

```math
\Sigma^W
=
\sum_q p_q\Sigma^{(q)}.
```

Here:

```math
u_q=m_q-\mu.
```

The between-quadrant component `Sigma^B` is built from rank-one spectral primitives:

```math
B_q=p_q u_q u_q^\top.
```

If `u_q` is nonzero, then:

```math
B_q u_q
=
p_q u_q u_q^\top u_q
=
p_q\|u_q\|^2u_q.
```

Therefore `u_q` is the nonzero eigenvector of its own quadrant contribution `B_q`, with eigenvalue:

```math
\ell_q=p_q\|u_q\|^2.
```

After normalization,

```math
v_q=\frac{u_q}{\|u_q\|}
```

is the corresponding unit eigenvector.

Thus, centered quadrant conditional means are rank-one spectral primitives.

The full between-quadrant covariance is:

```math
\Sigma^B
=
B_{\mathrm{CUPM}}
+
B_{\mathrm{CLPM}}
+
B_{\mathrm{DLPM}}
+
B_{\mathrm{DUPM}}.
```

Equivalently, define:

```math
C
=
\left(
\sqrt{p_{\mathrm{CUPM}}}u_{\mathrm{CUPM}}\;
\sqrt{p_{\mathrm{CLPM}}}u_{\mathrm{CLPM}}\;
\sqrt{p_{\mathrm{DLPM}}}u_{\mathrm{DLPM}}\;
\sqrt{p_{\mathrm{DUPM}}}u_{\mathrm{DUPM}}
\right).
```

Then:

```math
\Sigma^B=CC^\top.
```

The eigenvectors of `Sigma^B` are the left singular vectors of `C`, built entirely from weighted quadrant conditional mean displacements.

---


## 3. PCA Recovery Directly from Quadrant Conditional Means

This is the central static recovery chain:

```math
\{p_q,m_q,\Sigma^{(q)}\}_q
\longrightarrow
\Sigma^B+\Sigma^W
\longrightarrow
\Sigma
\longrightarrow
(\lambda_i,v_i).
```

The quadrant conditional means enter through:

```math
u_q=m_q-\mu.
```

From those displacements, construct the rank-one matrices:

```math
B_q=p_q u_q u_q^\top.
```

Then:

```math
\Sigma^B
=
\sum_q B_q
=
\sum_q p_q u_q u_q^\top.
```

If only the quadrant conditional means and probabilities are used, the recovered eigensystem is the eigensystem of the between-quadrant covariance `Sigma^B`.

To recover the full classical PCA eigensystem of the original covariance matrix, add the within-quadrant residual covariance terms:

```math
\Sigma
=
\Sigma^B+\Sigma^W
=
\sum_q p_q u_q u_q^\top
+
\sum_q p_q\Sigma^{(q)}.
```

Then diagonalize the recovered covariance matrix:

```math
\Sigma v_i=\lambda_i v_i.
```

This is not an approximation in the sample calculation. It is an exact covariance decomposition up to numerical machine precision.

The R example below verifies this directly by computing:

1. each quadrant conditional mean `m_q`;
2. each centered displacement `u_q`;
3. each rank-one matrix `B_q`;
4. the between-quadrant covariance `Sigma_B`;
5. the within-quadrant covariance `Sigma_W`;
6. the recovered covariance `Sigma_B + Sigma_W`;
7. the eigenvalues and eigenvectors of the recovered matrix.

In the reported run, the max absolute recovery error was approximately:

```math
1.08\times 10^{-19},
```

and the original and recovered leading eigenvectors had alignment:

```math
1.
```

The between-quadrant covariance alone had leading eigenvector alignment:

```math
0.9999992
```

with the full PCA leading eigenvector.

So the result is stronger than merely saying that PCA is related to quadrant means:

> The quadrant conditional means generate rank-one spectral primitives. Their sum recovers the between-quadrant eigensystem. Adding within-quadrant residual covariance recovers the full PCA eigensystem.


## 4. PCA Eigenvalue Attribution

Let `(lambda_i, v_i)` be a unit eigenpair of the full covariance matrix `Sigma`:

```math
\Sigma v_i=\lambda_i v_i,
\qquad
\|v_i\|=1.
```

Since:

```math
\lambda_i = v_i^\top\Sigma v_i,
```

and:

```math
\Sigma
=
\sum_q p_q u_q u_q^\top
+
\sum_q p_q\Sigma^{(q)},
```

we get the exact attribution identity:

```math
\lambda_i
=
\sum_q p_q(v_i^\top u_q)^2
+
\sum_q p_q v_i^\top\Sigma^{(q)}v_i.
```

This decomposes each classical PCA eigenvalue into:

1. between-quadrant conditional-mean displacement;
2. within-quadrant residual covariance.

> PCA diagonalizes covariance. Directional decomposition explains where that covariance came from.

---

## 5. Observable Directional Markov Transitions

Because `Q_t` is observed, transition probabilities can be estimated directly:

```math
\widehat P_{qq'}
=
P(Q_{t+1}=q'\mid Q_t=q)
=
\frac{
\#\{t:Q_t=q,\;Q_{t+1}=q'\}
}{
\#\{t:Q_t=q\}
}.
```

This gives a four-state observable Markov chain over:

```math
\{\mathrm{CUPM},\mathrm{CLPM},\mathrm{DLPM},\mathrm{DUPM}\}.
```

Examples:

```math
P(\mathrm{CLPM}_{t+1}\mid \mathrm{CLPM}_t)
```

is crash persistence.

```math
P(\mathrm{CUPM}_{t+1}\mid \mathrm{CLPM}_t)
```

is crash-to-rally reversal.

```math
P(\mathrm{CLPM}_{t+1}\mid \mathrm{CUPM}_t)
```

is rally-to-crash reversal.

This is not a latent Markov model. It is an observable directional Markov regime model.

---

## 6. One-Step Predictive Mixture

Given current quadrant `Q_t=q`, the one-step predictive mean is:

```math
E[Z_{t+1}\mid Q_t=q]
=
\sum_{q'}P_{qq'}m_{q'}.
```

The corresponding predictive covariance is:

```math
\mathrm{Cov}(Z_{t+1}\mid Q_t=q)
=
\sum_{q'}P_{qq'}
\left[
\Sigma^{(q')}
+
(m_{q'}-\mu_{q\to\cdot})(m_{q'}-\mu_{q\to\cdot})^\top
\right],
```

where:

```math
\mu_{q\to\cdot}
=
\sum_{q'}P_{qq'}m_{q'}.
```

This is the same mixture logic used by HMMs, but with observed directional states instead of inferred hidden states.

---

## 7. Dynamic Transition-Path Spectral Decomposition

The static quadrant decomposition can be extended to transition paths.

Define:

```math
p_{qq'}=P(Q_t=q,Q_{t+1}=q').
```

Let:

```math
m_{q\to q'}
=
E[Z_{t+1}\mid Q_t=q,Q_{t+1}=q'].
```

Let:

```math
\mu_{\mathrm{lead}}=E[Z_{t+1}].
```

Define the transition-path displacement:

```math
u_{q\to q'}
=
m_{q\to q'}-\mu_{\mathrm{lead}}.
```

Define the within-transition covariance:

```math
\Sigma^{(q\to q')}
=
\mathrm{Cov}(Z_{t+1}\mid Q_t=q,Q_{t+1}=q').
```

Then the lead covariance decomposes as:

```math
\Sigma_{\mathrm{lead}}
=
\Sigma^B_{\mathrm{dyn}}
+
\Sigma^W_{\mathrm{dyn}},
```

where:

```math
\Sigma^B_{\mathrm{dyn}}
=
\sum_{q,q'}p_{qq'}u_{q\to q'}u_{q\to q'}^\top,
```

and:

```math
\Sigma^W_{\mathrm{dyn}}
=
\sum_{q,q'}p_{qq'}\Sigma^{(q\to q')}.
```

Each transition path contributes a rank-one dynamic spectral primitive:

```math
B_{q\to q'}
=
p_{qq'}u_{q\to q'}u_{q\to q'}^\top.
```

If `u_{q->q'}` is nonzero, then:

```math
B_{q\to q'}u_{q\to q'}
=
p_{qq'}\|u_{q\to q'}\|^2u_{q\to q'}.
```

Thus transition paths are dynamic rank-one spectral primitives.

---


## 8. Dynamic PCA Recovery from Transition-Path Conditional Means

The dynamic recovery chain is:

```math
\{p_{qq'},m_{q\to q'},\Sigma^{(q\to q')}\}_{q,q'}
\longrightarrow
\Sigma^B_{\mathrm{dyn}}+\Sigma^W_{\mathrm{dyn}}
\longrightarrow
\Sigma_{\mathrm{lead}}
\longrightarrow
(\lambda_i^{\mathrm{lead}},v_i^{\mathrm{lead}}).
```

Transition-path conditional means enter through:

```math
u_{q\to q'}=m_{q\to q'}-\mu_{\mathrm{lead}}.
```

The dynamic between-transition covariance is:

```math
\Sigma^B_{\mathrm{dyn}}
=
\sum_{q,q'}p_{qq'}u_{q\to q'}u_{q\to q'}^\top.
```

Adding the within-transition residual covariance gives the full lead covariance:

```math
\Sigma_{\mathrm{lead}}
=
\Sigma^B_{\mathrm{dyn}}+\Sigma^W_{\mathrm{dyn}}.
```

Then diagonalizing `Sigma_lead` recovers the dynamic PCA eigensystem.

In the reported run, the dynamic recovery error was approximately:

```math
1.08\times 10^{-19},
```

and the dynamic recovered leading eigenvector had alignment:

```math
1
```

with the original lead covariance leading eigenvector.

The dynamic between-transition covariance alone had leading eigenvector alignment:

```math
0.9999993
```

with the full lead PC1.


## 10. Dynamic Eigenvalue Attribution

Let `(lambda_i_lead, v_i_lead)` be an eigenpair of `Sigma_lead`.

Then:

```math
\lambda_i^{\mathrm{lead}}
=
\sum_{q,q'}p_{qq'}(v_i^{\mathrm{lead}\top}u_{q\to q'})^2
+
\sum_{q,q'}p_{qq'}v_i^{\mathrm{lead}\top}\Sigma^{(q\to q')}v_i^{\mathrm{lead}}.
```

This gives dynamic spectral attribution by transition path.

Instead of saying:

> PC1 explains most of the variance,

we can say:

> PC1 is generated by specific observable transition paths, such as CLPM to CUPM, CUPM to CLPM, CLPM to CLPM, and CUPM to CUPM.

This is the observable-regime analogue of dynamic factor attribution.

---

## 10. Numerical Experiment

A bivariate time series of length `n = 5000` was simulated from two hidden volatility regimes:

- state 1: calm, low volatility, low correlation;
- state 2: turbulent, higher volatility, higher correlation.

The true hidden states followed a persistent two-state Markov process with:

```math
p_{\mathrm{stay}}=0.96.
```

The observed NNS states were then defined by quadrant membership relative to the mean of `X` and `Y`.

### Hidden State Frequencies

| State | Frequency |
|---:|---:|
| 1 calm | 0.4082 |
| 2 turbulent | 0.5918 |

### Observable Quadrant Frequencies

| Quadrant | Frequency |
|---|---:|
| CUPM | 0.3412 |
| CLPM | 0.3380 |
| DLPM | 0.1564 |
| DUPM | 0.1644 |

The concordant quadrants dominate:

```math
p_{\mathrm{CUPM}}+p_{\mathrm{CLPM}}
=
0.6792.
```

The divergent quadrants account for:

```math
p_{\mathrm{DLPM}}+p_{\mathrm{DUPM}}
=
0.3208.
```

### Hidden-State Composition by Observable Quadrant

| Hidden state | CUPM | CLPM | DLPM | DUPM |
|---:|---:|---:|---:|---:|
| 1 calm | 0.3388 | 0.3172 | 0.5921 | 0.5645 |
| 2 turbulent | 0.6612 | 0.6828 | 0.4079 | 0.4355 |

The concordant quadrants are mostly turbulent:

```math
P(\mathrm{turbulent}\mid \mathrm{CUPM})=0.6612,
```

```math
P(\mathrm{turbulent}\mid \mathrm{CLPM})=0.6828.
```

The divergent quadrants are more often calm:

```math
P(\mathrm{calm}\mid \mathrm{DLPM})=0.5921,
```

```math
P(\mathrm{calm}\mid \mathrm{DUPM})=0.5645.
```

Thus the observable quadrant states recover meaningful information about the latent volatility regime without fitting an HMM.

---

## 11. Static Results

### Covariance Recovery

The original covariance matrix was:

```math
\Sigma
=
\begin{pmatrix}
0.00078857 & 0.00055621\\
0.00055621 & 0.00079558
\end{pmatrix}.
```

The recovered covariance matrix was:

```math
\Sigma^B+\Sigma^W
=
\begin{pmatrix}
0.00078857 & 0.00055621\\
0.00055621 & 0.00079558
\end{pmatrix}.
```

The max absolute recovery error was:

```math
1.084202\times 10^{-19}.
```

This is numerical zero.

### Static Eigenvalue Recovery

The eigenvalues of the original covariance matrix were:

```math
\lambda_1=0.0013482978,
\qquad
\lambda_2=0.0002358557.
```

The eigenvalues of the recovered matrix were identical.

The leading eigenvector of the original covariance matrix was:

```math
v_1=(0.7048772,\;0.7093294).
```

The leading eigenvector of the recovered matrix was:

```math
v_1=(0.7048772,\;0.7093294).
```

The alignment was:

```math
1.
```

The leading eigenvector of the between-quadrant covariance `Sigma^B` alone was:

```math
v^B_1=(0.7039744,\;0.7102253).
```

Its alignment with the full PC1 was:

```math
0.9999992.
```

Thus, in this simulation, the quadrant conditional-mean geometry essentially determines PC1.

### Static Eigenvalue Attribution

For the first eigenvalue:

```math
\lambda_1=0.0013482978.
```

The decomposition was:

```math
\text{between-quadrant}=0.0008537907,
```

```math
\text{within-quadrant}=0.0004945070.
```

Therefore:

```math
\frac{0.0008537907}{0.0013482978}
\approx 63.3\%.
```

About 63.3 percent of the leading PCA eigenvalue came from between-quadrant conditional-mean displacement.

About 36.7 percent came from within-quadrant residual covariance.

### Per-Quadrant Contributions to `lambda_1`

| Quadrant | Total Contribution | Share of `lambda_1` |
|---|---:|---:|
| CLPM | 0.0006692665 | 49.6% |
| CUPM | 0.0006446941 | 47.8% |
| DUPM | 0.0000182523 | 1.4% |
| DLPM | 0.0000160849 | 1.2% |

The concordant quadrants contributed approximately:

```math
49.6\%+47.8\%=97.4\%
```

of the leading eigenvalue.

Therefore:

> The leading eigenvalue is a concordant co-movement eigenvalue.

---

## 12. Observable Transition Matrix

The estimated transition matrix was:

| Current `Q_t` | CUPM | CLPM | DLPM | DUPM |
|---|---:|---:|---:|---:|
| CUPM | 0.3464 | 0.3470 | 0.1512 | 0.1553 |
| CLPM | 0.3503 | 0.3485 | 0.1444 | 0.1568 |
| DLPM | 0.3171 | 0.3235 | 0.1739 | 0.1854 |
| DUPM | 0.3350 | 0.3118 | 0.1754 | 0.1778 |

Key transition probabilities:

```math
P(\mathrm{CLPM}_{t+1}\mid \mathrm{CLPM}_t)=0.3485.
```

```math
P(\mathrm{CUPM}_{t+1}\mid \mathrm{CLPM}_t)=0.3503.
```

```math
P(\mathrm{CLPM}_{t+1}\mid \mathrm{CUPM}_t)=0.3470.
```

In this simulation, crash persistence is not dominant. CLPM is about equally likely to persist or flip to CUPM.

This reveals an important distinction:

> HMM state persistence does not imply directional quadrant persistence.

The hidden volatility regime is persistent, but the directional quadrant state is not strongly persistent. The NNS transition matrix measures persistence and reversal in observable directional outcomes, not persistence in latent volatility.

---

## 13. One-Step Forecast from CLPM

For current state:

```math
Q_t=\mathrm{CLPM},
```

the one-step predictive mean was:

```math
E[Z_{t+1}\mid Q_t=\mathrm{CLPM}]
=
(-0.0003297,\;-0.0006816).
```

The one-step predictive covariance was:

```math
\mathrm{Cov}(Z_{t+1}\mid Q_t=\mathrm{CLPM})
=
\begin{pmatrix}
0.0008044702 & 0.0005762073\\
0.0005762073 & 0.0008115607
\end{pmatrix}.
```

This is an HMM-style predictive mixture, but no hidden-state filtering was used.

---

## 14. Dynamic Transition-Path Results

### Dynamic Covariance Recovery

The original lead covariance matrix was:

```math
\Sigma_{\mathrm{lead}}
=
\begin{pmatrix}
0.00078851 & 0.00055651\\
0.00055651 & 0.00079558
\end{pmatrix}.
```

The recovered dynamic covariance matrix was identical:

```math
\Sigma^B_{\mathrm{dyn}}+\Sigma^W_{\mathrm{dyn}}
=
\begin{pmatrix}
0.00078851 & 0.00055651\\
0.00055651 & 0.00079558
\end{pmatrix}.
```

The max absolute dynamic recovery error was:

```math
1.084202\times 10^{-19}.
```

Again, this is numerical zero.

### Dynamic Eigensystem Recovery

The eigenvalues of `Sigma_lead` were:

```math
\lambda_1^{\mathrm{lead}}=0.0013485649,
\qquad
\lambda_2^{\mathrm{lead}}=0.0002355231.
```

The recovered dynamic covariance matrix had the same eigenvalues.

The leading eigenvector of `Sigma_lead` was:

```math
v_1^{\mathrm{lead}}=(0.7048570,\;0.7093494).
```

The leading eigenvector of the recovered dynamic covariance was identical.

The dynamic between-transition covariance `Sigma^B_dyn` alone had leading eigenvector:

```math
v^{B,dyn}_1=(0.7040324,\;0.7101679).
```

Its alignment with the full lead PC1 was:

```math
0.9999993.
```

Thus the transition-path conditional-mean geometry essentially determines the dynamic leading eigenvector.

### Dynamic Eigenvalue Attribution

For the leading dynamic eigenvalue:

```math
\lambda_1^{\mathrm{lead}}=0.001348564855.
```

The decomposition was:

```math
\text{between-transition}=0.0008586159410,
```

```math
\text{within-transition}=0.0004899489139.
```

Therefore:

```math
\frac{0.0008586159410}{0.001348564855}
\approx 63.7\%.
```

About 63.7 percent of the dynamic leading eigenvalue came from transition-path conditional-mean displacement.

About 36.3 percent came from within-transition residual covariance.

---

## 15. Transition-Path Contributions to Dynamic `lambda_1`

The largest transition-path contributors to the leading dynamic eigenvalue were:

| Path | Total Contribution | Share of `lambda_1_lead` |
|---|---:|---:|
| CLPM to CUPM | 0.0002541752 | 18.85% |
| CUPM to CLPM | 0.0002489761 | 18.46% |
| CLPM to CLPM | 0.0002451239 | 18.18% |
| CUPM to CUPM | 0.0002324682 | 17.24% |
| DLPM to CLPM | 0.0000955944 | 7.09% |
| DUPM to CUPM | 0.0000838107 | 6.21% |
| DUPM to CLPM | 0.0000797236 | 5.91% |
| DLPM to CUPM | 0.0000743513 | 5.51% |

The top four paths are all concordant-to-concordant paths:

```math
\{\mathrm{CLPM},\mathrm{CUPM}\}
\to
\{\mathrm{CLPM},\mathrm{CUPM}\}.
```

Together they contributed approximately:

```math
18.85\%+18.46\%+18.18\%+17.24\%
=
72.73\%
```

of the leading dynamic eigenvalue.

Therefore:

> Dynamic PC1 is generated primarily by transitions among concordant regimes.

This is stronger than saying PC1 is a correlation factor. It identifies the observable transition paths that generate the dominant eigenstructure.

---

## 16. Interpretation

The results support three claims.

### Claim 1: NNS quadrant decomposition exactly recovers covariance and PCA.

The static covariance recovery error was approximately `1e-19`.

The recovered eigenvalues and eigenvectors matched the original covariance eigensystem exactly to numerical precision.

### Claim 2: Centered quadrant means are rank-one eigenvector primitives.

For every quadrant:

```math
B_q u_q=\ell_q u_q.
```

The numerical errors were zero to machine precision, and the alignment with each rank-one eigenvector was `1`.

### Claim 3: Time-indexed quadrant transitions produce observable Markov regimes with exact spectral attribution.

For every transition path:

```math
B_{q\to q'}u_{q\to q'}
=
\ell_{q\to q'}u_{q\to q'}.
```

The dynamic covariance decomposition recovered `Sigma_lead` exactly, and dynamic eigenvalue attribution decomposed the leading eigenvalue into interpretable transition-path contributions.

The main empirical finding is:

> Static PC1 is driven by CLPM and CUPM.  
> Dynamic PC1 is driven by transitions among CLPM and CUPM.

In other words, PC1 is not mysterious in this simulation. It is a concordant co-movement object.

---

## 17. Portfolio Interpretation

The directional decomposition changes the interpretation of hedging.

Classical PCA suggests hedging against PC1 exposure:

```math
w^\top v_1=0.
```

Directional NNS allows regime-targeted hedging. For example, lower-tail co-movement exposure can be neutralized by:

```math
w^\top u_{\mathrm{CLPM}}=0.
```

This condition means the portfolio is locally neutral to the CLPM conditional-mean displacement, not merely to an abstract variance-maximizing axis.

Dynamic regime hedging can target transition paths:

```math
w^\top u_{\mathrm{CLPM}\to\mathrm{CLPM}}=0,
```

or:

```math
w^\top u_{\mathrm{CLPM}\to\mathrm{CUPM}}=0.
```

These are economically interpretable exposures:

- CLPM to CLPM: crash persistence;
- CLPM to CUPM: crash-to-rally reversal;
- CUPM to CLPM: rally-to-crash reversal;
- CUPM to CUPM: rally persistence.

Thus, directional spectral decomposition provides not only attribution, but also targeted risk control.

---

## 18. Theoretical Significance

Classical PCA begins with the covariance matrix:

```math
\Sigma,
```

and extracts orthogonal eigenvectors:

```math
v_i.
```

The directional framework begins with interpretable regime primitives:

```math
p_q,\quad m_q,\quad u_q,\quad \Sigma^{(q)}.
```

Then it reconstructs:

```math
\Sigma
=
\sum_q p_q u_q u_q^\top
+
\sum_q p_q\Sigma^{(q)}.
```

Therefore the PCA eigensystem is downstream of directional structure:

```math
\{p_q,u_q,\Sigma^{(q)}\}_q
\longrightarrow
\Sigma
\longrightarrow
(\lambda_i,v_i).
```

The map does not generally reverse:

```math
(\lambda_i,v_i)
\not\longrightarrow
\{p_q,u_q,\Sigma^{(q)}\}_q.
```

PCA discards the regime-indexed origin of covariance.

NNS preserves it.

---

## 19. Bottom Line

The central result is:

> PCA diagonalizes covariance. Directional decomposition explains its origin.

The static result says:

> PCA eigenvalues and eigenvectors are recoverable from directional quadrant decomposition.

The dynamic result says:

> PCA eigenvalues and eigenvectors are recoverable from observable directional transition-path decomposition.

The portfolio result says:

> Risk can be attributed and hedged by named directional regimes, not only by abstract orthogonal factors.

The HMM comparison says:

> HMMs infer latent regimes and then interpret them. NNS defines observable directional regimes first and then measures their dynamics.

Together, these results show that directional NNS is not merely an alternative dependence statistic. It is a regime-indexed spectral genealogy for covariance, PCA, and dynamic risk.

---

## Appendix: Full R Code

```r
# =============================================================================
# Directional Markov Regimes
# Observable NNS Quadrant Analogy to Hidden Markov Models
# =============================================================================

set.seed(123)

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------

pop_cov <- function(M) {
  M <- as.matrix(M)
  n <- nrow(M)
  if (n <= 1) return(matrix(0, ncol(M), ncol(M)))
  cov(M) * (n - 1) / n
}

mvrnorm_base <- function(n, mu, Sigma) {
  p <- length(mu)
  Z <- matrix(rnorm(n * p), n, p)
  sweep(Z %*% chol(Sigma), 2, mu, "+")
}

assign_quadrants <- function(Z, center) {
  X <- Z[, 1]
  Y <- Z[, 2]
  cx <- center[1]
  cy <- center[2]

  Q <- rep(NA_character_, nrow(Z))

  Q[X >  cx & Y >  cy] <- "CUPM"
  Q[X <= cx & Y <= cy] <- "CLPM"
  Q[X >  cx & Y <= cy] <- "DLPM"
  Q[X <= cx & Y >  cy] <- "DUPM"

  factor(Q, levels = c("CUPM", "CLPM", "DLPM", "DUPM"))
}

cat_matrix <- function(M, digits = 8) {
  print(round(M, digits))
}

# -----------------------------------------------------------------------------
# 1. Simulate bivariate time series with two hidden volatility regimes
# -----------------------------------------------------------------------------

n <- 5000

# True hidden states:
# 1 = calm
# 2 = turbulent
state <- integer(n)
state[1] <- 1

p_stay <- 0.96

for (t in 2:n) {
  state[t] <- ifelse(runif(1) < p_stay, state[t - 1], 3 - state[t - 1])
}

mu_calm <- c(0.0003, 0.0003)
Sigma_calm <- matrix(
  c(0.0002, 0.00003,
    0.00003, 0.0002),
  nrow = 2,
  byrow = TRUE
)

mu_turb <- c(-0.0008, -0.0012)
Sigma_turb <- matrix(
  c(0.0012, 0.0009,
    0.0009, 0.0012),
  nrow = 2,
  byrow = TRUE
)

Z <- matrix(0, n, 2)

for (t in seq_len(n)) {
  if (state[t] == 1) {
    Z[t, ] <- mvrnorm_base(1, mu_calm, Sigma_calm)
  } else {
    Z[t, ] <- mvrnorm_base(1, mu_turb, Sigma_turb)
  }
}

colnames(Z) <- c("X", "Y")

cat("\n============================================================\n")
cat("SIMULATED DATA\n")
cat("============================================================\n")
cat("Number of observations:", n, "\n")
cat("True hidden state frequencies:\n")
print(prop.table(table(state)))

# -----------------------------------------------------------------------------
# 2. Observable NNS quadrant assignment by mean split
# -----------------------------------------------------------------------------

# For ex post decomposition, use full-sample mean.
# For live forecasting, replace this with a training, expanding, rolling,
# or externally fixed benchmark.
center <- colMeans(Z)

Q <- assign_quadrants(Z, center)

cat("\n============================================================\n")
cat("OBSERVABLE NNS QUADRANT STATES\n")
cat("============================================================\n")
cat("Mean benchmark:\n")
print(center)

cat("\nQuadrant frequencies:\n")
p_quad <- prop.table(table(Q))
print(round(p_quad, 6))

cat("\nCross-tab: true hidden state versus observable quadrant\n")
print(table(state, Q))

cat("\nConditional distribution of hidden states within each quadrant:\n")
print(round(prop.table(table(state, Q), margin = 2), 4))

# -----------------------------------------------------------------------------
# 3. Static NNS covariance decomposition
# -----------------------------------------------------------------------------

mu <- colMeans(Z)
Sigma <- pop_cov(Z)

levels_Q <- levels(Q)

quad_means <- matrix(NA_real_, nrow = 2, ncol = length(levels_Q))
colnames(quad_means) <- levels_Q
rownames(quad_means) <- colnames(Z)

quad_covs <- vector("list", length(levels_Q))
names(quad_covs) <- levels_Q

Sigma_B <- matrix(0, 2, 2)
Sigma_W <- matrix(0, 2, 2)

rank_one_table <- data.frame(
  quadrant = character(),
  p = numeric(),
  mean_X = numeric(),
  mean_Y = numeric(),
  u_X = numeric(),
  u_Y = numeric(),
  ell_rank1 = numeric(),
  max_abs_Bu_minus_ell_u = numeric(),
  alignment_with_Bq_eigenvector = numeric(),
  stringsAsFactors = FALSE
)

cat("\n============================================================\n")
cat("STATIC NNS QUADRANT DECOMPOSITION\n")
cat("============================================================\n")

for (q in levels_Q) {
  idx <- Q == q
  Zq <- Z[idx, , drop = FALSE]

  p_q <- nrow(Zq) / nrow(Z)
  m_q <- colMeans(Zq)
  u_q <- as.numeric(m_q - mu)

  Sigma_q <- pop_cov(Zq)

  B_q <- p_q * tcrossprod(u_q)
  W_q <- p_q * Sigma_q

  Sigma_B <- Sigma_B + B_q
  Sigma_W <- Sigma_W + W_q

  quad_means[, q] <- m_q
  quad_covs[[q]] <- Sigma_q

  ell_q <- p_q * sum(u_q^2)

  lhs <- as.numeric(B_q %*% u_q)
  rhs <- ell_q * u_q
  max_err <- max(abs(lhs - rhs))

  eig_Bq <- eigen(B_q, symmetric = TRUE)
  v_Bq <- eig_Bq$vectors[, 1]

  if (sum(u_q^2) > 0) {
    v_u <- u_q / sqrt(sum(u_q^2))
    align <- abs(sum(v_u * v_Bq))
  } else {
    align <- NA_real_
  }

  rank_one_table <- rbind(
    rank_one_table,
    data.frame(
      quadrant = q,
      p = p_q,
      mean_X = m_q[1],
      mean_Y = m_q[2],
      u_X = u_q[1],
      u_Y = u_q[2],
      ell_rank1 = ell_q,
      max_abs_Bu_minus_ell_u = max_err,
      alignment_with_Bq_eigenvector = align,
      stringsAsFactors = FALSE
    )
  )
}

print(rank_one_table, digits = 8)

Sigma_recovered <- Sigma_B + Sigma_W

cat("\nOriginal population covariance Sigma:\n")
cat_matrix(Sigma)

cat("\nRecovered covariance Sigma_B + Sigma_W:\n")
cat_matrix(Sigma_recovered)

cat("\nMax absolute recovery error:\n")
print(max(abs(Sigma - Sigma_recovered)))

cat("\nBetween-quadrant covariance Sigma_B:\n")
cat_matrix(Sigma_B)

cat("\nWithin-quadrant residual covariance Sigma_W:\n")
cat_matrix(Sigma_W)

# -----------------------------------------------------------------------------
# 4. Static eigensystem recovery and spectral attribution
# -----------------------------------------------------------------------------

eig_Sigma <- eigen(Sigma, symmetric = TRUE)
eig_recovered <- eigen(Sigma_recovered, symmetric = TRUE)
eig_B <- eigen(Sigma_B, symmetric = TRUE)

cat("\n============================================================\n")
cat("STATIC EIGENVALUE AND EIGENVECTOR RECOVERY\n")
cat("============================================================\n")

cat("\nEigenvalues of original Sigma:\n")
print(eig_Sigma$values)

cat("\nEigenvalues of recovered Sigma_B + Sigma_W:\n")
print(eig_recovered$values)

cat("\nLeading eigenvector of original Sigma:\n")
print(eig_Sigma$vectors[, 1])

cat("\nLeading eigenvector of recovered Sigma_B + Sigma_W:\n")
print(eig_recovered$vectors[, 1])

cat("\nLeading eigenvector of Sigma_B only:\n")
print(eig_B$vectors[, 1])

cat("\nAbsolute alignment between original PC1 and recovered PC1:\n")
print(abs(sum(eig_Sigma$vectors[, 1] * eig_recovered$vectors[, 1])))

cat("\nAbsolute alignment between original PC1 and Sigma_B PC1:\n")
print(abs(sum(eig_Sigma$vectors[, 1] * eig_B$vectors[, 1])))

cat("\n============================================================\n")
cat("STATIC SPECTRAL ATTRIBUTION ALONG FULL PCA EIGENVECTORS\n")
cat("============================================================\n")

static_attr <- data.frame(
  eigen_index = integer(),
  lambda = numeric(),
  between = numeric(),
  within = numeric(),
  between_plus_within = numeric(),
  recovery_error = numeric()
)

for (i in 1:2) {
  v <- eig_Sigma$vectors[, i]
  lambda_i <- eig_Sigma$values[i]

  between_i <- as.numeric(t(v) %*% Sigma_B %*% v)
  within_i <- as.numeric(t(v) %*% Sigma_W %*% v)

  static_attr <- rbind(
    static_attr,
    data.frame(
      eigen_index = i,
      lambda = lambda_i,
      between = between_i,
      within = within_i,
      between_plus_within = between_i + within_i,
      recovery_error = abs(lambda_i - between_i - within_i)
    )
  )
}

print(static_attr, digits = 10)

# Per-quadrant contribution to lambda1
v1_static <- eig_Sigma$vectors[, 1]

quad_lambda1_contrib <- data.frame(
  quadrant = character(),
  p = numeric(),
  between_contrib_to_lambda1 = numeric(),
  within_contrib_to_lambda1 = numeric(),
  total_contrib_to_lambda1 = numeric(),
  stringsAsFactors = FALSE
)

for (q in levels_Q) {
  p_q <- as.numeric(p_quad[q])
  u_q <- as.numeric(quad_means[, q] - mu)
  Sigma_q <- quad_covs[[q]]

  b_contrib <- p_q * sum(v1_static * u_q)^2
  w_contrib <- p_q * as.numeric(t(v1_static) %*% Sigma_q %*% v1_static)

  quad_lambda1_contrib <- rbind(
    quad_lambda1_contrib,
    data.frame(
      quadrant = q,
      p = p_q,
      between_contrib_to_lambda1 = b_contrib,
      within_contrib_to_lambda1 = w_contrib,
      total_contrib_to_lambda1 = b_contrib + w_contrib,
      stringsAsFactors = FALSE
    )
  )
}

cat("\nPer-quadrant contribution to classical lambda1:\n")
print(
  quad_lambda1_contrib[order(-quad_lambda1_contrib$total_contrib_to_lambda1), ],
  digits = 10
)

cat("\nCheck sum of per-quadrant contributions to lambda1:\n")
print(sum(quad_lambda1_contrib$total_contrib_to_lambda1))
cat("Classical lambda1:\n")
print(eig_Sigma$values[1])

# -----------------------------------------------------------------------------
# 5. Observable transition matrix
# -----------------------------------------------------------------------------

Q_current <- Q[-n]
Q_next <- Q[-1]

trans_counts <- table(Q_current, Q_next)

# Use unrounded transition probabilities for calculations.
trans_prob <- prop.table(trans_counts, margin = 1)

# Rounded matrix only for display.
trans_prob_print <- round(trans_prob, 4)

cat("\n============================================================\n")
cat("OBSERVABLE DIRECTIONAL MARKOV TRANSITION MATRIX\n")
cat("============================================================\n")

cat("\nTransition counts:\n")
print(trans_counts)

cat("\nTransition probabilities P(Q_{t+1} = q_next | Q_t = q_current):\n")
print(trans_prob_print)

cat("\nCrash persistence, CLPM to CLPM:\n")
print(trans_prob["CLPM", "CLPM"])

cat("\nCrash-to-rally reversal, CLPM to CUPM:\n")
print(trans_prob["CLPM", "CUPM"])

cat("\nRally-to-crash reversal, CUPM to CLPM:\n")
print(trans_prob["CUPM", "CLPM"])

# -----------------------------------------------------------------------------
# 6. One-step predictive mixture given current quadrant
# -----------------------------------------------------------------------------

forecast_mean <- function(q) {
  row <- as.numeric(trans_prob[q, ])
  names(row) <- levels_Q

  out <- rep(0, 2)
  names(out) <- colnames(Z)

  for (qq in levels_Q) {
    out <- out + row[qq] * quad_means[, qq]
  }

  out
}

forecast_cov <- function(q) {
  row <- as.numeric(trans_prob[q, ])
  names(row) <- levels_Q

  mu_fc <- forecast_mean(q)
  cov_sum <- matrix(0, 2, 2)

  for (qq in levels_Q) {
    m_qq <- quad_means[, qq]
    Sigma_qq <- quad_covs[[qq]]
    d <- as.numeric(m_qq - mu_fc)

    cov_sum <- cov_sum + row[qq] * (Sigma_qq + tcrossprod(d))
  }

  colnames(cov_sum) <- rownames(cov_sum) <- colnames(Z)
  cov_sum
}

cat("\n============================================================\n")
cat("ONE-STEP FORECAST FROM CURRENT OBSERVABLE QUADRANT\n")
cat("============================================================\n")

current_q <- "CLPM"

cat("\nForecast mean E[Z_{t+1} | Q_t = CLPM]:\n")
print(forecast_mean(current_q))

cat("\nForecast covariance Cov(Z_{t+1} | Q_t = CLPM):\n")
cat_matrix(forecast_cov(current_q), digits = 10)

# -----------------------------------------------------------------------------
# 7. Dynamic transition-path covariance decomposition
# -----------------------------------------------------------------------------

Z_lead <- Z[-1, , drop = FALSE]

mu_lead <- colMeans(Z_lead)
Sigma_lead <- pop_cov(Z_lead)

Sigma_B_dyn <- matrix(0, 2, 2)
Sigma_W_dyn <- matrix(0, 2, 2)

path_table <- data.frame(
  path = character(),
  q_current = character(),
  q_next = character(),
  p_path = numeric(),
  mean_X_lead = numeric(),
  mean_Y_lead = numeric(),
  u_X = numeric(),
  u_Y = numeric(),
  ell_rank1 = numeric(),
  max_abs_Bu_minus_ell_u = numeric(),
  alignment_with_Bpath_eigenvector = numeric(),
  stringsAsFactors = FALSE
)

path_covs <- list()

for (q in levels_Q) {
  for (qq in levels_Q) {
    idx <- which(Q_current == q & Q_next == qq)
    if (length(idx) == 0) next

    Z_path <- Z_lead[idx, , drop = FALSE]

    p_path <- nrow(Z_path) / nrow(Z_lead)
    m_path <- colMeans(Z_path)
    u_path <- as.numeric(m_path - mu_lead)
    Sigma_path <- pop_cov(Z_path)

    B_path <- p_path * tcrossprod(u_path)
    W_path <- p_path * Sigma_path

    Sigma_B_dyn <- Sigma_B_dyn + B_path
    Sigma_W_dyn <- Sigma_W_dyn + W_path

    ell_path <- p_path * sum(u_path^2)

    lhs <- as.numeric(B_path %*% u_path)
    rhs <- ell_path * u_path
    max_err <- max(abs(lhs - rhs))

    eig_Bpath <- eigen(B_path, symmetric = TRUE)
    v_Bpath <- eig_Bpath$vectors[, 1]

    if (sum(u_path^2) > 0) {
      v_u <- u_path / sqrt(sum(u_path^2))
      align <- abs(sum(v_u * v_Bpath))
    } else {
      align <- NA_real_
    }

    path_name <- paste(q, qq, sep = "->")
    path_covs[[path_name]] <- Sigma_path

    path_table <- rbind(
      path_table,
      data.frame(
        path = path_name,
        q_current = q,
        q_next = qq,
        p_path = p_path,
        mean_X_lead = m_path[1],
        mean_Y_lead = m_path[2],
        u_X = u_path[1],
        u_Y = u_path[2],
        ell_rank1 = ell_path,
        max_abs_Bu_minus_ell_u = max_err,
        alignment_with_Bpath_eigenvector = align,
        stringsAsFactors = FALSE
      )
    )
  }
}

Sigma_lead_recovered <- Sigma_B_dyn + Sigma_W_dyn

cat("\n============================================================\n")
cat("DYNAMIC TRANSITION-PATH COVARIANCE DECOMPOSITION\n")
cat("============================================================\n")

cat("\nLead-sample mean mu_lead:\n")
print(mu_lead)

cat("\nOriginal lead covariance Sigma_lead:\n")
cat_matrix(Sigma_lead)

cat("\nRecovered dynamic covariance Sigma_B_dyn + Sigma_W_dyn:\n")
cat_matrix(Sigma_lead_recovered)

cat("\nMax absolute dynamic recovery error:\n")
print(max(abs(Sigma_lead - Sigma_lead_recovered)))

cat("\nDynamic between-transition covariance Sigma_B_dyn:\n")
cat_matrix(Sigma_B_dyn)

cat("\nDynamic within-transition residual covariance Sigma_W_dyn:\n")
cat_matrix(Sigma_W_dyn)

cat("\nTransition-path rank-one primitive checks:\n")
print(path_table[order(-path_table$ell_rank1), ], digits = 8)

# -----------------------------------------------------------------------------
# 8. Dynamic eigensystem and transition-path spectral attribution
# -----------------------------------------------------------------------------

eig_lead <- eigen(Sigma_lead, symmetric = TRUE)
eig_dyn_recovered <- eigen(Sigma_lead_recovered, symmetric = TRUE)
eig_B_dyn <- eigen(Sigma_B_dyn, symmetric = TRUE)

cat("\n============================================================\n")
cat("DYNAMIC EIGENVALUE AND EIGENVECTOR RECOVERY\n")
cat("============================================================\n")

cat("\nEigenvalues of original Sigma_lead:\n")
print(eig_lead$values)

cat("\nEigenvalues of recovered dynamic covariance:\n")
print(eig_dyn_recovered$values)

cat("\nEigenvalues of dynamic between-transition covariance Sigma_B_dyn:\n")
print(eig_B_dyn$values)

cat("\nLeading eigenvector of Sigma_lead:\n")
print(eig_lead$vectors[, 1])

cat("\nLeading eigenvector of recovered dynamic covariance:\n")
print(eig_dyn_recovered$vectors[, 1])

cat("\nLeading eigenvector of Sigma_B_dyn only:\n")
print(eig_B_dyn$vectors[, 1])

cat("\nAbsolute alignment between Sigma_lead PC1 and recovered dynamic PC1:\n")
print(abs(sum(eig_lead$vectors[, 1] * eig_dyn_recovered$vectors[, 1])))

cat("\nAbsolute alignment between Sigma_lead PC1 and Sigma_B_dyn PC1:\n")
print(abs(sum(eig_lead$vectors[, 1] * eig_B_dyn$vectors[, 1])))

cat("\n============================================================\n")
cat("DYNAMIC SPECTRAL ATTRIBUTION ALONG FULL LEAD PCA EIGENVECTORS\n")
cat("============================================================\n")

dynamic_attr <- data.frame(
  eigen_index = integer(),
  lambda = numeric(),
  between_transition = numeric(),
  within_transition = numeric(),
  between_plus_within = numeric(),
  recovery_error = numeric()
)

for (i in 1:2) {
  v <- eig_lead$vectors[, i]
  lambda_i <- eig_lead$values[i]

  between_i <- as.numeric(t(v) %*% Sigma_B_dyn %*% v)
  within_i <- as.numeric(t(v) %*% Sigma_W_dyn %*% v)

  dynamic_attr <- rbind(
    dynamic_attr,
    data.frame(
      eigen_index = i,
      lambda = lambda_i,
      between_transition = between_i,
      within_transition = within_i,
      between_plus_within = between_i + within_i,
      recovery_error = abs(lambda_i - between_i - within_i)
    )
  )
}

print(dynamic_attr, digits = 10)

# Per-transition-path contribution to lead lambda1
v1_dyn <- eig_lead$vectors[, 1]

path_lambda1_contrib <- data.frame(
  path = character(),
  p_path = numeric(),
  between_contrib_to_lambda1 = numeric(),
  within_contrib_to_lambda1 = numeric(),
  total_contrib_to_lambda1 = numeric(),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(path_table))) {
  path_name <- path_table$path[i]
  p_path <- path_table$p_path[i]
  u_path <- c(path_table$u_X[i], path_table$u_Y[i])
  Sigma_path <- path_covs[[path_name]]

  b_contrib <- p_path * sum(v1_dyn * u_path)^2
  w_contrib <- p_path * as.numeric(t(v1_dyn) %*% Sigma_path %*% v1_dyn)

  path_lambda1_contrib <- rbind(
    path_lambda1_contrib,
    data.frame(
      path = path_name,
      p_path = p_path,
      between_contrib_to_lambda1 = b_contrib,
      within_contrib_to_lambda1 = w_contrib,
      total_contrib_to_lambda1 = b_contrib + w_contrib,
      stringsAsFactors = FALSE
    )
  )
}

cat("\nTop transition-path contributions to lead lambda1:\n")
print(
  path_lambda1_contrib[order(-path_lambda1_contrib$total_contrib_to_lambda1), ],
  digits = 10
)

cat("\nTop transition-path BETWEEN contributions to lead lambda1:\n")
print(
  path_lambda1_contrib[order(-path_lambda1_contrib$between_contrib_to_lambda1), ],
  digits = 10
)

cat("\nCheck sum of path contributions to lead lambda1:\n")
print(sum(path_lambda1_contrib$total_contrib_to_lambda1))
cat("Lead lambda1:\n")
print(eig_lead$values[1])

# -----------------------------------------------------------------------------
# 9. Conceptual contrast with classical HMM
# -----------------------------------------------------------------------------

cat("\n============================================================\n")
cat("CONCEPTUAL CONTRAST WITH CLASSICAL HMM\n")
cat("============================================================\n")

cat("\nClassical HMM workflow:\n")
cat("  1. Choose number of hidden states.\n")
cat("  2. Specify emission distributions.\n")
cat("  3. Estimate hidden states and parameters, usually by EM.\n")
cat("  4. Interpret the latent states after estimation.\n")

cat("\nNNS directional Markov regime workflow:\n")
cat("  1. Define observable states by directional quadrant membership.\n")
cat("  2. Estimate state probabilities by frequencies.\n")
cat("  3. Estimate transition probabilities by counts.\n")
cat("  4. Estimate state means, state covariances, and spectral contributions directly.\n")

cat("\nBottom line:\n")
cat("  HMMs infer hidden regimes and then interpret them.\n")
cat("  NNS defines interpretable directional regimes first and then measures their dynamics.\n")

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")

```
