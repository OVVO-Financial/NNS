# Numerical Example: Degree-0 Bayes and Degree-1 Hinge-Surface Recovery

This example illustrates the distinction between:

1. **Degree-0 Bayes**, which gives the exact conditional event probability directly, and
2. **Degree-1 hinge-surface recovery**, which reconstructs the same event probability indirectly through recovery of the joint CDF.

The point is not that degree 1 replaces degree 0 as a probability ratio. Rather, degree 1 supplies the hinge surface from which the joint law can be recovered.

## Setup

We simulate a dependent bivariate sample:

```r
library(NNS)

set.seed(123)
n <- 2000
x <- rnorm(n)
y <- rnorm(n) + 0.8 * x

t_x <- 0
t_y <- 0
```

Empirical sample quantities:

- $P(X>0) = 0.508$
- $P(Y>0) = 0.500$
- $P(X>0, Y>0) = 0.362$
- $P(Y>0 \mid X>0) = 0.7125984$

---

## 1. Degree-0 Bayes

At degree 0, the NNS operators coincide with quadrant probabilities:

$$ \mathrm{Co.UPM}(0,x,y;t_x,t_y)=P(X>t_x,\;Y>t_y), \qquad \mathrm{UPM}(0,t_x,x)=P(X>t_x). $$

So the conditional event probability is obtained exactly as

$$ P(Y>t_y \mid X>t_x) = \frac{\mathrm{Co.UPM}(0,x,y;t_x,t_y)}{\mathrm{UPM}(0,t_x,x)}. $$

Using $t_x=t_y=0$:

```r
joint_prob_deg0 <- Co.UPM(0, x, y, target_x = t_x, target_y = t_y)
marg_prob_x_deg0 <- UPM(0, t_x, x)

cond_eventprob_deg0 <- joint_prob_deg0 / marg_prob_x_deg0
cond_eventprob_deg0
```

Output:

```r
[1] 0.7125984
```

This matches the empirical conditional event probability exactly in-sample.

---

## 2. Degree-1 hinge-surface recovery

Define the raw lower hinge surface

$$
H(t_x,t_y)=E[(t_x-X)_+(t_y-Y)_+].
$$

The degree-1 recovery theorem states that

$$
\frac{\partial^2 H}{\partial t_x \partial t_y}(t_x,t_y)=F_{X,Y}(t_x,t_y).
$$

So we can recover the joint CDF at $(0,0)$ numerically by mixed finite differences.

```r
targets <- seq(-3, 3, length.out = 61)
h <- diff(targets)[1]

hinge_surface <- outer(
  targets,
  targets,
  Vectorize(function(tx, ty) {
    Co.LPM_nD(
      data = cbind(x, y),
      target = c(tx, ty),
      degree = 1,
      norm = FALSE
    )
  })
)

i0 <- which.min(abs(targets - 0))

joint_cdf_from_hinge <- (
  hinge_surface[i0 + 1, i0 + 1] -
  hinge_surface[i0 + 1, i0] -
  hinge_surface[i0,     i0 + 1] +
  hinge_surface[i0,     i0]
) / h^2
```

Recovered and empirical lower-left quadrant probabilities:

- Recovered $F(0,0)$ from hinge surface: `0.372003`
- Empirical $P(X \le 0, Y \le 0)$: `0.354000`

---

## 3. Reconstructing the same upper-right event probability

From the recovered joint CDF,

$$
P(X>0,Y>0)=1-F_X(0)-F_Y(0)+F_{X,Y}(0,0).
$$

Using degree-0 marginal CDF values at zero:

```r
FX0 <- LPM.ratio(0, t_x, x)
FY0 <- LPM.ratio(0, t_y, y)

joint_eventprob_from_hinge <- 1 - FX0 - FY0 + joint_cdf_from_hinge
marg_eventprob_x <- 1 - FX0

cond_eventprob_from_hinge <- joint_eventprob_from_hinge / marg_eventprob_x
cond_eventprob_from_hinge
```

Output:

```r
[1] 0.7480377
```

Comparison:

- Recovered $P(X>0, Y>0)$: `0.380003`
- Empirical $P(X>0, Y>0)$: `0.362000`
- Recovered $P(Y>0 \mid X>0)$: `0.7480377`
- Empirical $P(Y>0 \mid X>0)$: `0.7125984`

---

## 4. Interpretation

This example shows two different routes to the same event-level quantity $P(Y>0 \mid X>0)$:

- **Degree 0** gives it directly and exactly through quadrant probabilities.
- **Degree 1** gives it indirectly by recovering the joint CDF from the raw hinge surface and then applying inclusion-exclusion.

The degree-1 reconstruction is close, but not exact, because the mixed derivative is approximated numerically on a finite grid. The discrepancy is therefore numerical, not conceptual.

More importantly, degree 1 should be understood structurally:

- mixed second derivatives recover the joint CDF,
- further differentiation recovers the joint density when it exists,
- and only then does the full continuous analogue of Bayes arise through density ratios.

Thus:

- **degree 0** is the exact event-probability layer,
- **degree 1** is the law-recovery layer.

---

## 5. Summary table

| Quantity | Value |
|---|---:|
| Exact degree-0 conditional event probability | 0.7125984 |
| Degree-1 hinge-surface reconstructed event probability | 0.7480377 |
| Empirical event probability truth | 0.7125984 |

These results confirm that degree 0 and degree 1 do not play the same role. Degree 0 yields Bayes directly. Degree 1 recovers the joint law from which event probabilities, and ultimately continuous conditional densities, can be constructed.
