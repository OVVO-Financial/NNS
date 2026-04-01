# `NNS.diff` for Noisy Gradients

## Executive Summary

`NNS.diff` is a core numerical differentiation routine in the NNS package. It determines the derivative of a univariate function by first solving for an appropriate perturbation scale `h` using a geometric procedure based on projected secant lines onto the y-axis. These projected points are then used to infer a finite step size `h` that is fed into multiple derivative estimators.

The method works as follows:

1. Evaluate the function at the target `point` and at `point ± h` (initial `h` defaults to 0.1).
2. Compute the y-intercepts (`B₁` and `B₂`) of the two secant lines connecting those points.
3. Treat the interval between the two intercepts as a bracket and iteratively bisect it.
4. At each candidate midpoint `B`, solve via `uniroot` for the step size `h*` that would make the left-hand secant project exactly to that `B`.
5. The sign of `h*` tells the algorithm which side of the bracket to shrink. The process repeats until either `|h*| < tol` or `max.iter` is reached.

The final inferred `h*`, whether fully converged or stall-limited by noise, is then supplied to multiple downstream estimators. For **analytic, complex-compatible functions**, the strongest downstream output is the **Complex Step Derivative (Inferred h)**. For **non-analytic, piecewise, thresholded, saturated, quantized, or real-only black-box functions**, the core projected derivative itself, `DERIVATIVE`, becomes the most important output.

So `NNS.diff` has a two-regime interpretation:

- **analytic regime:** use the inferred step to drive the complex-step row
- **black-box real-only regime:** use the projected-secant derivative as the all-terrain estimate

**For analytic, complex-compatible functions, the strongest workflow is therefore:**

```r
out <- NNS.diff(
  f = f,
  point = x0,
  h = abs(point) * 0.1 + 0.01,
  tol = 1e-10,
  max.iter = 1000,
  digits = 12,
  print.trace = FALSE,
  plot = FALSE
)

grad <- out["Complex Step Derivative (Inferred h)", 1]
```

**For non-analytic or real-only black-box functions, the preferred workflow is instead:**

```r
out <- NNS.diff(
  f = f,
  point = x0,
  h = abs(point) * 0.1 + 0.01,
  tol = 1e-10,
  max.iter = 1000,
  digits = 12,
  print.trace = FALSE,
  plot = FALSE
)

grad <- out["DERIVATIVE", 1]
```

`NNS.diff` does **not** simply return “a derivative.” Its real contribution is that it first solves the classic perturbation-scale problem geometrically, then hands that stable scale to a derivative formula that is far less numerically fragile than ordinary finite differences, or, when complex inputs are invalid, it returns a projected derivative that remains stable and directionally useful on real-only surfaces.

---

## 1. Where `NNS.diff` sits inside NNS

`NNS.diff` is not an isolated add-on. In the official manual it appears as a named core routine within the NNS package, alongside other NNS procedures such as `NNS.dep`, `NNS.reg`, `NNS.ARMA`, and `NNS.copula`. The manual indexes `NNS.diff` as its own documented entry and defines it as **“NNS Numerical Differentiation.”**

The documented description is specific:

> “Determines numerical derivative of a given univariate function using projected secant lines on the y-axis. These projected points infer finite steps `h`, in the finite step method.”

The original interface was:

```r
NNS.diff(f, point, h = 0.1, tol = 1e-10, digits = 12, print.trace = FALSE)
```

The current implementation adds `max.iter` and `plot`, and returns a matrix containing:
- the projected derivative
- the inferred step size
- finite-step estimates at the initial and inferred `h`
- the complex-step estimate at the inferred `h`
- convergence diagnostics
- the termination code

That makes `NNS.diff` better understood as a **small diagnostic framework for local differentiation**, not just a single derivative formula.

---

## 2. The innovative idea: geometry before differencing

Most numerical differentiation methods start by **assuming** you already know a good perturbation size `h`.

### Standard finite differences
They begin with

$$
\hat{f}'(x;h)=\frac{f(x+h)-f(x-h)}{2h}
$$

and then force the user to decide what `h` should be.

### Richardson extrapolation
Richardson begins with finite differences too, then tries to cancel truncation error by shrinking `h` repeatedly.

### `NNS.diff`
`NNS.diff` **flips the order**.

It treats the choice of `h` as the primary unknown and solves for it **geometrically** before any derivative formula is applied:

1. Evaluate `f` at `point`, `point - h`, and `point + h`.
2. Compute the y-intercepts of the two secant lines: $B_1 = f(x) - \frac{f(x) - f(x-h)}{h} \cdot x$ and $B_2 = f(x) - \frac{f(x+h) - f(x)}{h} \cdot x$.
3. Bracket the interval `[min(B_1, B_2), max(B_1, B_2)]`.
4. For a candidate midpoint `B`, solve $B = f(x) - \frac{f(x) - f(x-h^*)}{h^*} \cdot x$ for `h*` via `uniroot`.
5. Use the sign of `h*` to decide which half of the bracket to discard.
6. Repeat until either `|h*| < tol` or `max.iter` is reached.

In clean analytic functions the search converges quickly to a tiny `h*`, often around $10^{-11}$. In noisy functions it naturally stalls at a **larger** `h*`, automatically moving away from the tiny steps that would amplify noise.

This is the key innovation: `NNS.diff` converts the step-size problem into a **stable geometric bisection**.

<img src="/examples/secants.png"  style="border: none; outline: none; margin: 0; padding: 0; display: block;"/>


## 3. What the updated implementation actually does

The current implementation adds several practical improvements.

### A. Explicit bisection with termination diagnostics
- `termination.code = 0` means clean convergence to tolerance
- `termination.code = 1` means the search hit `max.iter` and stalled in a noisy regime
- `termination.code = 2` means `uniroot` failed

### B. Multiple downstream estimators using the same inferred `h`
The returned matrix includes:
- Projected secant derivative (`DERIVATIVE`)
- Initial finite-step rows
- Inferred finite-step rows
- **Complex Step Derivative (Inferred h)**

### C. Noise-aware defaults
`max.iter` provides a hard upper bound so noisy functions cannot hang indefinitely.

### D. Correct handling of locally linear real-only regions
A crucial update is the fix for the case where the initial projected intercepts are identical, $B_1 = B_2$. That situation does **not** mean the derivative fails to exist. It usually means the local slope has already been identified exactly.

That correction matters for:
- `abs(x)` away from the kink
- `ReLU` away from the threshold
- clipped linear interiors

After the fix, `NNS.diff` correctly returns the local common secant slope instead of incorrectly treating the case as a failure.

This turns out to be essential for revealing the full value of `NNS_Proj` on piecewise and black-box surfaces.

---

## 4. Why this is preferred over standard methods for analytic functions

For analytic, complex-compatible functions, the preferred output is:

$$
\texttt{NNS.diff} \;\rightarrow\; \texttt{"Complex Step Derivative (Inferred h)"}
$$

because it combines the strongest feature of two different ideas.

### Existing methods leave one key problem unsolved

#### Fixed finite differences
They require the user to choose `h`. In noise, that is exactly the hard part. If `h` is too small, the variance grows like

$$
\mathrm{Var}(\hat f'(x;h)) \propto \frac{\sigma^2}{h^2}.
$$

If `h` is too large, truncation bias dominates.

#### Richardson extrapolation
Richardson is excellent in clean smooth settings because it shrinks `h` and cancels truncation error asymptotically. But in noisy settings, repeated refinement can amplify noise instead of helping.

#### Complex step alone
Complex step is numerically excellent because it avoids the subtractive cancellation that plagues real-axis differences. But in practice you still need a **useful perturbation scale**.

### What `NNS.diff -> Complex Step (Inferred h)` adds

This workflow solves both pieces at once:

1. `NNS.diff` finds a **geometry-informed perturbation scale**
2. complex step uses that scale with a derivative formula that avoids the real-axis subtraction problem

It is not that `NNS.diff` somehow replaces complex step. It is that `NNS.diff` supplies the missing piece that complex step usually does not address directly: **how to choose a practical perturbation scale without manual tuning.**

---

## 5. What problem this solves in practice

For analytic functions, this combination solves several concrete problems at once.

### Problem 1: manual step-size selection
Most users do not know the right `h` ahead of time. `NNS.diff` turns step-size selection into an inferred quantity rather than a hand-tuned hyperparameter.

### Problem 2: cancellation fragility in finite differences
Real-valued finite differences subtract nearby evaluations. Complex step avoids that subtraction route.

### Problem 3: Richardson’s weakness in noise
Richardson improves clean asymptotics but is not designed around noisy stability. `NNS.diff` instead searches for a practical scale region and explicitly exposes when the search has stalled by hitting `max.iter`.

### Problem 4: low-configuration gradient estimation
This is useful when you want a gradient estimate from a simulator or algorithmic model, but you do not want to tune a finite-difference step by hand, run a full step sweep, or trust aggressive extrapolation in noise.

---

## 6. Benchmark setup recap: analytic functions

We compared:

- `NNS_Proj`: projected secant derivative from `NNS.diff`
- `NNS_FinInit`: averaged finite step at the initial `h`
- `NNS_FinInf`: averaged finite step at the inferred `h`
- `NNS_CplxInf`: complex-step derivative at the inferred `h`
- `Richardson`
- `OracleFD`: oracle centered finite-difference sweep over an external `h` grid

Functions tested:

- `sin(x)` at `x = 1`
- `exp(x)` at `x = 1`
- `x^3` at `x = 2`

Noise levels:

- `sigma = 0`
- `1e-4`
- `1e-3`
- `1e-2`

Replications: `R = 100`

Important implementation detail: under noisy benchmarking, we used **frozen-noise surfaces per replication**, so `NNS.diff` sees a deterministic perturbed surface inside each call rather than a moving stochastic target.

---

## 7. Results table: analytic benchmark (`max.iter = 1000`)

| Function | Sigma | NNS Proj RRMSE | NNS FinInit RRMSE | NNS FinInf RRMSE | NNS CplxInf RRMSE | Richardson RRMSE | OracleFD RRMSE | Median NNS h | Median Iterations | OracleFD h | Termination |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| sin(x) at x=1 | 0 | 3.66e-07 | 1.67e-03 | 1.28e-06 | 2.59e-13 | 1.71e-14 | 1.65e-10 | 2.20e-11 | 30 | 3.16e-05 | 0 |
| sin(x) at x=1 | 1e-4 | 8.87e-03 | 2.16e-03 | 7.12e-01 | 3.88e-05 | 6.25e-02 | 2.06e-03 | 8.35e-03 | 1000 | 8.88e-02 | 1 |
| sin(x) at x=1 | 1e-3 | 3.14e-02 | 1.38e-02 | 7.78e+00 | 3.03e-04 | 5.90e-01 | 9.19e-03 | 2.01e-02 | 1000 | 1.68e-01 | 1 |
| sin(x) at x=1 | 1e-2 | 1.72e-01 | 1.43e-01 | 6.14e+01 | 1.63e-02 | 6.30e+00 | 3.74e-02 | 5.20e-02 | 1000 | 3.16e-01 | 1 |
| exp(x) at x=1 | 0 | 1.64e-07 | 1.67e-03 | 3.52e-07 | 1.67e-14 | 1.19e-14 | 1.70e-10 | 3.20e-11 | 29 | 3.16e-05 | 0 |
| exp(x) at x=1 | 1e-4 | 4.96e-03 | 1.70e-03 | 1.47e-01 | 1.71e-05 | 1.24e-02 | 6.96e-04 | 3.82e-03 | 1000 | 4.70e-02 | 1 |
| exp(x) at x=1 | 1e-3 | 9.22e-03 | 2.84e-03 | 5.39e-01 | 1.27e-04 | 1.17e-01 | 3.03e-03 | 1.37e-02 | 1000 | 8.88e-02 | 1 |
| exp(x) at x=1 | 1e-2 | 3.53e-02 | 2.86e-02 | 6.34e+00 | 1.58e-03 | 1.25e+00 | 1.48e-02 | 4.88e-02 | 1000 | 2.30e-01 | 1 |
| x^3 at x=2 | 0 | 7.43e-06 | 8.33e-04 | 9.71e-07 | 0.00e+00 | 5.34e-14 | 8.59e-11 | 8.60e-11 | 30 | 3.16e-05 | 0 |
| x^3 at x=2 | 1e-4 | 2.37e-03 | 8.37e-04 | 5.68e-02 | 1.58e-06 | 2.82e-03 | 2.19e-04 | 1.77e-03 | 1000 | 3.42e-02 | 1 |
| x^3 at x=2 | 1e-3 | 4.89e-03 | 9.25e-04 | 2.86e-01 | 1.55e-05 | 2.66e-02 | 9.00e-04 | 7.64e-03 | 1000 | 6.46e-02 | 1 |
| x^3 at x=2 | 1e-2 | 1.40e-02 | 6.55e-03 | 4.64e+00 | 1.44e-04 | 2.84e-01 | 3.85e-03 | 2.33e-02 | 1000 | 1.68e-01 | 1 |

---

## 8. What the analytic results show

### A. In clean settings, Richardson still wins
At `sigma = 0`, Richardson is best, exactly as one expects in clean smooth settings.

### B. In noisy settings, the inferred scale is valuable
The median inferred `h` increases as noise increases. That is what a sensible scale selector should do.

### C. The complex-step row is the standout result
`NNS_CplxInf` remains extremely accurate across the analytic examples, even when the projected derivative and Richardson worsen under noise.

So the benchmark supports a very specific conclusion:

> the most useful output of `NNS.diff` for analytic noisy functions is not necessarily the projected derivative itself, but the inferred step scale fed into the complex-step formula.

### D. Increasing `max.iter` does not materially change the answer
When `max.iter` increases from `100` to `1000`, the noisy results barely change while median iterations rise to `1000`. That means the search is **stall-limited**, not under-iterated.

---

## 9. Extending the story: why `NNS.diff` also matters for non-analytic and black-box functions

The analytic conclusion is only half the story.

The complex-step row is only valid when the function is analytic and accepts complex perturbations. But many practical objectives are not of that form. They may be:

- piecewise continuous
- thresholded
- clipped or saturated
- quantized
- rounded
- Monte Carlo based
- real-only black-box systems

In those cases, the complex-step row is invalid or meaningless. This is where the projected derivative itself, `NNS_Proj`, becomes the key output.

Rather than treating this as a fallback of lesser importance, the updated experiments show that `NNS_Proj` solves a **different problem**:

> how to obtain a stable, directionally meaningful local gradient proxy on a noisy real-only surface where complex perturbations are impossible or inappropriate.

That is why `NNS_Proj` is best described as the **all-terrain estimator**.

---

## 10. Black-box benchmark setup

To reveal that value, a second benchmark focused on non-analytic and black-box-like functions.

Categories included:

- **piecewise continuous**
  - `abs(x)` smooth side
  - `abs(x)` at kink
  - `ReLU` smooth side
  - `ReLU` at kink
- **saturation**
  - clipped linear interior
  - clipped linear threshold
- **discontinuous**
  - indicator threshold
- **quantized**
  - rounded quadratic
- **black-box payoff-like**
  - Monte Carlo option payoff style

For these functions, derivative existence at the evaluation point is not always the right evaluation criterion. So the benchmark also tracked:

- sign accuracy
- tangent prediction error
- estimator standard deviation
- iteration and termination diagnostics

That is important, because for black-box systems a **stable local directional estimate** is often more useful than a classical smooth derivative in the textbook sense.

---

## 11. Black-box results: grouped median summary

| Category | Sigma | NNS_Proj_RRMSE | NNS_FinInit_RRMSE | NNS_FinInf_RRMSE | Richardson_RRMSE | OracleFD_RRMSE | NNS_Proj_SignAcc | NNS_FinInit_SignAcc | NNS_FinInf_SignAcc | Richardson_SignAcc | OracleFD_SignAcc | NNS_Proj_TangentMAE | NNS_FinInit_TangentMAE | NNS_FinInf_TangentMAE | Richardson_TangentMAE | OracleFD_TangentMAE | Median_NNS_h | Median_NNS_Iterations | Median_Termination_Code |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| piecewise_continuous | 1e-4 | 1.0367e-03 | 7.4347e-04 | 3.5924e-01 | 3.3314e-02 | 2.9618e-04 | 1.00 | 1.00 | 0.9773 | 1.00 | 1.00 | 1.1129e-04 | 1.1303e-04 | 1.0471e-03 | 2.9742e-04 | 1.1356e-04 | 5.0932e-02 | 1000 | 1 |
| saturation | 1e-4 | 7.4398e-04 | 6.6351e-04 | 1.1489e+01 | 2.9752e-02 | 2.0921e-04 | 1.00 | 1.00 | 0.9701 | 1.00 | 1.00 | 1.1972e-04 | 1.1997e-04 | 1.6204e-02 | 2.4907e-04 | 1.1922e-04 | 4.4941e-02 | 1000 | 1 |
| piecewise_continuous | 1e-3 | 9.5846e-03 | 6.8462e-03 | 1.4168e+01 | 3.1223e-01 | 3.5651e-03 | 1.00 | 1.00 | 0.8647 | 1.00 | 1.00 | 1.1325e-03 | 1.1427e-03 | 3.7789e-02 | 2.4987e-03 | 1.1292e-03 | 3.1348e-02 | 1000 | 1 |
| saturation | 1e-3 | 8.0072e-03 | 7.0641e-03 | 8.6091e+00 | 3.4521e-01 | 2.3545e-03 | 1.00 | 1.00 | 0.9286 | 1.00 | 1.00 | 1.2044e-03 | 1.2107e-03 | 2.5865e-02 | 2.8792e-03 | 1.2091e-03 | 2.9158e-02 | 1000 | 1 |
| piecewise_continuous | 1e-2 | 7.7017e-02 | 7.5496e-02 | 4.9546e+01 | 3.1399e+00 | 2.9688e-02 | 1.00 | 1.00 | 0.8026 | 0.68 | 1.00 | 1.0507e-02 | 1.0556e-02 | 2.0036e-01 | 2.6432e-02 | 1.0441e-02 | 4.4360e-02 | 1000 | 1 |
| saturation | 1e-2 | 8.3284e-02 | 6.5639e-02 | 5.0264e+01 | 3.4355e+00 | 2.3065e-02 | 1.00 | 1.00 | 0.8718 | 0.57 | 1.00 | 1.2193e-02 | 1.2141e-02 | 1.7557e-01 | 3.0687e-02 | 1.2139e-02 | 3.7506e-02 | 1000 | 1 |

---

## 12. What the black-box results show

### A. `NNS_Proj` is exact on clean locally linear regions
After the implementation fix, the projected derivative is exact on:

- `abs(x)` smooth side
- `ReLU` smooth side
- clipped linear interior

at `sigma = 0`.

That is exactly the behavior one wants on locally affine black-box surfaces.

### B. `NNS_Proj` is highly directionally stable
Across the grouped non-analytic results:

- `NNS_Proj_SignAcc = 1.00`
- `NNS_FinInit_SignAcc = 1.00`
- `OracleFD_SignAcc = 1.00`

while Richardson degrades sharply at high noise.

For example, at `sigma = 1e-2`:

- piecewise continuous: `Richardson_SignAcc = 0.68`
- saturation: `Richardson_SignAcc = 0.57`

So `NNS_Proj` preserves local direction much more reliably than Richardson in noisy non-analytic settings.

### C. `NNS_Proj` is near-oracle on local tangent prediction
The tangent MAE columns are especially revealing.

For both piecewise and saturation categories, `NNS_Proj_TangentMAE` is extremely close to `OracleFD_TangentMAE`. In other words:

> `NNS_Proj` gives nearly oracle-quality local tangent prediction without requiring an external sweep over step sizes.

That is a major practical advantage.

### D. The inferred scale still matters, but the projected derivative is the useful output
In these black-box settings, the inferred `h` still adapts upward under noise, but the main value is not a complex-step row. The main value is that the geometric search stabilizes the local scale enough for the **projected derivative** to remain useful.

### E. Stall-limited does not mean useless
As in the analytic noisy case, median iterations frequently hit `1000` with termination code `1`.

That means the search is **stall-limited**, not fully converged. But the projected derivative remains highly usable. So the correct interpretation is:

> in noise, `NNS.diff` often settles into a practical scale region rather than converging to an infinitesimal one, and that practical scale is enough to support a stable projected estimate.

---

## 13. Why `NNS_Proj` is the all-terrain estimator

For black-box and non-analytic functions, `NNS_Proj` solves a different practical problem than the complex-step row.

### Complex-step row solves
How do I get an extremely accurate derivative for an analytic function once I have a good perturbation scale?

### `NNS_Proj` solves
How do I extract a stable local directional slope estimate from a noisy real-only system where complex inputs are invalid?

That is why `NNS_Proj` is all-terrain:

- it uses only real evaluations
- it does not require analytic continuation
- it behaves correctly on locally linear regions
- it remains directionally stable under noise
- it gives near-oracle local tangent prediction on black-box-like examples
- it avoids the noise fragility of Richardson in non-analytic settings

This makes it useful for:

- piecewise losses
- clipped controls
- threshold decision rules
- quantized simulators
- payoff functions
- Monte Carlo black-box objectives
- real-world physical systems that cannot accept complex perturbations

---

## 14. Practical recommendation by regime

### Analytic, complex-compatible functions
Use:

```r
out["Complex Step Derivative (Inferred h)", 1]
```

because this combines automatic geometry-based scale selection with a highly stable derivative formula.

### Non-analytic, piecewise, thresholded, quantized, or real-only black-box functions
Use:

```r
out["DERIVATIVE", 1]
```

because the projected derivative is the estimator that remains meaningful and robust when complex continuation is unavailable.

---

## 15. Why this is ultimately useful

The overall contribution of `NNS.diff` is broader than “it estimates derivatives.”

It gives a practical answer to a question that standard methods usually leave unresolved:

> **What local perturbation scale should I use to obtain a stable gradient estimate in noise?**

For analytic functions, that answer is:
- infer the scale geometrically
- then use the complex-step row

For black-box real-only functions, that answer is:
- infer the scale geometrically
- then rely on the projected derivative itself

This is useful because it avoids:
- hand-tuning `h`
- blind trust in tiny steps
- aggressive Richardson refinement in noisy settings
- full external step sweeps in routine workflows

---

## 16. Final caveat

The “standout winner” status of the complex-step row applies only when the function is **analytic and complex-compatible**.

Do **not** rely on the complex-step row when:

- `f(z)` is not meaningful for complex `z`
- the model is non-analytic or branch-heavy
- the function is piecewise, clipped, thresholded, rounded, or quantized in a way that breaks complex continuation
- the system is a real-world physical process with no meaningful complex perturbation interpretation

In those cases, the projected derivative is not merely a fallback. It is the correct general-purpose output.

---

## Bottom line

The historical and conceptual importance of `NNS.diff` is not that it is just another finite-difference routine.

Its innovation is that it treats numerical differentiation as a **geometry-driven scale-selection problem first**, and only then as a derivative-estimation problem.

That yields two distinct but complementary strengths:

- for **analytic noisy models**, `NNS.diff` is best used as an automatic inferred-step selector feeding  
  **`"Complex Step Derivative (Inferred h)"`**
- for **non-analytic, piecewise, thresholded, saturated, quantized, or black-box real-only models**, the most useful output is the projected derivative,  
  **`"DERIVATIVE"`**, that is, **`NNS_Proj`**

So the full conclusion is:

> `NNS.diff` is not just a derivative routine. It is a geometry-based local-scale inference method whose output should be chosen by regime. The complex-step row is the best endpoint for analytic functions, while `NNS_Proj` is the all-terrain estimator that makes the method broadly useful on real-world black-box surfaces.


### Experiments

- [Benchmarks](https://github.com/OVVO-Financial/NNS/blob/Data-and-Simulation-Routines/NNS-Simulation-Routines/Benchmarks%20for%20NNS_diff.R)
- [Non-Analytic Benchmarks](https://github.com/OVVO-Financial/NNS/blob/Data-and-Simulation-Routines/NNS-Simulation-Routines/Non-Analytic%20Benchmarks%20for%20NNS_diff.R)
