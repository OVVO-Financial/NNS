# Directional causation (book ch14), redesigned for #68.
#
# For a lag tau >= 0 the two directions are evaluated on explicitly lagged
# pairs:
#   x -> y on (X_{t-tau}, Y_t),   y -> x on (Y_{t-tau}, X_t).
# Each pair is mapped to the copula scale with per-column ranks (the degree-0
# lower partial moment, the empirical CDF), so the statistic is invariant to
# shifts, rescaling and any increasing transform of either variable. The lag
# supplies the temporal asymmetry; at tau = 0 the asymmetry comes from the
# directional measures themselves. The index measures directional lagged
# association: it does not condition on other variables.
#
# Each direction gets a raw score in [0, 1],
#   C~_{x -> y} = (P_{x -> y} + max(rho_{x -> y} - rho_{y -> x}, 0)) / 2,
# where
#   P_{x -> y}   share of Y_t's deviation about its mean resolved by
#                conditioning on the position of X_{t-tau} (equal-frequency
#                cells of the copula-scale cause). The deviation is the
#                degree-2 lower plus upper partial moment about the mean with
#                the m / (m - 1) factor: exactly unbiased under any
#                permutation, for any distribution (#96).
#   rho_{x -> y} the debiased directional dependence of Y_t on X_{t-tau}
#                (NNS.dep's cell measure on the lagged co-partial moments of
#                the pair), averaged over the levels of the partition on
#                X_{t-tau}. NNS.dep reports the strongest level; the average
#                has no selection step, so independent series are centred
#                whatever their margins (#96).
# The index is
#   C_{x -> y} = (C~_{x -> y} + k) / (C~_{x -> y} + C~_{y -> x} + 2k),
#   C_{y -> x} = 1 - C_{x -> y},
# with a neutral pseudo-count k = 1 / sqrt(n - tau), the order of the raw
# scores' sampling variation, so weak or absent evidence stays near 0.5.


# Copula scale: the degree-0 LPM of each observation against its own column.
.nns_caus_u <- function(v) LPM.ratio(0, v, v)


# Deviation of v about its mean: the degree-2 lower plus upper partial moment
# with the m / (m - 1) factor (the sample variance). As a U-statistic it has
# the same expectation in a random subset as in the whole sample, so the
# within-cell average below is unbiased under independence for any
# distribution, ties included -- no normal-theory small-sample factor. Degree
# 2 keeps the deviation about the cell's centre: knowing y = x^2 leaves x
# split between +-sqrt(y), which resolves none of x's deviation about its
# mean (a pairwise measure such as the Gini mean difference would credit it).
.nns_caus_dev <- function(v) {
  m <- length(v)
  if (m < 2L) return(NA_real_)
  mu <- mean(v)
  (LPM(2, mu, v) + UPM(2, mu, v)) * m / (m - 1)
}


# Share of the effect's deviation resolved by the cause's position.
# Both inputs are on the copula scale and have the same length.
.nns_caus_P <- function(cause_u, effect_u) {
  n <- length(cause_u)
  if (n < 3L) return(0)
  k <- max(2L, as.integer(floor(sqrt(n))))
  cell <- pmax(1L, as.integer(ceiling(cause_u * k)))
  total <- .nns_caus_dev(effect_u)
  if (!is.finite(total) || total <= 0) return(0)
  parts <- split(effect_u, cell)
  w <- vapply(parts, length, integer(1L))
  d <- vapply(parts, .nns_caus_dev, numeric(1L))
  ok <- is.finite(d)
  if (!any(ok)) return(0)
  within <- sum(w[ok] * d[ok]) / sum(w[ok])
  min(1, max(0, 1 - within / total))
}


# Directional dependence of b on a for the causal comparison: the unfloored
# debiased dependence of NNS.dep's cells, averaged over the levels of the
# recursive gravity split on a (same partition settings as NNS.dep).
.nns_caus_rho <- function(a, b) {
  l <- length(a)
  obs <- max(8L, as.integer(l / 8L))
  ord <- max(ceiling(log(l, 2)), 1L)
  z <- tryCatch({
    q <- NNS_part_cpp(a, b, type = "XONLY", order_in = as.integer(ord),
                      obs_req = as.integer(obs), min_obs_stop = FALSE,
                      noise_reduction = "off", quadrants_only = TRUE)$quadrant
    mean(NNS_dep_levels_cpp(as.numeric(a), as.numeric(b), as.character(q)))
  }, error = function(e) 0)
  if (is.finite(z)) z else 0
}


# Raw scores for both directions at lag tau.
.nns_caus_raw <- function(x, y, tau) {
  n <- length(x)
  lead <- (tau + 1L):n
  lag <- 1L:(n - tau)
  x_lag <- .nns_caus_u(x[lag]);  y_lead <- .nns_caus_u(y[lead])   # (X_{t-tau}, Y_t)
  y_lag <- .nns_caus_u(y[lag]);  x_lead <- .nns_caus_u(x[lead])   # (Y_{t-tau}, X_t)

  P_xy <- .nns_caus_P(x_lag, y_lead)
  P_yx <- .nns_caus_P(y_lag, x_lead)
  rho_xy <- .nns_caus_rho(x_lag, y_lead)
  rho_yx <- .nns_caus_rho(y_lag, x_lead)

  list(raw_xy = (P_xy + max(rho_xy - rho_yx, 0)) / 2,
       raw_yx = (P_yx + max(rho_yx - rho_xy, 0)) / 2,
       P = c(P_xy, P_yx), rho = c(rho_xy, rho_yx),
       n_pairs = n - tau,
       pairs = list(xy = cbind(x_lag, y_lead), yx = cbind(y_lag, x_lead)))
}


.nns_caus_index <- function(raw_xy, raw_yx, n_pairs) {
  k <- 1 / sqrt(n_pairs)
  c_xy <- (raw_xy + k) / (raw_xy + raw_yx + 2 * k)
  c("C(x--->y)" = c_xy, "C(y--->x)" = 1 - c_xy)
}


# tau = "ts": one lag for both directions, so the two raw scores are
# comparable: the mean of the first NNS.seas period <= sqrt(n) of each series.
.nns_caus_ts_lag <- function(x, y) {
  l <- length(x)
  x_tau <- NNS.seas(x, plot = FALSE)$periods
  y_tau <- NNS.seas(y, plot = FALSE)$periods
  x_tau <- x_tau[x_tau <= (l)^(1/2)][1]
  y_tau <- y_tau[y_tau <= (l)^(1/2)][1]
  if(is.na(x_tau) || is.na(y_tau))
    stop(paste0("tau = \"ts\": NNS.seas found no seasonal period <= sqrt(length(x)) for ",
                if(is.na(x_tau)) "x" else "y", "; supply a numeric tau instead."))
  min(as.integer(round(mean(c(x_tau, y_tau)))), l - 3L)
}


# Circular shift offsets at least max(tau + 1, n / 10) steps from alignment:
# the series keeps its own dynamics and only its alignment with the other
# series is broken. The permutation null for time series in NNS.caus.
.nns_caus_shift_offsets <- function(n, tau = 0L) {
  m <- max(as.integer(tau) + 1L, as.integer(ceiling(n / 10)))
  if (n - 2L * m >= 0L) m:(n - m) else seq_len(n - 1L)
}


# Distinct offsets only: a repeated one adds no information but shrinks the p-value.
.nns_caus_shift_draw <- function(n, tau, nperm) {
  K <- .nns_caus_shift_offsets(n, tau)
  if (length(K) <= nperm) K else K[sample.int(length(K), nperm)]
}


.nns_caus_shift <- function(v, k) {
  n <- length(v)
  c(v[(k + 1L):n], v[seq_len(k)])
}


# Raw directional causation of y on x: C~_{y -> x}. Uni.caus(response, X_j)
# is the causal weight of X_j on the response in NNS.reg and NNS.stack
# dimension reduction. Those models fit the aligned columns, so the weight is
# scored on aligned pairs for tau = "cs" and tau = "ts" alike (#96): a lagged
# score would describe pairs the model never fits.
Uni.caus <- function(x, y, tau, plot = FALSE){
  n <- min(length(x), length(y))
  r <- .nns_caus_raw(as.numeric(y[seq_len(n)]), as.numeric(x[seq_len(n)]), 0L)
  if (plot) .nns_caus_plot(r, 0L, .nns_caus_index(r$raw_xy, r$raw_yx, r$n_pairs))
  r$raw_xy
}


.nns_caus_plot <- function(r, tau, idx) {
  original.par <- par(no.readonly = TRUE)
  on.exit(par(original.par), add = TRUE)
  par(mfrow = c(1, 2))
  lab <- if (tau > 0) paste0("t-", tau) else "t"
  plot(r$pairs$xy, pch = 16, cex = .5, col = "steelblue",
       xlab = paste0("X(", lab, ") copula scale"), ylab = "Y(t) copula scale",
       main = sprintf("x ---> y: C = %.3f", idx[[1]]))
  plot(r$pairs$yx, pch = 16, cex = .5, col = "red",
       xlab = paste0("Y(", lab, ") copula scale"), ylab = "X(t) copula scale",
       main = sprintf("y ---> x: C = %.3f", idx[[2]]))
}


# Internal helper: signed causation x -> y from cp, in [-1, 1].
signed_from_cp <- function(cp){
  as.numeric(cp[["C(x--->y)"]] - cp[["C(y--->x)"]])
}


# Core causation computation (no permutation logic).
NNS.caus_core <- function(x, y = NULL,
                          factor.2.dummy = FALSE,
                          tau = 0,
                          plot = FALSE,
                          p.value = FALSE,
                          nperm = 100L,
                          permute = c("auto", "y", "x", "both", "shift"),
                          seed = NULL,
                          conf.int = 0.95){
  if(!is.null(y))  if(anyNA(x) || anyNA(y)) stop("You have some missing values, please address.")
  if(is.null(y))  if(anyNA(x)) stop("You have some missing values, please address.")

  orig.tau <- tau
  orig.plot <- plot

  if(is.data.frame(x) && ncol(x) == 1) x <- x[[1L]]
  if(any(class(x)%in%c("tbl","data.table"))) x <- as.data.frame(x)
  if(is.data.frame(y) && ncol(y) == 1) y <- y[[1L]]

  if(factor.2.dummy){
    if(!is.null(dim(x))){
      if(!is.numeric(x)){
        x <- do.call(cbind, lapply(x, factor_2_dummy_FR))
      } else {
        x <- apply(x, 2, as.double)
      }
      if(is.list(x)){
        x <- do.call(cbind, x)
        x <- apply(x, 2, as.double)
      }
    } else {
      x <- factor_2_dummy(x)
      if(is.null(dim(x))){
        x <- as.double(x)
      } else {
        x <- apply(x, 2, as.double)
      }
    }
  }

  if(is.null(y)) {
    return(NNS.caus.matrix(x, tau = orig.tau, factor.2.dummy = factor.2.dummy, plot = orig.plot,
                           p.value = p.value, nperm = nperm, conf.int = conf.int, seed = seed,
                           permute = match.arg(permute)))
  }

  if(is.factor(y)) y <- as.numeric(y)

  # ---- input validation ---------------------------------------------------
  if(NCOL(x) > 1 || NCOL(y) > 1)
    stop("x and y must each be a single variable (factor.2.dummy = TRUE turns a factor with k levels into k - 1 columns).")
  if(length(x) != length(y)) stop("x and y must have the same length.")
  if(any(is.infinite(x)) || any(is.infinite(y))) stop("x and y must be finite.")
  if(length(unique(x)) < 2 || length(unique(y)) < 2)
    stop("x and y must not be constant: causation is undefined for a variable with no variation.")
  if(!(identical(tau, "cs") || identical(tau, "ts"))){
    if(!is.numeric(tau) || length(tau) != 1 || !is.finite(tau) || tau < 0 ||
       tau != round(tau) || tau >= length(x) - 2)
      stop("tau must be \"cs\", \"ts\" or a single non-negative integer smaller than length(x) - 2.")
  }

  x <- as.numeric(x); y <- as.numeric(y)

  if(identical(tau, "cs")) {
    tau_use <- 0L
  } else if(identical(tau, "ts")) {
    tau_use <- .nns_caus_ts_lag(x, y)
  } else {
    tau_use <- as.integer(tau)
  }

  r <- .nns_caus_raw(x, y, tau_use)
  cp <- .nns_caus_index(r$raw_xy, r$raw_yx, r$n_pairs)
  if(plot) .nns_caus_plot(r, tau_use, cp)
  cp
}
