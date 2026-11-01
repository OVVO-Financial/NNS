# Directional causation (book ch14), redesigned for #68.
#
# For a lag tau >= 0 the two directions are evaluated on explicitly lagged
# pairs:
#   x -> y on (X_{t-tau}, Y_t),   y -> x on (Y_{t-tau}, X_t).
# Each pair is mapped to the copula scale with per-column ranks (the degree-0
# lower partial moment, the empirical CDF), so the statistic is invariant to
# shifts, rescaling and any monotone transform of either variable. The lag
# supplies the temporal asymmetry; at tau = 0 the asymmetry comes from the
# directional measures themselves.
#
# Each direction gets a raw score in [0, 1],
#   C~_{x -> y} = (P_{x -> y} + max(rho_{x -> y} - rho_{y -> x}, 0)) / 2,
# where
#   P_{x -> y}   share of Y_t's degree-1 deviation (L1 + U1 about its mean)
#                resolved by conditioning on the position of X_{t-tau}
#                (equal-frequency cells of the copula-scale cause), and
#   rho_{x -> y} NNS.dep(X_{t-tau}, Y_t, asym = TRUE)$Dependence, the
#                directional dependence of Y_t on X_{t-tau}, built from the
#                lagged co-partial moments of the pair.
# The index is
#   C_{x -> y} = (C~_{x -> y} + k) / (C~_{x -> y} + C~_{y -> x} + 2k),
#   C_{y -> x} = 1 - C_{x -> y},
# with a neutral pseudo-count k = 1 / sqrt(n - tau), the order of the raw
# scores' sampling variation, so weak or absent evidence stays near 0.5.


# Copula scale: the degree-0 LPM of each observation against its own column.
.nns_caus_u <- function(v) LPM.ratio(0, v, v)


# Share of the effect's degree-1 deviation resolved by the cause's position.
# Both inputs are on the copula scale and have the same length.
.nns_caus_P <- function(cause_u, effect_u) {
  n <- length(cause_u)
  if (n < 3L) return(0)
  k <- max(2L, as.integer(floor(sqrt(n))))
  cell <- pmax(1L, as.integer(ceiling(cause_u * k)))
  # Degree-1 deviation about the mean, L1 + U1, with a small-sample factor so
  # independence scores near 0 instead of a positive bias of order 1 / m.
  dev <- function(v) {
    m <- length(v)
    if (m < 2L) return(NA_real_)
    mu <- mean(v)
    (LPM(1, mu, v) + UPM(1, mu, v)) * sqrt(m / (m - 1))
  }
  total <- dev(effect_u)
  if (!is.finite(total) || total <= 0) return(0)
  parts <- split(effect_u, cell)
  w <- vapply(parts, length, integer(1L))
  d <- vapply(parts, dev, numeric(1L))
  ok <- is.finite(d)
  if (!any(ok)) return(0)
  within <- sum(w[ok] * d[ok]) / sum(w[ok])
  min(1, max(0, 1 - within / total))
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
  dep <- function(a, b) {
    z <- tryCatch(NNS.dep(a, b, asym = TRUE)$Dependence, error = function(e) 0)
    if (is.finite(z)) z else 0
  }
  rho_xy <- dep(x_lag, y_lead)
  rho_yx <- dep(y_lag, x_lead)

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


.nns_caus_tau <- function(tau) {
  if (identical(tau, "cs")) return(0L)
  if (identical(tau, "ts")) return(3L)
  as.integer(tau)
}


# Raw directional causation of y on x: C~_{y -> x}. Uni.caus(response, X_j)
# is the causal weight of X_j on the response in NNS.reg and NNS.stack
# dimension reduction. tau = "ts" uses a lag of 3 there.
Uni.caus <- function(x, y, tau, plot = FALSE){
  tau <- .nns_caus_tau(tau)
  n <- min(length(x), length(y))
  if (tau >= n - 2L) tau <- max(0L, n - 3L)
  r <- .nns_caus_raw(as.numeric(y[seq_len(n)]), as.numeric(x[seq_len(n)]), tau)
  if (plot) .nns_caus_plot(r, tau, .nns_caus_index(r$raw_xy, r$raw_yx, r$n_pairs))
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
                          permute = c("y","x","both"),
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
                           p.value = p.value, nperm = nperm, conf.int = conf.int, seed = seed))
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
    # One lag for both directions, so the two raw scores are comparable: the
    # mean of the first NNS.seas period <= sqrt(n) of each series.
    l <- length(x)
    x_tau <- NNS.seas(x, plot = FALSE)$periods
    y_tau <- NNS.seas(y, plot = FALSE)$periods
    x_tau <- x_tau[x_tau <= (l)^(1/2)][1]
    y_tau <- y_tau[y_tau <= (l)^(1/2)][1]
    if(is.na(x_tau) || is.na(y_tau))
      stop(paste0("tau = \"ts\": NNS.seas found no seasonal period <= sqrt(length(x)) for ",
                  if(is.na(x_tau)) "x" else "y", "; supply a numeric tau instead."))
    tau_use <- min(as.integer(round(mean(c(x_tau, y_tau)))), l - 3L)
  } else {
    tau_use <- as.integer(tau)
  }

  r <- .nns_caus_raw(x, y, tau_use)
  cp <- .nns_caus_index(r$raw_xy, r$raw_yx, r$n_pairs)
  if(plot) .nns_caus_plot(r, tau_use, cp)
  cp
}
