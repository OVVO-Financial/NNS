# Pairwise causation matrix: entry [i, j] is the directional causation index
# C(column i ---> column j) in [0, 1], so [i, j] + [j, i] = 1 off the diagonal
# and 0.5 is neutral. The diagonal is 0 (no self-causation).
# With p.value = TRUE, the null of each index is built by permutation with
# each pair's lag fixed at its observed value: for cross-sections the columns
# are shuffled independently; for time series (tau = "ts" or a positive lag)
# each pair keeps its first column and circularly shifts the second, so both
# series keep their own dynamics and the pair is always at least
# max(tau + 1, n / 10) steps from alignment. lower_CI /
# upper_CI hold the null's partial-moment interval and p.value_matrix the
# two-sided permutation p-value of C(i ---> j) - C(j ---> i).
NNS.caus.matrix <- function(x, tau = 0, factor.2.dummy = FALSE, plot = FALSE, p.value = FALSE, nperm = 100, conf.int = 0.95, seed = NULL){
  if(is.null(ncol(x))){
    stop("supply both 'x' and 'y' or a matrix-like 'x'")
  }
  n <- ncol(x)
  nm <- colnames(x)
  causes <- matrix(0, n, n, dimnames = list(nm, nm))
  pairs <- utils::combn(n, 2)

  pair_index <- function(xm, i, j, plot = FALSE, tau_ij = tau) {
    NNS.caus(xm[, i], xm[, j], plot = plot, tau = tau_ij, factor.2.dummy = factor.2.dummy)
  }

  for(k in seq_len(ncol(pairs))){
    i <- pairs[1, k]; j <- pairs[2, k]
    cp <- pair_index(x, i, j, plot = plot)
    causes[i, j] <- cp[["C(x--->y)"]]
    causes[j, i] <- cp[["C(y--->x)"]]
  }
  diag(causes) <- 0
  causes[is.na(causes)] <- 0

  if(!p.value) return(causes)

  if(!is.null(seed)) set.seed(seed)
  lower_CI <- matrix(0, n, n, dimnames = list(nm, nm))
  upper_CI <- matrix(0, n, n, dimnames = list(nm, nm))
  p.value_matrix <- matrix(0, n, n, dimnames = list(nm, nm))

  # Each pair's lag is resolved once on the observed columns (#96), so the
  # permuted columns are compared at the same lag.
  pair_tau <- vapply(seq_len(ncol(pairs)), function(k) {
    if (identical(tau, "ts"))
      .nns_caus_ts_lag(as.numeric(x[, pairs[1, k]]), as.numeric(x[, pairs[2, k]]))
    else if (identical(tau, "cs")) 0L else as.integer(tau)
  }, integer(1L))
  time_series <- any(pair_tau > 0L)
  # Permute column by column so a data frame keeps its column types: apply()
  # would turn a factor column into character (#96).
  x_cols <- if (is.data.frame(x)) x else as.data.frame(x)
  null_mat <- array(NA_real_, dim = c(nperm, ncol(pairs)))
  for(b in seq_len(nperm)){
    if (time_series) {
      # Shift one series per pair: shifting every column independently would
      # leave pairs whose two shifts nearly coincide still aligned.
      for(k in seq_len(ncol(pairs))){
        i <- pairs[1, k]; j <- pairs[2, k]
        x_pair <- x_cols[, c(i, j)]
        x_pair[[2]] <- .nns_caus_shift(x_cols[[j]], pair_tau[k])
        null_mat[b, k] <- pair_index(x_pair, 1L, 2L, tau_ij = pair_tau[k])[["C(x--->y)"]]
      }
    } else {
      x_perm <- x_cols
      for(j in seq_len(n)) x_perm[[j]] <- sample(x_cols[[j]])
      for(k in seq_len(ncol(pairs))){
        null_mat[b, k] <- pair_index(x_perm, pairs[1, k], pairs[2, k], tau_ij = pair_tau[k])[["C(x--->y)"]]
      }
    }
  }

  p <- (1 - conf.int)/2
  for(k in seq_len(ncol(pairs))){
    i <- pairs[1, k]; j <- pairs[2, k]
    null_ij <- null_mat[, k]
    lower_CI[i, j] <- LPM.VaR(p, 0, null_ij)
    upper_CI[i, j] <- UPM.VaR(p, 0, null_ij)
    lower_CI[j, i] <- LPM.VaR(p, 0, 1 - null_ij)
    upper_CI[j, i] <- UPM.VaR(p, 0, 1 - null_ij)
    obs_signed <- causes[i, j] - causes[j, i]
    null_signed <- 2 * null_ij - 1
    pval <- (1 + sum(abs(null_signed) >= abs(obs_signed))) / (1 + nperm)
    p.value_matrix[i, j] <- pval
    p.value_matrix[j, i] <- pval
  }
  diag(lower_CI) <- diag(upper_CI) <- diag(p.value_matrix) <- 0
  lower_CI[is.na(lower_CI)] <- 0
  upper_CI[is.na(upper_CI)] <- 0
  p.value_matrix[is.na(p.value_matrix)] <- 0

  return(.NNS.out(list(
    causality = causes,
    lower_CI = lower_CI,
    upper_CI = upper_CI,
    p.value_matrix = p.value_matrix
  )))
}
