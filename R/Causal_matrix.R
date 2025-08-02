# Efficient antisymmetric causal matrix using pairwise signed net causation.
# Optionally returns permutation-based lower and upper CI matrices if p.value = TRUE.
NNS.caus.matrix <- function(x, tau = 0, factor.2.dummy = FALSE, plot = FALSE, p.value = FALSE, nperm = 100, conf.int = 0.95, seed = NULL){
  if(is.null(ncol(x))){
    stop("supply both 'x' and 'y' or a matrix-like 'x'")
  }
  n <- ncol(x)
  causes <- matrix(0, n, n, dimnames = list(colnames(x), colnames(x)))
  pairs <- utils::combn(n, 2)
  
  for(k in seq_len(ncol(pairs))){
    i <- pairs[1, k]
    j <- pairs[2, k]
    cp <- NNS.caus(x[, i], x[, j], plot = plot, tau = tau, factor.2.dummy = factor.2.dummy)
    val_ij <- if(names(cp)[3] == "C(x--->y)"){
      as.numeric(cp[3])
    } else if(names(cp)[3] == "C(y--->x)"){
      -as.numeric(cp[3])
    } else {
      as.numeric(cp[3])
    }
    causes[i, j] <- -val_ij  
    causes[j, i] <- val_ij
  }
  diag(causes) <- 0
  causes[is.na(causes)] <- 0
  
  if(!p.value) return(causes)
  
  if(!is.null(seed)) set.seed(seed)
  lower_CI <- matrix(0, n, n, dimnames = list(colnames(x), colnames(x)))
  upper_CI <- matrix(0, n, n, dimnames = list(colnames(x), colnames(x)))
  
  null_mat <- array(NA, dim = c(nperm, ncol(pairs)))
  
  for(b in seq_len(nperm)){
    x_perm <- apply(x, 2, sample)
    for(k in seq_len(ncol(pairs))){
      i <- pairs[1, k]; j <- pairs[2, k]
      cp_perm <- NNS.caus(x_perm[, i], x_perm[, j], plot = plot, tau = tau, factor.2.dummy = factor.2.dummy)
      third_name <- names(cp_perm)[3]
      net_val <- as.numeric(cp_perm[3])  # already normalized signed log-ratio
      val_ij <- if(third_name == "C(x--->y)") net_val else if(third_name == "C(y--->x)") -net_val else net_val
      null_mat[b, k] <- val_ij
    }
  }
  
  for(k in seq_len(ncol(pairs))){
    i <- pairs[1, k]; j <- pairs[2, k]
    null_vals <- null_mat[, k]  # normalized
    p <- (1 - conf.int)/2
    lower <- LPM.VaR(p, 0, null_vals)
    upper <- UPM.VaR(p, 0, null_vals)
    lower_CI[j, i] <- lower
    upper_CI[j, i] <- upper
    lower_CI[i, j] <- lower
    upper_CI[i, j] <- upper
  }
  
  diag(lower_CI) <- diag(upper_CI) <- 0
  lower_CI[is.na(lower_CI)] <- 0
  upper_CI[is.na(upper_CI)] <- 0
  
  p.value_matrix <- matrix(0, n, n, dimnames = list(colnames(x), colnames(x)))
  for(k in seq_len(ncol(pairs))){
    i <- pairs[1, k]; j <- pairs[2, k]
    null_vals_trans <- null_mat[, k]  # normalized
    obs_ij <- causes[i, j]
    pval <- (1 + sum(abs(null_vals_trans) >= abs(obs_ij))) / (1 + nperm)
    p.value_matrix[i, j] <- pval
    p.value_matrix[j, i] <- pval
  }
  diag(p.value_matrix) <- 0
  p.value_matrix[is.na(p.value_matrix)] <- 0
  
  return(list(
    causality = causes,
    lower_CI = lower_CI,
    upper_CI = upper_CI,
    p.value_matrix = p.value_matrix
  ))
}