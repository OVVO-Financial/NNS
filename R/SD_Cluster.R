#' NNS SD-based Clustering
#'
#' Clusters a set of variables by iteratively extracting Stochastic Dominance (SD)-efficient sets,
#' subject to a minimum cluster size.
#'
#' @param data A numeric matrix or data frame of variables to be clustered.
#' @param degree Numeric options: (1, 2, 3). Degree of stochastic dominance test.
#' @param type Character, either  \code{"discrete"} (default) or \code{"continuous"}; specifies the type of CDF.
#' @param min_cluster Integer. The minimum number of elements required for a valid cluster.
#' @param dendrogram Logical; \code{FALSE} (default). If \code{TRUE}, a dendrogram is produced based on a simple "distance" measure between clusters.
#'
#' @return
#' A list with the following components:
#' \itemize{
#'   \item \code{Clusters}: A named list of cluster memberships where each element is the set of variable names belonging to that cluster.
#'   \item \code{Dendrogram} (optional): If \code{dendrogram = TRUE}, an \code{hclust} object is also returned.
#' }
#'
#' @details
#' The function applies \code{\link{NNS.SD.efficient.set}} iteratively, peeling off the SD-efficient set at each step
#' if it meets or exceeds \code{min_cluster} in size, until no more subsets can be extracted or all variables are exhausted.
#' Variables in each SD-efficient set form a cluster, with any remaining variables aggregated into the final cluster if it meets
#' the \code{min_cluster} threshold.
#'
#' @author Fred Viole, OVVO Financial Systems
#'
#' @references Viole, F. and Nawrocki, D. (2016) "LPM Density Functions for the Computation of the SD Efficient Set." Journal of Mathematical Finance, 6, 105-126.  \doi{10.4236/jmf.2016.61012}.
#'
#' Viole, F. (2017) "A Note on Stochastic Dominance." \doi{10.2139/ssrn.3002675}
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' x <- rnorm(100)
#' y <- rnorm(100)
#' z <- rnorm(100)
#' A <- cbind(x, y, z)
#'
#' # Perform SD-based clustering (degree 1), requiring at least 2 elements per cluster
#' results <- NNS.SD.cluster(data = A, degree = 1, min_cluster = 2)
#' print(results$Clusters)
#'
#' # Produce a dendrogram as well
#' results_with_dendro <- NNS.SD.cluster(data = A, degree = 1, min_cluster = 2, dendrogram = TRUE)
#' }
#'
#' @export
 

NNS.SD.cluster <- function(data, degree = 1, type = "discrete", min_cluster = 1, dendrogram = FALSE) {
  clusters <- list()
  iteration <- 1
  n <- ncol(data)
  
  if(is.null(colnames(data))) colnames(data) <- paste0("X_",1:ncol(data))
  original_names <- colnames(data)
  
  # Ensure the input data is a matrix
  remaining_data <- as.matrix(data)
  
  
  # Continue clustering until the number of remaining columns is less than or equal to min_cluster
  while (ncol(remaining_data) > min_cluster) {
    # Use the original NNS.SD.efficient.set call as provided
    SD_set <- NNS.SD.efficient.set(remaining_data, degree = degree, type = type, status = FALSE)
    
    if (length(SD_set) == 0) {
      break
    }
    
    # Store the SD-efficient set as a cluster
    clusters[[paste0("Cluster_", iteration)]] <- SD_set
    
    # Remove the identified SD set from remaining_data
    remaining_data <- remaining_data[, !(colnames(remaining_data) %in% SD_set), drop = FALSE]
    
    # Ensure remaining_data remains a matrix
    remaining_data <- as.matrix(remaining_data)
    
    iteration <- iteration + 1
    
    # If the number of remaining columns is now less than or equal to min_cluster, add them as the final cluster
    if (ncol(remaining_data) <= min_cluster) {
      clusters[[paste0("Cluster_", iteration)]] <- colnames(remaining_data)
      break
    }
  }
  
  # If there are still variables left (and not already added), add them as the final cluster
  if (ncol(remaining_data) > min_cluster && !paste0("Cluster_", iteration) %in% names(clusters)) {
    clusters[[paste0("Cluster_", iteration)]] <- colnames(remaining_data)
  }
  
  # Check if the final cluster has fewer elements than min_cluster; if so, merge it with the previous cluster (if one exists)
  final_cluster_name <- paste0("Cluster_", length(clusters))
  if (length(clusters[[final_cluster_name]]) < min_cluster && length(clusters) > 1) {
    previous_cluster_name <- paste0("Cluster_", length(clusters) - 1)
    clusters[[previous_cluster_name]] <- c(clusters[[previous_cluster_name]], clusters[[final_cluster_name]])
    clusters[[final_cluster_name]] <- NULL
  }
  
  # Flatten the clusters into a single vector and generate cluster labels
  all_vars <- unlist(clusters)
  
  
  
  cluster_labels <- unlist(lapply(seq_along(clusters), function(i) rep(i, length(clusters[[i]]))))
  
  
  if(dendrogram){
    # Ensure there are at least two variables for hierarchical clustering
    if (length(all_vars) < 2) {
      warning("Not enough variables for hierarchical clustering. Returning clusters only.")
      return(.NNS.out(list("Clusters" = clusters, "Order" = NULL)))
    }
    
    # Use the extraction order inherent in all_vars as a tie-breaker.
    extraction_order <- seq_along(all_vars)

    if(length(clusters)==1) epsilon <- 0 else epsilon <- 1e-3  # small tie-breaker weight
    dist_matrix <- as.dist(
      outer(cluster_labels, cluster_labels, function(a, b) n * abs(a - b)) +
        epsilon * outer(extraction_order, extraction_order, function(i, j) abs(i - j))
    )
    attr(dist_matrix, "Labels") <- all_vars

    # Perform hierarchical clustering
    hc <- hclust(dist_matrix, method = "complete")
    
    plot(hc,
         main = paste0("Hierarchical Clustering of Stochastic Dominance Sets \nSD Degree: ", degree),
         xlab = "Variables",
         ylab = "SD Distance",
         sub = ""
    )

    hc$order <- match(hc$labels, original_names)
    
    return(.NNS.out(list("Clusters" = clusters, "Dendrogram" = hc)))
  } else return(.NNS.out(list("Clusters" = clusters)))
}

