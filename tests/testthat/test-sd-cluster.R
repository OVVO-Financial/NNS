test_that("the returned dendrogram keeps hclust's leaf order", {
  set.seed(123)
  A <- matrix(rnorm(200), 10, 20)
  colnames(A) <- paste0("V", 1:20)
  pdf(NULL); on.exit(dev.off())
  res <- NNS.SD.cluster(A, degree = 1, min_cluster = 2, dendrogram = TRUE)

  cl <- res$Clusters
  expect_gt(length(cl), 1)
  v <- unlist(cl, use.names = FALSE)
  lab <- rep(seq_along(cl), lengths(cl))
  k <- seq_along(v)
  d <- as.dist(outer(lab, lab, function(a, b) 20 * abs(a - b)) + 1e-3 * outer(k, k, function(i, j) abs(i - j)))
  attr(d, "Labels") <- v
  hc <- hclust(d, method = "complete")

  expect_identical(res$Dendrogram$order, hc$order)
  expect_identical(order.dendrogram(as.dendrogram(res$Dendrogram)), res$Dendrogram$order)
  expect_identical(unname(labels(as.dendrogram(res$Dendrogram))), v[hc$order])
})

test_that("variables left when peeling stops at once form the final cluster", {
  A <- cbind(A = 1:5, B = 2:6)
  expect_identical(NNS.SD.cluster(A, min_cluster = 2)$Clusters, list(Cluster_1 = c("A", "B")))
  expect_identical(NNS.SD.cluster(A, min_cluster = 3)$Clusters, list(Cluster_1 = c("A", "B")))
  expect_identical(NNS.SD.cluster(A[, "A", drop = FALSE])$Clusters, list(Cluster_1 = "A"))

  pdf(NULL); on.exit(dev.off())
  res <- NNS.SD.cluster(A, min_cluster = 2, dendrogram = TRUE)
  expect_identical(res$Clusters, list(Cluster_1 = c("A", "B")))
  expect_identical(unname(res$Dendrogram$labels), c("A", "B"))
  res <- NNS.SD.cluster(cbind(A = 1:5, B = 11:15), dendrogram = TRUE)
  expect_identical(res$Clusters, list(Cluster_1 = "B", Cluster_2 = "A"))
  expect_warning(res <- NNS.SD.cluster(A[, "A", drop = FALSE], dendrogram = TRUE), "Not enough variables")
  expect_identical(res$Clusters, list(Cluster_1 = "A"))
})
