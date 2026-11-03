# Recreate Chapter 18's mean-reducer partition figures from the repository root:
# Rscript book/scripts/render_foundation_figures.R
library(NNS)

x <- seq(-5, 5, 0.05)
y <- x^3

render_partition <- function(file, type = NULL) {
  png(file, width = 1600, height = 1450, res = 180, type = "cairo")
  on.exit(dev.off())
  par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3.5, 1.2), cex = 0.85)
  for (i in 1:4) {
    NNS.part(x, y, order = i, type = type, Voronoi = TRUE,
             obs.req = 0, noise.reduction = "mean")
  }
}

dir.create("book/images", showWarnings = FALSE)
render_partition("book/images/ch18_part_1.png")
render_partition("book/images/ch18_part_2.png", type = "XONLY")
