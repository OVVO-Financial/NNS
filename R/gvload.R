# Import calls and globalvariable calls

#' @importFrom grDevices adjustcolor rainbow rgb
#' @importFrom graphics abline boxplot legend lines par plot points segments text matplot title axis mtext
#' @importFrom stats coef cor lm na.omit sd median complete.cases resid uniroot aggregate density
#' @import rgl
#' @import data.table
#' @import stringr
#' @importFrom utils globalVariables


.onLoad <- function(libname = find.package("NNS"), pkgname = "NNS"){

  # CRAN Note avoidance

  utils::globalVariables(
    c("quadrant","prior.quadrant",".","min_x_seg","max_x_seg","min_y_seg","max_y_seg",
      "mean_y_seg","mean_x_seg","sub.clpm",'sub.cupm','sub.dlpm','sub.dupm','weight',
      "Coefficient","X.Lower.Range","X.Upper.Range","y.hat",
      "NNS.ID","max.x1","max.x2","min.x1","min.x2","counts",
      "Period","Coefficient.of.Variance","Variable.Coefficient.of.Variance"
      ))
  invisible()
}
