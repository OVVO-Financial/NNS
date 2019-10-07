# Import calls and globalvariable calls

#' @importFrom grDevices adjustcolor rainbow rgb
#' @importFrom graphics abline boxplot legend lines par plot points segments text matplot title axis mtext barplot hist strwidth
#' @importFrom stats coef cor lm na.omit sd median complete.cases resid uniroot aggregate density hat qnorm model.matrix fivenum
#' @importFrom utils globalVariables head tail combn flush.console
#' @importFrom data.table data.table %chin% .I .N .SD := as.data.table fwrite is.data.table rbindlist set setcolorder setnames setorderv as.IDate as.ITime
#' @import data.table
#' @import doParallel
#' @import rgl
#' @import stringr


.onLoad <- function(libname = find.package("NNS"), pkgname = "NNS"){

  # CRAN Note avoidance

  utils::globalVariables(
    c("quadrant","quadrant.new","prior.quadrant",".","tmp.x","tmp.y","min_x_seg","max_x_seg","min_y_seg","max_y_seg",
      "mean_y_seg","mean_x_seg","sub.clpm",'sub.cupm','sub.dlpm','sub.dupm','weight','mean.x','mean.y',
      "Coefficient","X.Lower.Range","X.Upper.Range","y.hat","interval",
      "NNS.ID","max.x1","max.x2","min.x1","min.x2","counts",'old.counts',
      "Period","Coefficient.of.Variance","Variable.Coefficient.of.Variance", "Sum", "j","lpm","upm",
      "i.x","i.y","q_new","x.x","x.y","standard.errors",
      "detectCores","makeCluster","registerDoSEQ","clusterExport",
      "%dopar%","foreach","stopCluster",
      "%do%", "k"
    ))

  requireNamespace("data.table")
  requireNamespace("rgl")
  requireNamespace("doParallel")
  requireNamespace("stringr")

  .datatable.aware = TRUE

  invisible()


}
