NNS.M.reg <- function (X_n, Y, factor.2.dummy = TRUE, order = NULL, n.best = NULL, type = NULL, point.est = NULL, point.only = FALSE,
                       plot = FALSE, residual.plot = TRUE, location = NULL, noise.reduction = 'off', dist = "L2",
                       return.values = FALSE, plot.regions = FALSE, ncores = NULL, confidence.interval = NULL){
  
  dist <- tolower(dist)
  
  ### For Multiple regressions
  ###  Turn each column into numeric values
  original.IVs <- X_n
  original.DV <- Y
  n <- ncol(original.IVs)
  
  feature.names <- colnames(original.IVs)
  if(is.null(feature.names)){
    feature.names <- paste0("x", 1:n)
  } else {
    missing.feature.names <- is.na(feature.names) | feature.names == ""
    feature.names[missing.feature.names] <- paste0("x", which(missing.feature.names))
    feature.names <- make.unique(feature.names, sep = ".")
  }
  colnames(original.IVs) <- feature.names
  
  if(is.null(ncol(X_n))) X_n <- t(t(X_n))
  
  if(is.null(names(Y))){
    y.label <- "Y"
  } else {
    y.label <- names(Y)
  }
  
  np <- nrow(point.est)
  
  if(is.null(np) && !is.null(point.est)) point.est <- t(point.est)
  
  if(!is.null(point.est)){
    if(ncol(point.est) != n){
      stop("Please ensure 'point.est' is of compatible dimensions to 'x'")
    }
    colnames(point.est) <- feature.names
  }
  
  original.matrix <- cbind.data.frame(original.DV, original.IVs)
  norm.matrix <- apply(original.matrix, 2, function(z) NNS.rescale(z, 0, 1))
  
  minimums <- apply(original.IVs, 2, min)
  maximums <- apply(original.IVs, 2, max)
  
  ###  Regression Point Matrix
  if(is.numeric(order) || is.null(order)){
    reg.points <- lapply(1:ncol(original.IVs), function(b) NNS.reg(original.IVs[, b], original.DV, factor.2.dummy = factor.2.dummy, order = order, type = type, noise.reduction = noise.reduction, plot = FALSE, multivariate.call = TRUE, ncores = 1)$x)
    
    if(length(unique(sapply(reg.points, length))) != 1){
      reg.points.matrix <- do.call(cbind, lapply(reg.points, `length<-`, max(lengths(reg.points))))
    } else {
      reg.points.matrix <- do.call(cbind, reg.points)
    }
  } else {
    reg.points.matrix <- original.IVs
  }
  
  ### If regression points are error (not likely)...
  if(length(reg.points.matrix[ , 1]) == 0  || is.null(reg.points.matrix)){
    stn <- .95
    for(i in 1 : n){
      part.map <- NNS.part(original.IVs[ , i], original.DV, order = order, type = type, noise.reduction = noise.reduction, obs.req = 0)
      dep <- NNS.dep(original.IVs[ , i], original.DV)$Dependence
      char_length_order <- dep * max(nchar(part.map$df$quadrant))
      if(dep > stn){
        reg.points[[i]] <- NNS.part(original.IVs[ , i], original.DV, order = ifelse(char_length_order%%1 < .5, floor(char_length_order), ceiling(char_length_order)), type = type, noise.reduction = 'off', obs.req = 0)$regression.points$x
      } else {
        reg.points[[i]] <- NNS.part(original.IVs[ , i], original.DV, order = ifelse(char_length_order%%1 < .5, floor(char_length_order), ceiling(char_length_order)), noise.reduction = noise.reduction, type = "XONLY", obs.req = 1)$regression.points$x
      }
    }
    reg.points.matrix <- do.call('cbind', lapply(reg.points, `length<-`, max(lengths(reg.points))))
  }
  
  colnames(reg.points.matrix) <- feature.names
  
  if(is.numeric(order) || is.null(order)) reg.points.matrix <- unique(reg.points.matrix)
  
  if(!is.null(order) && order=="max" && is.null(n.best)) n.best <- 1
  
  ### Determine core configuration for native C++ multi-threading
  if(is.null(ncores)){
    num_cores <- as.integer(max(1L, parallel::detectCores(), na.rm = TRUE)) - 1
    if(num_cores < 1L) num_cores <- 1L
  } else {
    num_cores <- as.integer(ncores)
  }
  
  NNS.ID <- lapply(1:n, function(j) findInterval(original.IVs[ , j], vec = na.omit(sort(reg.points.matrix[ , j])), left.open = FALSE))
  
  NNS.ID <- do.call(cbind, NNS.ID)
  
  ### Create unique identifier of each observation's interval
  NNS.ID <- gsub(do.call(paste, as.data.frame(NNS.ID)), pattern = " ", replacement = ".")
  
  ### Grouped regression-point matrix and fitted values by NNS.ID, computed in
  ### C++ (NNS_mreg_reduce_cpp): each group's central tendency (per
  ### noise.reduction) of the IV columns and the DV, mapped back to every
  ### observation, then the grouped residual-bias correction.
  order_is_numeric <- is.numeric(order) || is.null(order)
  reducer <- if(!order_is_numeric) 5L else switch(noise.reduction,
                                                  "mean"       = 1L,
                                                  "median"     = 2L,
                                                  "mode"       = 3L,
                                                  "mode_class" = 4L,
                                                  0L)   # "off" / default -> gravity

  red <- NNS_mreg_reduce_cpp(as.matrix(original.IVs), as.numeric(original.DV),
                             as.character(NNS.ID), as.integer(reducer), !is.null(type))

  ### REGRESSION.POINT.MATRIX: one row per unique NNS.ID (sorted-id order)
  REGRESSION.POINT.MATRIX <- as.data.frame(red$rpm, stringsAsFactors = FALSE)
  colnames(REGRESSION.POINT.MATRIX) <- c(feature.names, "y.hat")

  ### Fitted y.hat in original observation order (bias-corrected)
  y.hat <- red$yhat

  fitted.matrix <- as.data.frame(original.IVs, stringsAsFactors = FALSE)
  colnames(fitted.matrix) <- feature.names
  fitted.matrix$y         <- original.DV
  fitted.matrix$y.hat     <- red$yhat
  fitted.matrix$NNS.ID    <- NNS.ID
  fitted.matrix$residuals <- red$residuals
  
  if(is.null(n.best)){
    dependence <- NNS.copula(cbind(original.IVs, original.DV))
    n.best <- max(1, floor((1-dependence)*sqrt(n)))
  }
  
  ### Clamp n.best to available RPM rows.
  ### Oversized n.best means "use all available RPM".
  rpm_n <- nrow(REGRESSION.POINT.MATRIX)
  
  if (identical(n.best, "all") ||
      (is.numeric(n.best) && length(n.best) == 1L && is.infinite(n.best))) {
    n.best <- rpm_n
  } else {
    n.best <- suppressWarnings(as.integer(n.best[1L]))
    if (is.na(n.best)) n.best <- rpm_n
    n.best <- max(1L, min(n.best, rpm_n))
  }
  
  # OPTIMIZED: Bulk prediction calculation bypasses row-by-row mapping loops.
  # Use the single-k path kernel because only column n.best was consumed.
  if(n.best > 1 && !point.only){
    fitted.matrix$y.hat <- as.numeric(NNS.distance.path.single.bulk(
      rpm = REGRESSION.POINT.MATRIX,
      Xtest = original.IVs,
      k = n.best,
      class = type,
      ncores = num_cores
    ))
    
    y.hat <- fitted.matrix$y.hat
    if(!is.null(type)) y.hat <- ifelse(y.hat %% 1 < 0.5, floor(y.hat), ceiling(y.hat))
  }
  
  ### Point Estimates
  if (!is.null(point.est)) {
    # Calculate central points
    central.points <- apply(REGRESSION.POINT.MATRIX[, 1:n, drop = FALSE], 2, gravity)
    
    predict.fit <- numeric()
    outsiders <- point.est < minimums | point.est > maximums
    outsiders[is.na(outsiders)] <- 0
    
    # Single point estimation
    if (is.null(np)) {
      if (!any(outsiders)) {
        predict.fit <- NNS::NNS.distance(
          rpm = REGRESSION.POINT.MATRIX,
          dist.estimate = point.est,
          k = n.best,
          class = type
        )
      } else {
        boundary.points <- pmin(pmax(point.est, minimums), maximums)
        mid.points <- (boundary.points + central.points) / 2
        mid.points_2 <- (boundary.points + mid.points) / 2
        
        last.known.distances <- c(
          sqrt(sum((boundary.points - central.points) ^ 2)),
          sqrt(sum((boundary.points - mid.points) ^ 2)),
          sqrt(sum((boundary.points - mid.points_2) ^ 2))
        )
        
        boundary.estimates <- NNS::NNS.distance(
          rpm = REGRESSION.POINT.MATRIX,
          dist.estimate = boundary.points,
          k = n.best,
          class = type
        )
        
        gradients <- sapply(1:3, function(i) {
          compare.points <- list(central.points, mid.points, mid.points_2)[[i]]
          (boundary.estimates - NNS::NNS.distance(
            rpm = REGRESSION.POINT.MATRIX,
            dist.estimate = compare.points,
            k = n.best,
            class = type
          )) / last.known.distances[i]
        })
        
        last.known.gradient <- sum(gradients * c(3, 2, 1)) / 6
        last.distance <- sqrt(sum((point.est - boundary.points) ^ 2))
        
        predict.fit <- last.distance * last.known.gradient + boundary.estimates
      }
    }
    
    # Multiple point estimation
    if (!is.null(np)) {
      # OPTIMIZED: Replaced row-by-row distance operations with a single-k bulk call.
      DISTANCES <- as.numeric(NNS.distance.path.single.bulk(
        rpm = REGRESSION.POINT.MATRIX,
        Xtest = point.est,
        k = n.best,
        class = type,
        ncores = num_cores
      ))
      
      # OPTIMIZED: Fully vectorized matrix handling for out-of-bounds outliers
      if (any(rowSums(outsiders) > 0)) {
        outsider.indices <- which(rowSums(outsiders) > 0)
        outside.points_matrix <- as.matrix(point.est[outsider.indices, , drop = FALSE])
        
        boundary.points_matrix <- outside.points_matrix
        for (j in 1:ncol(boundary.points_matrix)) {
          boundary.points_matrix[, j] <- pmin(pmax(boundary.points_matrix[, j], minimums[j]), maximums[j])
        }
        
        # Add/subtract the central point by column using R's matrix recycling.
        # This is equivalent to column-wise central-point offsets while
        # avoiding helper dispatch in this high-level orchestration path.
        central.offset_matrix <- rep(central.points, each = nrow(boundary.points_matrix))
        mid.points_matrix <- (boundary.points_matrix + central.offset_matrix) / 2
        mid.points_2_matrix <- (boundary.points_matrix + mid.points_matrix) / 2
        
        last.known.distances_1 <- sqrt(rowSums((boundary.points_matrix - central.offset_matrix)^2))
        last.known.distances_2 <- sqrt(rowSums((boundary.points_matrix - mid.points_matrix)^2))
        last.known.distances_3 <- sqrt(rowSums((boundary.points_matrix - mid.points_2_matrix)^2))
        
        boundary.estimates <- as.numeric(NNS.distance.path.single.bulk(rpm = REGRESSION.POINT.MATRIX, Xtest = boundary.points_matrix, k = n.best, class = type, ncores = num_cores))
        mid.estimates <- as.numeric(NNS.distance.path.single.bulk(rpm = REGRESSION.POINT.MATRIX, Xtest = mid.points_matrix, k = n.best, class = type, ncores = num_cores))
        mid_2.estimates <- as.numeric(NNS.distance.path.single.bulk(rpm = REGRESSION.POINT.MATRIX, Xtest = mid.points_2_matrix, k = n.best, class = type, ncores = num_cores))
        
        central.estimate_single <- NNS.distance(rpm = REGRESSION.POINT.MATRIX, dist.estimate = central.points, k = n.best, class = type)[1]
        
        g1 <- (boundary.estimates - central.estimate_single) / pmax(last.known.distances_1, 1e-10)
        g2 <- (boundary.estimates - mid.estimates) / pmax(last.known.distances_2, 1e-10)
        g3 <- (boundary.estimates - mid_2.estimates) / pmax(last.known.distances_3, 1e-10)
        
        last.known.gradient <- (g1 * 3 + g2 * 2 + g3 * 1) / 6
        last.distance <- sqrt(rowSums((outside.points_matrix - boundary.points_matrix)^2))
        
        DISTANCES[outsider.indices] <- last.distance * last.known.gradient + boundary.estimates
      }
      predict.fit <- DISTANCES
    }
    
    if (point.only) {
      return(list(Point.est = predict.fit, RPM = REGRESSION.POINT.MATRIX))
    }
  } else {
    predict.fit <- NULL
  } # is.null point.est
  
  if(!is.null(type)){
    fitted.matrix$y.hat <- ifelse(fitted.matrix$y.hat %% 1 < 0.5, floor(fitted.matrix$y.hat), ceiling(fitted.matrix$y.hat))
    fitted.matrix$y.hat <- pmin(max(original.DV), pmax(min(original.DV), fitted.matrix$y.hat))
    if(!is.null(predict.fit)){
      predict.fit <- ifelse(predict.fit %% 1 < 0.5, floor(predict.fit), ceiling(predict.fit))
      predict.fit <- pmin(max(original.DV), pmax(min(original.DV), predict.fit))
    }
  }
  
  rhs.partitions <- as.data.frame(reg.points.matrix)
  fitted.matrix$residuals <-   fitted.matrix$y.hat - original.DV
  
  if(!is.null(type) && type=="class"){
    R2 <- as.numeric(format(mean(fitted.matrix$y.hat==fitted.matrix$y), digits = 4))
  } else {
    y.mean <- mean(fitted.matrix$y)
    R2 <- (sum((fitted.matrix$y - y.mean)*(fitted.matrix$y.hat - y.mean))^2)/(sum((fitted.matrix$y - y.mean)^2)*sum((fitted.matrix$y.hat - y.mean)^2))
  }
  
  lower.pred.int <- NULL
  upper.pred.int <- NULL
  pred.int <- NULL
  
  if(is.numeric(confidence.interval)){
    ci <- abs(UPM.VaR((1-confidence.interval)/2, degree = 1, fitted.matrix$residuals))
    fitted.matrix$conf.int.pos <- fitted.matrix$y.hat + ci
    fitted.matrix$conf.int.neg <- fitted.matrix$y.hat - ci

    if(!is.null(point.est)){
      lower.pred.int = predict.fit - ci
      upper.pred.int = ci + predict.fit

      pred.int = data.frame(lower.pred.int, upper.pred.int)
    }
  }
  
  ### 3d plot
  if(plot && n == 2){
    .nns_require_rgl()
    rgl::plot3d(x = original.IVs[ , 1], y = original.IVs[ , 2], z = original.DV, box = FALSE, size = 3, col='steelblue', xlab = colnames(reg.points.matrix)[1], ylab = colnames(reg.points.matrix)[2], zlab = y.label )

    if(plot.regions){
      region.yhat.fun <- switch(noise.reduction,
                                "mean"       = mean,
                                "median"     = median,
                                "mode"       = mode,
                                "mode_class" = mode,
                                gravity)   # "off" / default
      for(id in unique(NNS.ID)){
        m   <- NNS.ID == id
        x1v <- original.IVs[m, 1]; x2v <- original.IVs[m, 2]
        min.x1 <- min(x1v); max.x1 <- max(x1v)
        min.x2 <- min(x2v); max.x2 <- max(x2v)
        yh     <- region.yhat.fun(original.DV[m])
        rgl::quads3d(x = c(min.x1, min.x1, max.x1, max.x1),
                     y = c(min.x2, max.x2, max.x2, min.x2),
                     z = c(yh, yh, yh, yh), col = "pink", alpha = 1)
        if(identical(min.x1, max.x1) || identical(min.x2, max.x2)){
          rgl::segments3d(x = c(min.x1, max.x1),
                          y = c(min.x2, max.x2),
                          z = c(yh, yh), col = "pink", alpha = 1)
        }
      }
    }#plot.regions = T

    rgl::points3d(x = as.numeric(REGRESSION.POINT.MATRIX[[1]]), y = as.numeric(REGRESSION.POINT.MATRIX[[2]]), z = as.numeric(REGRESSION.POINT.MATRIX[[3]]), col = 'red', size = 5)
    if(!is.null(point.est)){
      if(is.null(np)){
        rgl::points3d(x = point.est[1], y = point.est[2], z = predict.fit, col = 'green', size = 5)
      } else {
        rgl::points3d(x = point.est[,1], y = point.est[,2], z = predict.fit, col = 'green', size = 5)
      }
    }
  }
  
  ### Residual plot
  if(residual.plot){
    resids <- cbind(original.DV, y.hat)
    r2.leg <- bquote(bold(R ^ 2 == .(format(R2, digits = 4))))
    if(!is.null(type) && type=="class") r2.leg <- paste("Accuracy: ", R2)
    plot(seq_along(original.DV), original.DV, pch = 1, lwd = 2, col = "steelblue", xlab = "Index", ylab = expression(paste("y (blue)   ", hat(y), " (red)")), cex.lab = 1.5, mgp = c(2, .5, 0))
    lines(seq_along(fitted.matrix$y.hat), fitted.matrix$y.hat, col = 'red', lwd = 2, lty = 1)
    
    if(is.numeric(confidence.interval)){
      polygon(c(seq_along(y.hat), rev(seq_along(y.hat))), c(na.omit(fitted.matrix$conf.int.pos), rev(na.omit(fitted.matrix$conf.int.neg))),
              col = rgb(1, 192/255, 203/255, alpha = 0.375),
              border = NA)
    }
    
    title(main = paste0("NNS Order = multiple"), cex.main = 2)
    legend(location, legend = r2.leg, bty = 'n')
  }
  
  ### Return Values
  if(return.values){
    return(list(R2 = R2,
                rhs.partitions = .NNS.df(rhs.partitions),
                RPM = .NNS.df(REGRESSION.POINT.MATRIX) ,
                Point.est = predict.fit,
                pred.int = .NNS.df(pred.int),
                Fitted.xy = .NNS.df(fitted.matrix)))
  } else {
    invisible(list(R2 = R2,
                   rhs.partitions = .NNS.df(rhs.partitions),
                   RPM = .NNS.df(REGRESSION.POINT.MATRIX),
                   Point.est = predict.fit,
                   pred.int = .NNS.df(pred.int),
                   Fitted.xy = .NNS.df(fitted.matrix)))
  }
}