NNS.M.reg <- function (X_n, Y, factor.2.dummy = TRUE, order = NULL, stn = NULL, n.best = NULL, type = NULL, point.est = NULL, point.only = FALSE,
                       plot = FALSE, residual.plot = TRUE, location = NULL, noise.reduction = 'off', dist = "L2",
                       return.values = FALSE, plot.regions = FALSE, ncores = NULL, confidence.interval = NULL){
  
  dist <- tolower(dist)
  
  
  ### For Multiple regressions
  ###  Turn each column into numeric values
  original.IVs <- X_n
  original.DV <- Y
  n <- ncol(original.IVs)
  
  if(is.null(ncol(X_n))) X_n <- t(t(X_n))
  
  if(is.null(names(Y))){
    y.label <- "Y"
  } else {
    y.label <- names(Y)
  }
  
  np <- nrow(point.est)
  
  if(is.null(np) & !is.null(point.est)){
    point.est <- t(point.est)
  } else {
    point.est <- point.est
  }
  
  if(!is.null(point.est)){
    if(ncol(point.est) != n){
      stop("Please ensure 'point.est' is of compatible dimensions to 'x'")
    }
  }
  
  original.matrix <- cbind.data.frame(original.IVs, original.DV)
  norm.matrix <- apply(original.matrix[,1:n], 2, function(z) NNS.rescale(z, 0, 1))
  
  minimums <- apply(original.IVs, 2, min)
  maximums <- apply(original.IVs, 2, max)
  
  dependence <- max(c(NNS.copula(original.matrix), NNS.copula(cbind(norm.matrix, original.DV))))
  
  if(is.null(order)) order <- max(1, ceiling(dependence*10))


  ###  Regression Point Matrix
  if(is.numeric(order)){
    reg.points <- lapply(1:ncol(original.IVs), function(b) NNS.reg(original.IVs[, b], original.DV, factor.2.dummy = factor.2.dummy, order = order, stn = stn, type = type, noise.reduction = noise.reduction, plot = FALSE, multivariate.call = TRUE, ncores = 1)$x)
    
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
  
  if(is.null(colnames(original.IVs))){
    colnames.list <- lapply(1 : ncol(original.IVs), function(i) paste0("x", i))
    colnames(reg.points.matrix) <- as.character(colnames.list)
  }
  
  if(is.numeric(order) || is.null(order)) reg.points.matrix <- unique(reg.points.matrix)

  
  if(!is.null(order) && order=="max" && is.null(n.best)) n.best <- 1
  
  ### Find intervals in regression points for each variable, use left.open T and F for endpoints.
  ### PARALLEL
  
  if(is.null(ncores)) {
    num_cores <- as.integer(max(2L, parallel::detectCores(), na.rm = TRUE)) - 1
  } else {
    num_cores <- ncores
  }
  
  if(num_cores > 1){
    cl <- parallel::makeCluster(num_cores)
    doParallel::registerDoParallel(cl)
    invisible(data.table::setDTthreads(1))
  } else {
    foreach::registerDoSEQ()
    invisible(data.table::setDTthreads(0, throttle = NULL))
  }

  NNS.ID <- lapply(1:n, function(j) findInterval(original.IVs[ , j], vec = na.omit(sort(reg.points.matrix[ , j])), left.open = FALSE))

  NNS.ID <- do.call(cbind, NNS.ID)
  
  ### Create unique identifier of each observation's interval
  NNS.ID <- gsub(do.call(paste, as.data.frame(NNS.ID)), pattern = " ", replacement = ".")
  
  
  ### Match y to unique identifier
  obs <- c(1 : length(Y))
  
  mean.by.id.matrix <- data.table::data.table(original.IVs, original.DV, NNS.ID, obs)
  
  data.table::setkey(mean.by.id.matrix, 'NNS.ID', 'obs')
  
  if(is.numeric(order) || is.null(order)){
    if(noise.reduction == 'off'){
      mean.by.id.matrix <- mean.by.id.matrix[ , c(paste("RPM", 1:n), "y.hat") := lapply(.SD, function(z) gravity(as.numeric(z))), .SDcols = seq_len(n+1) ,by = 'NNS.ID']
    }
    if(noise.reduction == 'mean'){
      mean.by.id.matrix <- mean.by.id.matrix[ , c(paste("RPM", 1:n), "y.hat") := lapply(.SD, function(z) mean(as.numeric(z))), .SDcols = seq_len(n+1), by = 'NNS.ID']
    }
    if(noise.reduction == 'median'){
      mean.by.id.matrix <- mean.by.id.matrix[ , c(paste("RPM", 1:n), "y.hat") := lapply(.SD, function(z) median(as.numeric(z))), .SDcols = seq_len(n+1), by = 'NNS.ID']
    }
    if(noise.reduction == 'mode'){
      mean.by.id.matrix <- mean.by.id.matrix[ , c(paste("RPM", 1:n), "y.hat") := lapply(.SD, function(z) mode(as.numeric(z))), .SDcols = seq_len(n+1), by = 'NNS.ID']
    }
    if(noise.reduction == 'mode_class'){
      mean.by.id.matrix <- mean.by.id.matrix[ , c(paste("RPM", 1:n), "y.hat") := lapply(.SD, function(z) mode_class(as.numeric(z))), .SDcols = seq_len(n+1), by = 'NNS.ID']
    }
  } else {
    mean.by.id.matrix <- mean.by.id.matrix[ , c(paste("RPM", 1:n), "y.hat") := .SD , .SDcols = seq_len(n+1), by = 'NNS.ID']
  }
  
  
  ###Order y.hat to order of original Y
  resid.plot <- mean.by.id.matrix[]
  data.table::setkey(resid.plot, 'obs')
  
  
  y.hat <- unlist(mean.by.id.matrix[ , .(y.hat)])
  
  if(!is.null(type)) y.hat <- ifelse(y.hat %% 1 < 0.5, floor(y.hat), ceiling(y.hat))
  
  
  fitted.matrix <- data.table::data.table(original.IVs, y = original.DV, y.hat, mean.by.id.matrix[ , .(NNS.ID)])
  
  fitted.matrix$residuals <- fitted.matrix$y - fitted.matrix$y.hat
  fitted.matrix[, bias := gravity(residuals),  by = NNS.ID]
  fitted.matrix$y.hat <- fitted.matrix$y.hat - fitted.matrix$bias
  fitted.matrix$bias <- NULL
  
  
  data.table::setkey(mean.by.id.matrix, 'NNS.ID')
  REGRESSION.POINT.MATRIX <- mean.by.id.matrix[ , c("obs") := NULL]
  
  REGRESSION.POINT.MATRIX <- REGRESSION.POINT.MATRIX[, .SD[1], by = NNS.ID]
  
  
  REGRESSION.POINT.MATRIX <- REGRESSION.POINT.MATRIX[, .SD, .SDcols = colnames(mean.by.id.matrix)%in%c(paste("RPM", 1:n), "y.hat")]
  
  data.table::setnames(REGRESSION.POINT.MATRIX, 1:n, colnames(mean.by.id.matrix)[1:n])
  
  RPM_CLASS <- apply(do.call(cbind, lapply(REGRESSION.POINT.MATRIX[ ,1:n], FUN = function(z) ifelse(z%%1 < .5, floor(z), ceiling(z)))), 2, as.integer)
  
  if(is.null(n.best)) n.best <- max(1, floor((1-dependence)*sqrt(n)))

  
  if(n.best > 1 && !point.only){
    if(num_cores > 1){
      fitted.matrix$y.hat <- parallel::parApply(cl, original.IVs, 1, function(z) NNS.distance(rpm = REGRESSION.POINT.MATRIX, rpm_class = RPM_CLASS, dist.estimate = z, type = dist, k = n.best, class = type)[1])
    } else {
      fits <- data.table::data.table(original.IVs)
      
      fits <- fits[, DISTANCES :=  NNS.distance(rpm = REGRESSION.POINT.MATRIX, rpm_class = RPM_CLASS, dist.estimate = .SD, type = dist, k = n.best, class = type)[1], by = 1:nrow(original.IVs)]
      
      fitted.matrix$y.hat <- as.numeric(unlist(fits$DISTANCES))
    }
    
    y.hat <- fitted.matrix$y.hat
    
    if(!is.null(type)) y.hat <- ifelse(y.hat %% 1 < 0.5, floor(y.hat), ceiling(y.hat))
  }
  
  
  
  ### Point estimates
  if(!is.null(point.est)){
    
    ### Point estimates
    central.points <- apply(REGRESSION.POINT.MATRIX[, .SD, .SDcols = 1:n], 2, function(x) gravity(x))
    
    predict.fit <- numeric()
    
    outsiders <- point.est<minimums | point.est>maximums
    outsiders[is.na(outsiders)] <- 0
    
    if(is.null(np)){
      l <- length(point.est)
      
      if(!any(outsiders)){
        predict.fit <- NNS::NNS.distance(rpm = REGRESSION.POINT.MATRIX, rpm_class = RPM_CLASS, dist.estimate = point.est, type = dist, k = n.best, class = type)
      } else {
        boundary.points <- pmin(pmax(point.est, minimums), maximums)
        mid.points <- (boundary.points + central.points) / 2
        mid.points_2 <- (boundary.points + mid.points) / 2
        last.known.distance_1 <- sqrt(sum((boundary.points - central.points) ^ 2))
        last.known.distance_2 <- sqrt(sum((boundary.points - mid.points) ^ 2))
        last.known.distance_3 <- sqrt(sum((boundary.points - mid.points_2) ^ 2))
        
        boundary.estimates <- NNS::NNS.distance(rpm = REGRESSION.POINT.MATRIX, rpm_class = RPM_CLASS, dist.estimate = boundary.points, type = dist, k = n.best, class = type)
        
        last.known.gradient_1 <- (boundary.estimates - NNS::NNS.distance(rpm = REGRESSION.POINT.MATRIX, rpm_class = RPM_CLASS, dist.estimate = central.points, type = dist, k = n.best, class = type)) / last.known.distance_1
        last.known.gradient_2 <- (boundary.estimates - NNS::NNS.distance(rpm = REGRESSION.POINT.MATRIX, rpm_class = RPM_CLASS, dist.estimate = mid.points, type = dist, k = n.best, class = type)) / last.known.distance_2
        last.known.gradient_3 <- (boundary.estimates - NNS::NNS.distance(rpm = REGRESSION.POINT.MATRIX, rpm_class = RPM_CLASS, dist.estimate = mid.points_2, type = dist, k = n.best, class = type)) / last.known.distance_3
        
        last.known.gradient <- (last.known.gradient_1*3 + last.known.gradient_2*2 + last.known.gradient_3) / 6
        
        last.distance <- sqrt(sum((point.est - boundary.points) ^ 2))
        
        predict.fit <- last.distance * last.known.gradient + boundary.estimates
      }
    }
    
    if(!is.null(np)){
      DISTANCES <- vector(mode = "list", np)
      distances <- data.table::data.table(point.est)
      if(num_cores > 1){
        DISTANCES <- parallel::parApply(cl, distances, 1, function(z) NNS.distance(rpm = REGRESSION.POINT.MATRIX, rpm_class = RPM_CLASS, dist.estimate = z, type = dist, k = n.best, class = type)[1])
        
        parallel::stopCluster(cl)
        rm(cl)
        foreach::registerDoSEQ()
        invisible(data.table::setDTthreads(0, throttle = NULL))
      } else {
        distances <- distances[, DISTANCES :=  NNS.distance(rpm = REGRESSION.POINT.MATRIX, rpm_class = RPM_CLASS, dist.estimate = .SD, type = dist, k = n.best, class = type)[1], by = 1:nrow(point.est)]
        
        DISTANCES <- as.numeric(unlist(distances$DISTANCES))
      }

      if(any(outsiders > 0)){
        outsiders <- rowSums(outsiders)
        outside.index <- as.numeric(which(outsiders>0))
        
        for(i in outside.index){
          outside.points <- point.est[i,]
          boundary.points <- pmin(pmax(outside.points, minimums), maximums)
          mid.points <- (boundary.points + central.points) / 2
          mid.points_2 <- (boundary.points + mid.points) / 2
          last.known.distance_1 <- sqrt(sum((boundary.points - central.points) ^ 2))
          last.known.distance_2 <- sqrt(sum((boundary.points - mid.points) ^ 2))
          last.known.distance_3 <- sqrt(sum((boundary.points - mid.points_2) ^ 2))
          
          boundary.estimates <- NNS::NNS.distance(rpm = REGRESSION.POINT.MATRIX, rpm_class = RPM_CLASS,
                                                  dist.estimate = boundary.points,
                                                  type = dist, k = n.best, class = type)
          
          last.known.gradient_1 <- (boundary.estimates - NNS::NNS.distance(rpm = REGRESSION.POINT.MATRIX, rpm_class = RPM_CLASS, dist.estimate = central.points, type = dist, k = n.best, class = type)) / last.known.distance_1
          last.known.gradient_2 <- (boundary.estimates - NNS::NNS.distance(rpm = REGRESSION.POINT.MATRIX, rpm_class = RPM_CLASS, dist.estimate = mid.points, type = dist, k = n.best, class = type)) / last.known.distance_2
          last.known.gradient_3 <- (boundary.estimates - NNS::NNS.distance(rpm = REGRESSION.POINT.MATRIX, rpm_class = RPM_CLASS, dist.estimate = mid.points_2, type = dist, k = n.best, class = type)) / last.known.distance_3
          
          last.known.gradient <- (last.known.gradient_1*3 + last.known.gradient_2*2 + last.known.gradient_3) / 6
          
          last.distance <- sqrt(sum((outside.points - boundary.points) ^ 2))
          
          
          DISTANCES[i] <- last.distance * last.known.gradient + boundary.estimates
        }
      }
      
      predict.fit <- DISTANCES
      
      if(point.only) return(list(Point.est = predict.fit,  RPM = REGRESSION.POINT.MATRIX[] ))
    }
    
  } else {
    predict.fit <- NULL
  } # is.null point.est
  
  if(!is.null(type)){
    fitted.matrix$y.hat <- ifelse(fitted.matrix$y.hat %% 1 < 0.5, floor(fitted.matrix$y.hat), ceiling(fitted.matrix$y.hat))
    if(!is.null(predict.fit)) predict.fit <- ifelse(predict.fit %% 1 < 0.5, floor(predict.fit), ceiling(predict.fit))
  }
  
  rhs.partitions <- data.table::data.table(reg.points.matrix)
  fitted.matrix$residuals <-  original.DV - fitted.matrix$y.hat
  
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
    fitted.matrix[, `:=` ( 'conf.int.pos' = abs(UPM.VaR((1-confidence.interval)/2, degree = 1, residuals)) + y.hat)]
    fitted.matrix[, `:=` ( 'conf.int.neg' = y.hat - abs(LPM.VaR((1-confidence.interval)/2, degree = 1, residuals)))]
    
    if(!is.null(point.est)){
      lower.pred.int = predict.fit - abs(LPM.VaR((1-confidence.interval)/2, degree = 1, fitted.matrix$residuals))
      upper.pred.int = abs(UPM.VaR((1-confidence.interval)/2, degree = 1, fitted.matrix$residuals)) + predict.fit
      
      pred.int = data.table::data.table(lower.pred.int, upper.pred.int)
    }
  }
  
  ### 3d plot
  if(plot && n == 2){
    region.1 <- mean.by.id.matrix[[1]]
    region.2 <- mean.by.id.matrix[[2]]
    region.3 <- mean.by.id.matrix[ , y.hat]
    
    rgl::plot3d(x = original.IVs[ , 1], y = original.IVs[ , 2], z = original.DV, box = FALSE, size = 3, col='steelblue', xlab = colnames(reg.points.matrix)[1], ylab = colnames(reg.points.matrix)[2], zlab = y.label )
    
    if(plot.regions){
      region.matrix <- data.table::data.table(original.IVs, original.DV, NNS.ID)
      region.matrix[ , `:=` (min.x1 = min(.SD), max.x1 = max(.SD)), by = NNS.ID, .SDcols = 1]
      region.matrix[ , `:=` (min.x2 = min(.SD), max.x2 = max(.SD)), by = NNS.ID, .SDcols = 2]
      if(noise.reduction == 'off'){
        region.matrix[ , `:=` (y.hat = gravity(original.DV)), by = NNS.ID]
      }
      if(noise.reduction =="mean"){
        region.matrix[ , `:=` (y.hat = mean(original.DV)), by = NNS.ID]
      }
      if(noise.reduction =="median"){
        region.matrix[ , `:=` (y.hat = median(original.DV)), by = NNS.ID]
      }
      if(noise.reduction=="mode"|| noise.reduction=="mode_class"){
        region.matrix[ , `:=` (y.hat = mode(original.DV)), by = NNS.ID]
      }
      
      data.table::setkey(region.matrix, NNS.ID, min.x1, max.x1, min.x2, max.x2)
      region.matrix[ ,{
        rgl::quads3d(x = .(min.x1[1], min.x1[1], max.x1[1], max.x1[1]),
                     y = .(min.x2[1], max.x2[1], max.x2[1], min.x2[1]),
                     z = .(y.hat[1], y.hat[1], y.hat[1], y.hat[1]), col="pink", alpha=1)
        if(identical(min.x1[1], max.x1[1]) || identical(min.x2[1], max.x2[1])){
          rgl::segments3d(x = .(min.x1[1], max.x1[1]),
                          y = .(min.x2[1], max.x2[1]),
                          z = .(y.hat[1], y.hat[1]), col = "pink", alpha = 1)
        }
      }
      , by = NNS.ID]
    }#plot.regions = T
    
    
    rgl::points3d(x = as.numeric(unlist(REGRESSION.POINT.MATRIX[ , .SD, .SDcols = 1])), y = as.numeric(unlist(REGRESSION.POINT.MATRIX[ , .SD, .SDcols = 2])), z = as.numeric(unlist(REGRESSION.POINT.MATRIX[ , .SD, .SDcols = 3])), col = 'red', size = 5)
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
    lines(seq_along(y.hat), y.hat, col = 'red', lwd = 2, lty = 1)
    
    if(is.numeric(confidence.interval)){
      lines(seq_along(y.hat), na.omit(fitted.matrix$conf.int.pos), col = 'pink')
      lines(seq_along(y.hat), na.omit(fitted.matrix$conf.int.neg), col = 'pink')
      
      polygon(c(seq_along(y.hat), rev(seq_along(y.hat))), c(na.omit(fitted.matrix$conf.int.pos), rev(na.omit(fitted.matrix$conf.int.neg))), col = "pink", border = NA)
    
      points(seq_along(original.DV), original.DV, pch = 1, lwd = 2, col = "steelblue")
      lines(seq_along(y.hat), y.hat, col = 'red', lwd = 2, lty = 1)
    }
    
    title(main = paste0("NNS Order = multiple"), cex.main = 2)
    legend(location, legend = r2.leg, bty = 'n')
  }
  
  
  
  
  ### Return Values
  if(return.values){
    return(list(R2 = R2,
                rhs.partitions = rhs.partitions,
                RPM = REGRESSION.POINT.MATRIX[] ,
                Point.est = predict.fit,
                pred.int = pred.int,
                Fitted.xy = fitted.matrix[]))
  } else {
    invisible(list(R2 = R2,
                   rhs.partitions = rhs.partitions,
                   RPM = REGRESSION.POINT.MATRIX[],
                   Point.est = predict.fit,
                   pred.int = pred.int,
                   Fitted.xy = fitted.matrix[]))
  }
  
}
