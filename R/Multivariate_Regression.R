NNS.M.reg <- function (X_n, Y, factor.2.dummy = TRUE, order = NULL, stn = NULL, n.best = NULL, type = NULL, point.est = NULL, point.only = FALSE,
                       plot = FALSE, residual.plot = TRUE, location = NULL, noise.reduction = 'off', dist = "L2",
                       return.values = FALSE, plot.regions = FALSE, ncores = NULL){
  
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
  
  minimums <- apply(original.IVs, 2, min)
  maximums <- apply(original.IVs, 2, max)
  
  reg.points <- list()
  sections <- list()
  
  if(is.null(order)) order <- floor(max(1,(NNS.copula(original.matrix)* 10)))
  
  
  ###  Regression Point Matrix
  if(is.numeric(order)){
    reg.points <- apply(original.IVs, 2, function(b) NNS.reg(b, original.DV, factor.2.dummy = factor.2.dummy, order = order, stn = stn, type = type, noise.reduction = noise.reduction, plot = FALSE, multivariate.call = TRUE, ncores = 1)$x)
    
    if(length(unique(sapply(reg.points, length))) != 1){
      reg.points.matrix <- do.call('cbind', lapply(reg.points, `length<-`, max(lengths(reg.points))))
    } else {
      reg.points.matrix <- reg.points
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
    colnames.list <- list()
    for(i in 1 : n){
      colnames.list[i] <- paste0("X", i)
    }
    colnames(reg.points.matrix) <- as.character(colnames.list)
  }
  
  if(is.numeric(order) || is.null(order)){
    reg.points.matrix <- unique(reg.points.matrix)
  }
  
  if(!is.null(order) && order=="max" && is.null(n.best)) n.best <- 1
  
  ### Find intervals in regression points for each variable, use left.open T and F for endpoints.
  NNS.ID <- list()
  
  ### PARALLEL
  
  if(is.null(ncores)) {
    num_cores <- as.integer(parallel::detectCores()) - 1
  } else {
    num_cores <- ncores
  }
  
  if(num_cores<=1){
    for(j in 1:n){
      sorted.reg.points <- na.omit(sort(reg.points.matrix[ , j]))
      NNS.ID[[j]] <- findInterval(original.IVs[ , j], vec = sorted.reg.points, left.open = FALSE)
    }
  } else {
    cl <- parallel::makeCluster(num_cores)
    doParallel::registerDoParallel(cl)
    invisible(data.table::setDTthreads(1))
    NNS.ID <- foreach(j = 1:n)%dopar%{
      sorted.reg.points <- na.omit(sort(reg.points.matrix[ , j]))
      return(findInterval(original.IVs[ , j], vec = sorted.reg.points, left.open = FALSE))
    }
  }
  
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
  
  if(!is.null(type)){
    y.hat <- ifelse(y.hat %% 1 < 0.5, floor(y.hat), ceiling(y.hat))
  }
  
  
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
  
  if(is.null(n.best)) n.best <- floor(sqrt(n))
  
  if(n.best > 1 && !point.only){
    if(num_cores > 1){
      fitted.matrix$y.hat <- parallel::parApply(cl, original.IVs, 1, function(z) NNS::NNS.distance(rpm = REGRESSION.POINT.MATRIX, rpm_class = RPM_CLASS, dist.estimate = z, type = dist, k = n.best, class = type)[1])
    } else {
      fits <- data.table::data.table(original.IVs)
      
      fits <- fits[, DISTANCES :=  NNS::NNS.distance(rpm = REGRESSION.POINT.MATRIX, rpm_class = RPM_CLASS, dist.estimate = .SD, type = dist, k = n.best, class = type)[1], by = 1:nrow(original.IVs)]
      
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
      DISTANCES <- list()

      if(num_cores>1){
        DISTANCES <- parallel::parApply(cl, distances, 1, function(z) NNS::NNS.distance(rpm = REGRESSION.POINT.MATRIX, rpm_class = RPM_CLASS, dist.estimate = z, type = dist, k = n.best, class = type)[1])
        
        parallel::stopCluster(cl)
        registerDoSEQ()
        invisible(data.table::setDTthreads(0, throttle = NULL))
      } else {
        distances <- data.table::data.table(point.est)
        distances <- distances[, DISTANCES :=  NNS::NNS.distance(rpm = REGRESSION.POINT.MATRIX, rpm_class = RPM_CLASS, dist.estimate = .SD, type = dist, k = n.best, class = type)[1], by = 1:nrow(point.est)]
        
        DISTANCES <- as.numeric(unlist(distances$DISTANCES))
      }

      if(any(outsiders>0)){
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
    R2 <- format(mean(fitted.matrix$y.hat==original.DV), digits = 4) 
  } else {
    R2num <- sum((fitted.matrix$y.hat - mean(original.DV))*(original.DV - mean(original.DV)))^ 2
    R2den <- sum((fitted.matrix$y.hat - mean(original.DV)) ^ 2) * sum((original.DV - mean(original.DV)) ^ 2)
    R2 <- max(0, min(1, R2num / R2den ))
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
    matplot(resids, type = 'l', xlab = "Index", ylab = expression(paste("y (black)   ", hat(y), " (red)")), cex.lab = 1.5, mgp = c(2, .5, 0))
    
    title(main = paste0("NNS Order = multiple"), cex.main = 2)
    legend(location, legend = r2.leg, bty = 'n')
  }
  
  ### Return Values
  if(return.values){
    return(list(R2 = R2,
                rhs.partitions = rhs.partitions,
                RPM = REGRESSION.POINT.MATRIX[] ,
                Point.est = predict.fit,
                Fitted.xy = fitted.matrix[]))
  } else {
    invisible(list(R2 = R2,
                   rhs.partitions = rhs.partitions,
                   RPM = REGRESSION.POINT.MATRIX[],
                   Point.est = predict.fit,
                   Fitted.xy = fitted.matrix[]))
  }
  
}
