#' NNS Normalization
#'
#' Normalizes a matrix of variables based on nonlinear scaling normalization method.
#'
#' @param X a numeric matrix or data frame, or a list.
#' @param linear logical; \code{FALSE} (default) Performs a linear scaling normalization, resulting in equal means for all variables.
#' @param chart.type  options: ("l", "b"); \code{NULL} (default).  Set \code{(chart.type = "l")} for line,
#' \code{(chart.type = "b")} for boxplot.
#' @param location Sets the legend location within the plot, per the \code{x} and \code{y} co-ordinates used in base graphics \link{legend}.
#' @return Returns a \link{data.frame} of normalized values.
#' @note Unequal vectors provided in a list will only generate \code{linear=TRUE} normalized values.
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments" (ISBN: 1490523995)
#' @examples
#' \dontrun{
#' set.seed(123)
#' x <- rnorm(100) ; y<-rnorm(100)
#' A <- cbind(x, y)
#' NNS.norm(A)
#' 
#' ### Normalize list of unequal vector lengths
#' 
#' vec1 <- c(1, 2, 3, 4, 5, 6, 7)
#' vec2 <- c(10, 20, 30, 40, 50, 60)
#' vec3 <- c(0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3)
#' 
#' vec_list <- list(vec1, vec2, vec3)
#' NNS.norm(vec_list)
#' }
#' @export

NNS.norm <- function(X,
                     linear = FALSE,
                     chart.type = NULL,
                     location = "topleft"){
  
  if(sum(is.na(X)) > 0) stop("You have some missing values, please address.")
  
  if(any(class(X)%in%c("tbl","data.table"))) X <- as.data.frame(X)

  if(any(class(X)%in%"list")){
    if(sum(diff(sapply(X, length)))>0) linear <- TRUE
    m <- sapply(X, mean)
  } else { 
    X <- apply(X, 2, unlist)
    m <- Rfast::colmeans(X)
  }
  
  
  m[m==0] <- 1e-10
  RG <- m %o% (1 / m)
  
  if(!linear){
    if(length(m) < 10){
      scale.factor <- abs(cor(X))
    } else {
      scale.factor <- abs(NNS.dep(X)$Dependence)
    }
    scales <- Rfast::colmeans(RG * scale.factor)
  } else {
    scales <- Rfast::colmeans(RG)
  }
  

  if(any(class(X)%in%"list")) X_Normalized <- mapply('*', X, scales) else X_Normalized <- t(t(X) * scales)
  
  if(any(class(X_Normalized)%in%"list")) n <- length(X_Normalized) else n <- ncol(X_Normalized)
  
  i <- seq_len(n)
  
  if(any(class(X)%in%"list")){
    if(is.null(names(X))){
      new.names <- list()
      for(i in 1 : n){
        new.names[[i]] <- paste0("x_", i)
      }
      names(X) <- unlist(new.names)
    }
  } else {
    if(is.null(colnames(X))){
      new.names <- list()
      for(i in 1 : n){
        new.names[[i]] <- paste0("x_", i)
      }
      colnames(X) <- unlist(new.names)
    }
  }
     
  if(any(class(X_Normalized)%in%"list")){
    names(X_Normalized) <- paste0(names(X), " Normalized")
  } else {
    labels <- c(colnames(X), paste0(colnames(X), " Normalized"))
    colnames(X_Normalized) <- labels[(n + 1) : (2 * n)]
    rows <- rownames(X_Normalized)
  }
  

  if(!is.null(chart.type) && !any(class(X)%in%"list")){
    left_label_size <- max(strwidth(cbind(X, X_Normalized), units = "inches"))*3
    bottom_label_size <- max(strwidth(colnames(X_Normalized), units = "inches"))*8
    
    original.par <- par(no.readonly = TRUE)
    if(chart.type == 'b' ){
      par(mar = c(bottom_label_size, left_label_size, 1, 1))
      boxplot(cbind(X, X_Normalized), las = 2, names = labels, col = c(rep("grey", n), rainbow(n)))
    }
    
    if(chart.type == 'l' ){
      par(mfrow = c(2, 1))   
      par(mar = c(ifelse((class(rows)!="numeric" || !is.null(rows)),4,2), left_label_size , 1, 1))
      
      matplot(X, type = 'l', col = c('steelblue', rainbow(n)), ylab = '', xaxt = 'n', lwd = 2, las = 1)
      legend(location, inset = c(0,0), c(colnames(X)), lty = 1, col = c('steelblue', rainbow(n)), bty = 'n', ncol = floor(n/sqrt(n)), lwd = 2, cex = n/sqrt(n)^exp(1))
      axis(1, at = seq(length(X_Normalized[ , 1]), 1, -floor(sqrt(length(X_Normalized[ , 1])))),
           labels = rownames(X_Normalized[seq(length(X_Normalized[ , 1]), 1, -floor(sqrt(length(X_Normalized[ , 1])))),]), las = 1,
           cex.axis = ifelse((class(rows)!="numeric" || !is.null(rows)),.75,1),
           las = ifelse((class(rows)!="numeric" || !is.null(rows)),3,1),srt=45)
      
      matplot(X_Normalized, type = 'l', col = c('steelblue', rainbow(n)), ylab = '', xaxt = 'n', lwd = 2, las = 1)
      axis(1, at = seq(length(X_Normalized[ , 1]), 1, -floor(sqrt(length(X_Normalized[ , 1])))),
           labels = rownames(X_Normalized[seq(length(X_Normalized[ , 1]), 1, -floor(sqrt(length(X_Normalized[ , 1])))),]), las = 1,
           cex.axis = ifelse((class(rows)!="numeric" || !is.null(rows)),.75,1),
           las = ifelse((class(rows)!="numeric" || !is.null(rows)),3,1),srt=45)
      
      legend(location, c(paste0(colnames(X), " Normalized")), lty = 1, col = c('steelblue', rainbow(n)), bty = 'n', ncol = ceiling(n/sqrt(n)), lwd = 2, cex = n/sqrt(n)^exp(1))
    }
    
    par(original.par)
    
  }
  
  
  
  return(X_Normalized)
  
}
