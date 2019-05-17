### Mode of a distribution
mode = function(x){
  if(length(na.omit(x)) > 1){
    d <- density(na.omit(x))
    d$x[which.max(d$y)]
  } else {
    x
  }
}

### Factor to dummy variable
factor_2_dummy = function(x){
  if(class(x) == "factor"){
    output = model.matrix(~x -1, x)[,-1]
  } else {
    output = x
  }
  output
}

