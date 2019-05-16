mode = function(x){
  if(length(na.omit(x)) > 1){
    d <- density(na.omit(x))
    d$x[which.max(d$y)]
  } else {
    x
  }
}


factor_2_dummy = function(x){
  if(class(x) == "factor"){
    n=length(unique(x))
    output = model.matrix(~x -1, x)[,-1]
  } else {
    output = x
  }
  output
}

dim.red_factor_2_dummy = function(x){
  if(class(x) == "factor"){
    n=length(unique(x))
    output = model.matrix(~x -1, x)[,-(n+1)]
  } else {
    output = x
  }
  output
}

