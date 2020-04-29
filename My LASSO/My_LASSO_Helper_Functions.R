### Soft-thresholds the number x at the value t
st = function(x, t){
  a = abs(x)
  if(a <= t){
    return(0)
  } else{
    b = a - t
    c = b * sign(x)
    return(c)
  }
}

### Compute the smallest lambda value that returns a null model
get.lambda.max = function(X, Y){
  corrs = apply(X, 2, function(x){
    return(x %*% Y)
  })
  return(max(abs(corrs)))
}

### Fit a LASSO model to predict Y using X with parameter l
### Optionally, start coordinate descent at beta.0
one.lasso = function(X, Y, l, beta.0 = rep(0, times = ncol(X))){
  ### Iterate through coordinate descent,
  ### possibly using beta.0 as a warm-start
  return(0)
}

one.beta = function()