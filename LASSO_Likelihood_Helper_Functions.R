profile.lik = function(gamma, X, Y) {
  n = length(Y)
  Y.new = BC(Y, gamma)
  fit = lm(Y.new ~ X)
  mse = sum((fit$residuals) ^ 2) / n
  lik = -n * log(mse) / 2
  Jacob = (gamma - 1) * sum(log(Y))
  lik = lik + Jacob
  return(lik)
}

profile.lik.formula = function(Y.obs, Y.BC, Y.BC.hat, gamma){
  n = length(Y.obs)
  MSE = mean((Y.BC - Y.BC.hat)^2)
  lik = -n * log(MSE) /2
  Jacob = (gamma -1) * sum(log(Y.obs))
  lik = lik + Jacob
  return(lik)
}

### Compute `profile likelihood' at gamma based on LASSO
get.profile.lik = function(gamma, X, Y, fit, lambda) {
  Y.hat = predict(fit, X, s = lambda)
  Y.obs = inv.BC(Y, gamma)
  lik = profile.lik.formula(Y.obs, Y, Y.hat, gamma)
  return(lik)
}

### Computes the BC transformation of Y.obs with the specified gamma,
### fits a lasso model and gets the profile likelihood at the
### specified lambda(s)
### Output: an unnamed vector of the same length as lambda
profile.lik.lasso = function(gamma, X, Y.obs, lambda){
  Y.BC = BC(Y.obs, gamma)
  fit = glmnet(X, Y.BC)
  Y.BC.hat = predict(fit, X, s = lambda)
  
  profile.likelihoods = sapply(seq_along(lambda), function(i){
    this.Y.hat = Y.BC.hat[,i]
    this.lik = profile.lik.formula(Y.obs, Y.BC, this.Y.hat, gamma)
    return(this.lik)
  })
  
  return(profile.likelihoods)
}


### Computes the BC transformation of Y.obs with the specified gamma,
### fits a LS model and gets the profile likelihood
### Output: a number
profile.lik.ls = function(gamma, X, Y.obs){
  Y.BC = BC(Y.obs, gamma)
  
  data = data.frame(X, Y.BC)
  fit = lm(Y.BC ~ X, data=data)
  Y.BC.hat = predict(fit, data.frame(X))
  
  profile.lik = profile.lik.formula(Y.obs, Y.BC, Y.BC.hat, gamma)
  return(profile.lik)
}

### Old version
# profile.lik.lasso = function(gamma, X, Y, folds = NULL) {
#   n = length(Y)
#   Y.new = BC(Y, gamma)
#   if (is.null(folds)) {
#     fit.cv = cv.glmnet(X, Y.new)
#   } else{
#     fit.cv = cv.glmnet(X, Y.new, foldid = folds)
#   }
#   Y.hat.min = predict(fit.cv, X, s = "lambda.min")
#   Y.hat.1se = predict(fit.cv, X, s = "lambda.1se")
#   resid.min = Y.new - Y.hat.min
#   resid.1se = Y.new - Y.hat.1se
#   sse.min = sum((resid.min) ^ 2) / n
#   sse.1se = sum((resid.1se) ^ 2) / n
#   lik.min = -n * log(sse.min) / 2
#   lik.1se = -n * log(sse.1se) / 2
#   Jacob = (gamma - 1) * sum(log(Y))
#   lik.min = lik.min + Jacob
#   lik.1se = lik.1se + Jacob
#   this.liks = c(lik.min, lik.1se)
#   this.lambdas = c(fit.cv$lambda.min, fit.cv$lambda.1se)
#   this.vars = list(
#     predict(fit.cv, X, s = "lambda.min", type = "nonzero"),
#     predict(fit.cv, X, s = "lambda.1se", type = "nonzero")
#   )
#   output = list(this.liks, this.lambdas, this.vars)
#   return(output)
# }

increment.counts = function(counts, inds.list) {
  inds = unlist(inds.list)
  counts[inds] = counts[inds] + 1
  return(counts)
}

#Returns a vector that the BC transform (with par. gamma) maps to Y
#Note: Does not account for negative values
inv.BC = function(Y, gamma, type = "good") {
  if (type == "good") {
    if (gamma == 0) {
      return(exp(Y))
    } else{
      Z = 1 + gamma * Y
      Z = Z ^ (1 / gamma)
      return(Z)
    }
  } else if (type == "simple") {
    if (gamma == 0) {
      return(exp(Y))
    } else{
      Z = Y ^ (1 / gamma)
      return(Z)
    }
  }
}

#Returns the BC transformation of Y with par. gamma
#Note: Does not account for negative values
BC = function(Y, gamma, type = "good") {
  if (type == "good") {
    if (gamma == 0) {
      return(log(Y))
    } else{
      Z = Y ^ gamma
      Z = (Z - 1) / gamma
      return(Z)
    }
  } else if (type == "simple") {
    if (gamma == 0) {
      return(log(Y))
    } else{
      Z = Y ^ gamma
      return(Z)
    }
  }
}

