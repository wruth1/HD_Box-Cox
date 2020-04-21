### Compute the profile likelihood for gamma
### DO NOT USE OUTSIDE OF OTHER METHODS!!!!!!!!!!!!!!!!!!
### Does not adjust Z or Z.hat to match gamma
profile.lik.formula = function(Y, Z, Z.hat, gamma) {
  n = length(Y)
  MSE = mean((Z - Z.hat) ^ 2)
  lik = -n * log(MSE) / 2
  Jacob = (gamma - 1) * sum(log(Y))
  lik = lik + Jacob
  return(lik)
}

### Computes the BC transformation of Y with the specified gamma,
### Fits a LS model, then returns the profile likelihood and its
### gradient
prof.lik.ls = function(gamma, X, Y, grad = F) {
  Z = BC(Y.obs, gamma)
  data = data.frame(X, Z)
  
  fit = lm(Z ~ X, data = data)
  Z.hat = predict(fit, data.frame(X))
  profile.lik = profile.lik.formula(Y, Z, Z.hat, gamma)
  
  if (!grad) {
    return(profile.lik)
  } else{
    grad.profile.lik = grad.prof.lik.formula(Y, Z, Z.hat, gamma)
    
    output = list(lik = profile.lik, grad = grad.profile.lik)
    return(output)
  }
}

### Computes the BC transformation of Y.obs with the specified gamma,
### fits a lasso model and gets the profile likelihood at the
### specified lambda(s)
### Output: an unnamed vector of the same length as lambda
prof.lik.lasso = function(gamma, X, Y.obs, lambda) {
  Y.BC = BC(Y.obs, gamma)
  fit = glmnet(X, Y.BC)
  Y.BC.hat = predict(fit, X, s = lambda)
  
  profile.likelihoods = sapply(seq_along(lambda), function(i) {
    this.Y.hat = Y.BC.hat[, i]
    this.lik = profile.lik.formula(Y.obs, Y.BC, this.Y.hat, gamma)
    return(this.lik)
  })
  
  return(profile.likelihoods)
}


### Computes the gradient of the profile likelihood
### DO NOT USE OUTSIDE OF OTHER METHODS!!!!!!!!!!!!!!!!!!
### Does not adjust Z or Z.hat to match gamma
grad.prof.lik.formula = function(Y, Z, Z.hat, gamma) {
  r = Z - Z.hat
  g = grad.BC.gamma(Y, gamma)
  J = sum(log(Y))
  a = t(r) %*% g
  b = t(r) %*% r
  output = -(a / b) + J
  return(output)
}