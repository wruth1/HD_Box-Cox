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

### Computes the gradient of the profile likelihood
### DO NOT USE OUTSIDE OF OTHER METHODS!!!!!!!!!!!!!!!!!!
### Does not adjust Z or Z.hat to match gamma
grad.prof.lik.formula.ls = function(X, Y, Z, Z.hat, gamma) {
  r = Z - Z.hat
  g0 = grad.BC.gamma(Y, gamma)
  g1 = t(X) %*% g0
  g2 = solve(t(X) %*% X, g1)
  g3 = X %*% g2
  g = g0 - g3
  J = sum(log(Y))
  a = t(r) %*% g
  b = t(r) %*% r
  output = -n*(a / b) + J
  return(output)
}

### Computes the BC transformation of Y with the specified gamma,
### Fits a LS model, then returns the profile likelihood and its
### gradient
prof.lik.ls = function(gamma, X, Y, grad = F) {
  Z = BC(Y, gamma)
  data = data.frame(X, Z)
  
  fit = lm(Z ~ X, data = data)
  Z.hat = predict(fit, data.frame(X))
  profile.lik = profile.lik.formula(Y, Z, Z.hat, gamma)
  
  if (!grad) {
    return(profile.lik)
  } else{
    grad.profile.lik = grad.prof.lik.formula.ls(X, Y, Z, 
                                                Z.hat, gamma)
    
    output = list(lik = profile.lik, grad = grad.profile.lik)
    return(output)
  }
}

### Computes the BC transformation of Y with the specified gamma,
### fits a lasso model and gets the profile likelihood at the
### specified lambda(s)
### Optionally also computes the gradient of the prof lik 
### at gamma and returns a list with the lik and the grad
prof.lik.lasso = function(gamma, X, Y, lambda, grad=F) {
  Z = BC(Y, gamma)
  fit = glmnet(X, Z)
  Z.hat = predict(fit, X, s = lambda)
  
  profile.likelihoods = sapply(seq_along(lambda), function(i) {
    this.Z.hat = Z.hat[, i]
    this.lik = profile.lik.formula(Y, Z, this.Z.hat, gamma)
    return(this.lik)
  })
  
  if(!grad){
    return(profile.likelihoods)
  } else{
    gradients = sapply(seq_along(lambda), function(i){
      this.Z.hat = Z.hat[, i]
      this.grad = grad.prof.lik.formula(Y, Z, this.Z.hat, gamma)
      return(this.grad)
    })
    output = list(lik = profile.likelihoods,
                  grad = gradients)
    return(output)
  }
}



