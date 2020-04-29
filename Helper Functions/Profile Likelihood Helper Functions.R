####################
### LS Functions ###
####################

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
prof.lik.ls = function(gamma, X, Y) {
  Z = BC(Y, gamma)
  data = data.frame(X, Z)
  
  fit = lm(Z ~ X, data = data)
  Z.hat = predict(fit, data.frame(X))
  profile.lik = profile.lik.formula(Y, Z, Z.hat, gamma)
  
  return(profile.lik)
}

### Computes the gradient of the profile likelihood for gamma
### in an LS model. Redundant, but necessary for "optimx()"
grad.prof.lik.ls = function(gamma, X, Y){
  Z = BC(Y, gamma)
  data = data.frame(X, Z)
  fit = lm(Z ~ X, data = data)
  Z.hat = predict(fit, data.frame(X))
  grad.profile.lik = grad.prof.lik.formula.ls(X, Y, Z, 
                                              Z.hat, gamma)
  return(grad.profile.lik)
}


### Finding the root of this function solves 
### the equation log-lik = val
lik.root.ls = function(gamma, val, X, Y){
  this.lik = prof.lik.ls(gamma, X, Y)
  to.root = this.lik - val
  return(to.root)
}






#######################
### LASSO Functions ###
#######################

### Computes the BC transformation of Y with the specified gamma,
### fits a lasso model and gets the profile likelihood at the
### specified lambda(s)
prof.lik.lasso = function(gamma, X, Y, lambda, penal = F) {
  Z = BC(Y, gamma)
  fit = glmnet(X, Z)
  Z.hat = predict(fit, X, s = lambda)
  
  profile.likelihoods = sapply(seq_along(lambda), function(i) {
    this.Z.hat = Z.hat[, i]
    this.lik = profile.lik.formula(Y, Z, this.Z.hat, gamma)
    return(this.lik)
  })
  
  # Subtract L1 penalty from profile likelihood if requested
  if(penal){
    beta.hats = predict(fit, type = "coefficients", s = lambda)[-1,]
    l1.pens = sapply(seq_along(lambda), function(j){
      this.beta.hat = beta.hats[,j]
      this.norm = sum(abs(this.beta.hat))
      this.l = lambda[j]
      this.pen = this.l * this.norm
      return(this.pen)
    })
    profile.likelihoods = profile.likelihoods - l1.pens
  }
  
  return(profile.likelihoods)
  
}

### Computes the profile likelihood for gamma, with beta fit using CV lasso
### penal controls whether the l1 penalty should be added to the likelihood
prof.lik.CV.lasso = function(gamma, X, Y, penal = F, folds = NULL){
  Z = BC(Y, gamma)
  fit = cv.glmnet(X, Z, foldid=folds)
  Z.hat = predict(fit, X, s = "lambda.1se")
  
  lik = profile.lik.formula(Y, Z, Z.hat, gamma)
  
  if(penal){
    l = fit$lambda.1se
    beta.hat = predict(fit, type = "coefficients", s = "lambda.1se")[-1]
    l1.pen = l * sum(abs(beta.hat))
    lik = lik - l1.pen
  }
  return(lik)
}

### Finding the root of this function solves 
### the equation log-lik = val
lik.root.CV.lasso = function(gamma, val, X, Y, folds = NULL){
  this.lik = prof.lik.CV.lasso(gamma, X, Y, folds = folds)
  to.root = this.lik - val
  return(to.root)
}


### Finding the root of this function solves 
### the equation log-lik = val
lik.root.lasso = function(gamma, val, X, Y){
  this.lik = prof.lik.lasso(gamma, X, Y)
  to.root = this.lik - val
  return(to.root)
}

### Computes the BC transformation of Y with the specified gamma,
### fits a lasso model and gets the profile likelihood at the
### specified lambda(s), then subtracts an l1 penalty in each beta
penal.prof.lik.lasso = function(gamma, X, Y, lambda) {
  Z = BC(Y, gamma)
  fit = glmnet(X, Z)
  Z.hat = predict(fit, X, s = lambda)
  beta.hat = predict(fit, X, s = lambda, type = "coefficients")

  profile.likelihoods = sapply(seq_along(lambda), function(i) {
    this.Z.hat = Z.hat[, i]
    this.beta.hat = beta.hat[-1,i]
    this.lik = profile.lik.formula(Y, Z, this.Z.hat, gamma)
    this.penal = lambda[i] * sum(abs(this.beta.hat)) / sd(Z)
    return(this.lik - this.penal)
  })
  
  return(profile.likelihoods)
  
}

