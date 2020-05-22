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

### Compute the penalized likelihood for gamma
### DO NOT USE OUTSIDE OF OTHER METHODS!!!!!!!!!!!!!!!!!!
### Does not adjust Z or Z.hat to match gamma
pen.lik.formula = function(Y, Z, Z.hat, gamma, penalty) {
  n = length(Y)
  SSE = sum((Z - Z.hat) ^ 2)
  lik = -n * log(SSE + penalty) / 2
  Jacob = (gamma - 1) * sum(log(Y))
  lik = lik + Jacob
  return(lik)
}



### Computes the BC transformation of Y with the specified gamma,
### fits a lasso model and gets the profile likelihood at the
### specified lambda(s)
prof.lik.lasso = function(gamma, X, Y, lambda) {
  Z = BC(Y, gamma)
  fit = glmnet(X, Z)
  Z.hat = predict(fit, X, s = lambda)
  
  profile.likelihoods = sapply(seq_along(lambda), function(i) {
    this.Z.hat = Z.hat[, i]
    this.lik = profile.lik.formula(Y, Z, this.Z.hat, gamma)
    return(this.lik)
  })
  
  return(profile.likelihoods)
  
}

### Computes the profile likelihood for gamma, with beta fit using CV lasso
### penal controls whether the l1 penalty should be added to the likelihood
prof.lik.CV.lasso = function(gamma, X, Y, folds = NULL,
                             lambda.type = "lambda.1se", all.lambdas=NULL){
  # browser()
  Z = BC(Y, gamma)
  fit = cv.glmnet(X, Z, foldid=folds, lambda = all.lambdas)
  Z.hat = predict(fit, X, s = lambda.type)

  lik = profile.lik.formula(Y, Z, Z.hat, gamma)
  
  return(lik)
}

### Finding the root of this function solves 
### the equation log-lik = val
lik.root.CV.lasso = function(gamma, val, X, Y, folds = NULL, 
                             lambda.type = "lambda.1se",
                             all.lambdas = NULL){
  this.lik = prof.lik.CV.lasso(gamma, X, Y, folds = folds, 
                               lambda.type = lambda.type,
                               all.lambdas = all.lambdas)
  to.root = this.lik - val
  return(to.root)
}


### Finding the root of this function solves 
### the equation log-lik = val
lik.root.lasso = function(gamma, val, X, Y, lambda){
  this.lik = prof.lik.lasso(gamma, X, Y, lambda)
  to.root = this.lik - val
  return(to.root)
}

### Computes the penalized profile likelihood for gamma, with beta fit using CV lasso
### penal controls whether the l1 penalty should be added to the likelihood
pen.lik.CV.lasso = function(gamma, X, Y, folds = NULL, details = F,
                             lambda.type = "lambda.1se", all.lambdas=NULL){
  # browser()
  Z = BC(Y, gamma)
  fit = cv.glmnet(X, Z, foldid=folds, lambda = all.lambdas)
  l = fit[lambda.type]
  b = coef(fit, s=lambda.type)
  Z.hat = predict(fit, X, s = lambda.type)
  
  penalty = l1.pen(l, b)
  
  lik = pen.lik.formula(Y, Z, Z.hat, gamma, penalty)
  
  if(details){
    active = predict(fit, s=lambda.type, type="nonzero")
    output = list(lik, active, l, sum(abs(b)))
    return(output)
  }
  return(lik)
}

### Finding the root of this function solves 
### the equation penalized prof lik = val
pen.lik.root.CV.lasso = function(gamma, val, X, Y, folds = NULL, 
                             lambda.type = "lambda.1se",
                             all.lambdas = NULL){
  this.lik = pen.lik.CV.lasso(gamma, X, Y, folds = folds, 
                               lambda.type = lambda.type,
                               all.lambdas = all.lambdas)
  to.root = this.lik - val
  return(to.root)
}
