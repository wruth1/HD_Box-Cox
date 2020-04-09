# ## Generate Y (i.e. generate epsilon)
# Y.lin = mu.Y + rnorm(n, 0, sigma)
# 
# ### Transform Y so that gamma.0 is correct BC parameter
# Z = inv.BC(Y.lin, gamma.0)


all.lambdas.raw = lapply(Gammas, function(this.gamma) {
  this.Z = BC(Z, this.gamma)
  this.fit = glmnet(X.std, this.Z)
  this.lambdas = this.fit$lambda
  return(this.lambdas)
})

all.lambdas = sort(unlist(all.lambdas.raw))


sim.output = lapply(Gammas, function(this.gamma){
  this.Z = BC(Z, this.gamma)
  sd.Z = sd(this.Z)
  this.fit = cv.glmnet(X.std, this.Z, lambda = all.lambdas,
                       nfolds = n.folds, foldid = folds)
  this.Z.hat = predict(this.fit, X.std, s=lambda.type)
  this.lambda = this.fit[lambda.type]

  #Profile likelihood
  prof.lik = profile.lik.formula(Z, this.Z, this.Z.hat, this.gamma)

  #Active set
  A.hat = predict(this.fit, s=lambda.type, type = "nonzero")

  output = list(prof.lik, A.hat)
  return(output)
})
