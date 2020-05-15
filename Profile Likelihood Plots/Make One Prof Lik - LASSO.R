all.lambdas.raw = lapply(Gammas, function(this.gamma) {
  this.Z = BC(Y, this.gamma)
  this.fit = glmnet(X.std, this.Z)
  this.lambdas = this.fit$lambda
  return(this.lambdas)
})

all.lambdas = sort(unlist(all.lambdas.raw))

sim.output = lapply(Gammas, function(this.gamma){
  this.Z = BC(Y, this.gamma)
  this.fit = cv.glmnet(X.std, this.Z, lambda = all.lambdas, foldid = folds)
  this.Z.hat = predict(this.fit, X.std, s=lambda.type)
  this.lambda = this.fit[lambda.type]
  # this.data = data.frame(this.Z, X.std)
  # this.fit = lm(this.Z ~ ., data=this.data)
  # this.Z.hat = predict(this.fit, this.data)
  
  #Error on the transformed data scale
  err.BC = mean((this.Z - this.Z.hat)^2)
  
  #Error on the observed data scale
  Z.hat = inv.BC(this.Z.hat, this.gamma)
  err.original = mean((Z - Z.hat)^2)
  
  #Profile likelihood
  prof.lik = profile.lik.formula(Z, this.Z, this.Z.hat, this.gamma)

  #Jacobian
  Jacob = (this.gamma - 1) * sum(log(Z))
  
  #Active set
  A.hat = predict(this.fit, s=lambda.type, type = "nonzero")
  
  #Geometric mean of Z
  gm.Z = geom_mean(this.Z)
  
  #Return all computations
  this.err = data.frame(transformed = err.BC,
                        observed = err.original,
                        lik = prof.lik,
                        Jacob = Jacob,
                        lambda = this.lambda,
                        sd.Z = sd.Z,
                        gm.Z = gm.Z)
  output = list(this.err, A.hat)
  return(output)
})
