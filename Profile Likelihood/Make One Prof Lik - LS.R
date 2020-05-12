sim.output = lapply(Gammas, function(this.gamma){
  this.Z = BC(Y, this.gamma)
  # sd.Z = sd(this.Z)
  this.fit = lm(this.Z ~ X.std)
  this.Z.hat = predict(this.fit)
  
  #Profile likelihood
  prof.lik = profile.lik.formula(Y, this.Z, this.Z.hat, this.gamma)
  
  return(prof.lik)
})
