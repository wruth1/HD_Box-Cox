sim.output = lapply(Gammas, function(this.gamma){
  this.Z = BC(Z, this.gamma)
  # sd.Z = sd(this.Z)
  this.fit = lm(this.Z ~ X.std)
  this.Z.hat = predict(this.fit)

  #Profile likelihood
  prof.lik = profile.lik.formula(Z, this.Z, this.Z.hat, this.gamma)
  
  output = list(prof.lik)
  return(output)
})
