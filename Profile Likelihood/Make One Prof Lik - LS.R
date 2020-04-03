# 
# #Initialize parallelization
# nclust = as.numeric(Sys.getenv("SLURM_NTASKS"))
# nclust = ifelse(is.na(nclust), detectCores(), nclust)
# nclust = 5
# cl = makeCluster(nclust)
# registerDoParallel(cl)
# 
# clusterSetRNGStream(cl=cl, iseed = 57857221)
# 
# #Pass info to cluster
# clusterExport(cl, c("Z", "X.std"))
# clusterEvalQ(cl, {
#   source("LASSO_Likelihood_Helper_Functions.R")
# })

sim.output = pblapply(Gammas, function(this.gamma){
  this.Z = BC(Z, this.gamma)
  sd.Z = sd(this.Z)
  this.fit = lm(this.Z ~ X.std)
  this.Z.hat = predict(this.fit, data.frame(X.std))
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
  
  
  #Return all computations
  this.err = data.frame(transformed = err.BC,
                        observed = err.original,
                        lik = prof.lik,
                        Jacob = Jacob,
                        sd.Z = sd.Z)
  output = list(this.err)
})

# stopCluster(cl)