
#Initialize parallelization
nclust = as.numeric(Sys.getenv("SLURM_NTASKS"))
nclust = ifelse(is.na(nclust), detectCores(), nclust)
# nclust = 5
cl = makeCluster(nclust)
registerDoParallel(cl)

clusterSetRNGStream(cl=cl, iseed = 57857221)

#Pass info to cluster
clusterExport(cl, c("Z", "X.std", "all.lambdas", "lambda.type",
                    "n.folds", "folds"))
clusterEvalQ(cl, {
  library(glmnet)
  source("LASSO_Likelihood_Helper_Functions.R")
})

sim.output = pblapply(Gammas, function(this.gamma){
  this.Z = BC(Z, this.gamma)
  sd.Z = sd(this.Z)
  this.fit = cv.glmnet(X.std, this.Z, lambda = all.lambdas,
                       nfolds = n.folds, foldid = folds)
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
  
  #Return all computations
  this.err = data.frame(transformed = err.BC,
                        observed = err.original,
                        lik = prof.lik,
                        Jacob = Jacob,
                        lambda = this.lambda,
                        sd.Z = sd.Z)
  output = list(this.err, A.hat)
  return(output)
}, cl=cl)

#Terminate parallelization
stopCluster(cl)
