library(MASS)
library(glmnet)
library(pbapply)
library(tidyverse)
library(reshape2)
library(egg)
library(doParallel)


time = Sys.time()

set.seed(74799272)

source("LASSO_Likelihood_Helper_Functions.R")

K = 12 #Number of CV folds

n = 100
p = 100
q = 10
sigma = 1
gamma.0 = 2

#Smallest and largest gamma candidates
gamma.min = 1
gamma.max = 3
#Step size for gamma candidates
gamma.step = 0.1
Gammas = seq(gamma.min, gamma.max, gamma.step)
len.G = length(Gammas)

n.lambda.raw = 1000



### Generate data outside loop
X = matrix(rnorm(n * p, 10, 10), nrow = n, ncol = p)
beta = c(rep(1, times = q), rep(0, times = p - q))
mu.Y.raw = X %*% beta
mu.Y = mu.Y.raw + 10 * abs(min(mu.Y.raw))
Y = mu.Y + rnorm(n, 0, sigma)

Z = inv.BC(Y, gamma.0)



fit.cv.glmnet = cv.glmnet(X, Z, nfolds = K, keep = T,
                          nlambda = n.lambda.raw)
folds = fit.cv.glmnet$foldid
all.lambdas = fit.cv.glmnet$lambda
len.L = length(all.lambdas)


#Initialize parallelization
nclust = as.numeric(Sys.getenv("SLURM_NTASKS"))
nclust = ifelse(is.na(nclust), detectCores(), nclust)
cl = makeCluster(nclust)
registerDoParallel(cl)

#Pass info to cluster
clusterExport(cl, c("X", "folds", "Z", "all.lambdas",
                    "len.G", "len.L", "Gammas", "n"))
clusterEvalQ(cl, {
  library(glmnet)

  source("LASSO_Likelihood_Helper_Functions.R")
})


sim.output = pblapply(seq_len(K), function(i){
  set.seed(30795084 + 1000*i)

  X.train = X[folds != i,]
  Z.train = Z[folds != i]
  X.test = X[folds == i,]
  Z.test = Z[folds == i]

  all.gamma.hats = rep(0, times = len.L)
  all.errs = rep(0, times = len.L)

  for (j in seq_along(all.lambdas)) {
    # print(paste0(i, ":", j, " of ",
    #              K, ":", len.L))
    this.lambda = all.lambdas[j]

    this.likelihoods = rep(0, times = len.G)

    for (k in seq_along(Gammas)) {
      this.gamma = Gammas[k]
      Y.train = BC(Z.train, this.gamma) #Y.test is not used

      fit.cv = glmnet(X.train, Y.train)

      ### To Do ###
      #Choose gamma within the loop
      #Only get predictions for the 'best' gamma


      this.lik = get.profile.lik(this.gamma, X.train, Y.train,
                                 fit.cv, this.lambda)
      this.likelihoods[k] = this.lik
    }

    ind.gamma = which.min(this.likelihoods)
    gamma.hat = Gammas[ind.gamma]
    all.gamma.hats[j] = gamma.hat

    Y.train = BC(Z.train, gamma.hat)
    fit.cv = glmnet(X.train, Y.train)
    Y.hat = predict(fit.cv, X.test, s = this.lambda)
    Z.hat = inv.BC(Y.hat, gamma.hat)
    this.err = mean((Z.hat - Z.test)^2)
    all.errs[j] = this.err
  }

  output = list(gammas = all.gamma.hats, errs = all.errs)
  return(output)
}, cl=cl)


#Terminate parallelization
stopCluster(cl)

### Dimensions are:
### Fold number
### lambda
sim.errs = array(0, dim = c(K, len.L))
sim.gammas = array(0, dim = c(K, len.L))
for(i in 1:K){
  this.sim = sim.output[[i]]
  sim.errs[i,] = this.sim$errs
  sim.gammas[i,] = this.sim$gammas
}


#Get average CV errors and their SEs
errs.mean = apply(sim.errs, c(2), mean)
errs.se = apply(sim.errs, c(2), sd) / sqrt(K)

#Get lambda min and lambda 1se
ind.min = which.min(errs.mean)
err.min = errs.mean[ind.min]
se.min = errs.se[ind.min]
ind.1se = min(which(errs.mean < err.min + se.min))
lambda.min = all.lambdas[ind.min]
lambda.1se = all.lambdas[ind.1se]

#Prepare errors for plotting and construct plot
data.errs = data.frame(lambda = all.lambdas,
                       err = errs.mean,
                       upper = errs.mean + errs.se)

plot.errs = ggplot(data.errs, aes(x = lambda)) +
  geom_line(aes(y = err), size = 1.5) +
  geom_line(aes(y = upper), colour = "red")
plot(plot.errs)


#Prepare gamma estimates for plotting and construct plot
mean.gamma = apply(sim.gammas, 2, mean)
se.gamma = apply(sim.gammas, 2, sd)/sqrt(K)
data.gamma = data.frame(lambda = all.lambdas,
                        gamma = mean.gamma,
                        upper = mean.gamma + se.gamma,
                        lower = mean.gamma - se.gamma)
plot.gamma = ggplot(data.gamma, aes(x = lambda)) +
  geom_line(aes(y = gamma)) +
  geom_line(aes(y = upper), colour = "red") +
  geom_line(aes(y = lower), colour = "red") +
  ggtitle("Estimates of Gamma. Averaged Across CV Iterations.") +
  geom_rug()

plot(plot.gamma)



print(Sys.time() - time)