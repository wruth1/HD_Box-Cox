##############################################################
### Examine profile likelihood for gamma at various lambda ###
##############################################################

library(ggplot2)
library(glmnet)

source("LASSO_Likelihood_Helper_Functions.R")

set.seed(74799272)


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


### Generate data inside loop
X = matrix(runif(n * p, 0, 10), nrow = n, ncol = p)
beta = c(rep(1, times = q), rep(0, times = p - q))
mu.Y.raw = X %*% beta
mu.Y = mu.Y.raw + 10 * abs(min(mu.Y.raw))
Y = mu.Y + rnorm(n, 0, sigma)

Z = inv.BC(Y, gamma.0)

fit.raw = glmnet(X, Z)
all.lambdas = fit.raw$lambda

all.gamma.hats = rep(0, times = length(all.lambdas))

for (j in seq_along(all.lambdas)) {
  # print(paste0(i, ":", j, " of ",
  #              K, ":", len.L))
  this.lambda = all.lambdas[j]
  
  this.likelihoods = rep(0, times = len.G)
  
  
  for (k in seq_along(Gammas)) {
    this.gamma = Gammas[k]
    this.Z = BC(Z, this.gamma) #Y.test is not used
    
    this.fit = glmnet(X, this.Z)
    
    
    
    this.lik = get.profile.lik(this.gamma, X, this.Z,
                               this.fit, this.lambda)
    this.likelihoods[k] = this.lik
  }
  
  ind.gamma = which.max(this.likelihoods)
  gamma.hat = Gammas[ind.gamma]
  all.gamma.hats[j] = gamma.hat
}


all.lambda.mins = rep(0, times = len.G)
all.lambda.1ses = rep(0, times = len.G)

for (k in seq_along(Gammas)) {
  this.gamma = Gammas[k]
  this.Z = BC(Z, this.gamma) #Y.test is not used
  
  this.fit = cv.glmnet(X, this.Z)
  
  this.lambda.min = this.fit$lambda.min
  all.lambda.mins[k] = this.lambda.min
  this.lambda.1se = this.fit$lambda.1se
  all.lambda.1ses[k] = this.lambda.1se
}


output = list(
  gammas = all.gamma.hats,
  lambdas = all.lambdas,
  lambda.mins = all.lambda.mins,
  lambda.1ses = all.lambda.1ses
)
