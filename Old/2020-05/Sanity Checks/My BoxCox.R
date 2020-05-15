library(MASS)
library(glmnet)
library(pbapply)
library(tidyverse)
library(reshape2)
library(egg)
library(doParallel)


time = Sys.time()

set.seed(88894930)

source("LASSO_Likelihood_Helper_Functions.R")

K = 10 #Number of CV folds

n = 100
p = 10
q = 2
sigma = 1
gamma.0 = 0.5
b = 1 #Size of each nonzero coefficient
mu.X = 100
sigma.X = 20

#Smallest and largest gamma candidates
gamma.min = 0
gamma.max = 1
#Step size for gamma candidates
gamma.step = 0.1
Gammas = seq(gamma.min, gamma.max, gamma.step)
len.G = length(Gammas)

n.lambda.raw = 1000



### Generate data outside loop
# X = matrix(rnorm(n * p, mu.X, sigma.X), nrow = n, ncol = p)
X = matrix(runif(n * p, 0, 10), nrow = n, ncol = p)
beta = c(rep(b, times = q), rep(0, times = p - q))
mu.Y.raw = X %*% beta
mu.Y = mu.Y.raw + 10 * abs(min(mu.Y.raw))
Y.lin = mu.Y + rnorm(n, 0, sigma)

Y = inv.BC(Y.lin, gamma.0)

fit.bc = boxcox(Y ~ X)
BC.Gammas = fit.bc$x

my.liks = sapply(BC.Gammas, function(this.gamma){
  Z = BC(Y, this.gamma)

  this.fit = lm(Z ~ X)
  Z.hat = predict(this.fit, data.frame(X))

  MSE = mean((Z - Z.hat)^2)
  lik = -n * log(MSE)/2
  Jacob = (this.gamma - 1)*sum(log(Y))
  lik = lik + Jacob
  return(lik)
})

plot(BC.Gammas, my.liks)


###################### Likelihoods aren't coming out right for ls!!!!




