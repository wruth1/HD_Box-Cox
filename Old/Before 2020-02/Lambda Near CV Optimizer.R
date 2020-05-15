library(MASS)
library(glmnet)
library(pbapply)
library(tidyverse)
library(reshape2)
library(egg)

time = Sys.time()

set.seed(92794619)

source("LASSO_Likelihood_Helper_Functions.R")


n = 100
p = 100
q = 10
sigma = 1
gamma.0 = 1

#Smallest and largest gamma candidates
gamma.min = -3
gamma.max = 5
#Step size for gamma candidates
gamma.step = 0.05
Gammas = seq(gamma.min, gamma.max, gamma.step)
len.G = length(Gammas)


### Generate data outside loop
X = matrix(rnorm(n * p, 10, 10), nrow = n, ncol = p)
beta = c(rep(1, times = q), rep(0, times = p - q))
mu.Y.raw = X %*% beta
mu.Y = mu.Y.raw + 10 * abs(min(mu.Y.raw))
Y = mu.Y + rnorm(n, 0, sigma)
e = Y - mu.Y

Z = inv.BC(Y, gamma.0)


fit.cv.glmnet = cv.glmnet(X, Z, keep = T)
folds = fit.cv.glmnet$foldid

this.gamma = gamma.0
this.info = profile.lik.lasso(this.gamma, X, Z, folds)
this.lik = this.info[[1]]
this.lambdas = this.info[[2]]
this.vars = this.info[[3]]

all.lambdas = fit.cv.glmnet$lambda

all.beta.hats = coef(fit.cv.glmnet, all.lambdas)


coef.sums = apply(all.beta.hats, 1, function(estimates) {
  return(!all(estimates == 0))
})
ever.active = which(coef.sums)
ever.active = ever.active[-1] ### Ignore intercept estimates


### Find when the active set changes
changes = c()
for (i in seq_along(all.lambdas)[-1]) {
  j = i - 1
  this.set = which(all.beta.hats[, i] != 0)
  last.set = which(all.beta.hats[, j] != 0)
  additions = setdiff(this.set, last.set)
  deletions = setdiff(last.set, this.set)
  same.set = (length(additions) == 0) & (length(deletions) == 0)
  if (!same.set) {
    a = all.lambdas[i]
    b = all.lambdas[j]
    #changes = c(changes, mean(c(a, b)))
    changes = c(changes, a, b)
  }
}

for (i in seq_along(ever.active)) {
  this.ind = ever.active[i]
  this.beta.hats = all.beta.hats[this.ind, ]
  this.data = data.frame(l = all.lambdas,
                         b.hat = this.beta.hats)
  this.plot = ggplot(data = this.data,
                     aes(x = l, y = b.hat)) +
    geom_line(size = 1.5) + ggtitle(paste0("X", this.ind - 1)) +
    geom_vline(xintercept = changes, col = "red")
  plot(this.plot)
}
