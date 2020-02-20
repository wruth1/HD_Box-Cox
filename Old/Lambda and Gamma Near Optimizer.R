####################################################################
### Find the neighbourhood of lambda values around the CV        ###
### optimizer within which the fitted active set doesn't change  ###
### I.e. The locations of the two closest knots to lambda.hat.CV ###
####################################################################


library(MASS)
library(glmnet)
library(pbapply)
library(tidyverse)
library(reshape2)
library(egg)
library(glmpath)

time = Sys.time()

set.seed(74380493)

source("LASSO_Likelihood_Helper_Functions.R")


n = 100
p = 100
q = 10
sigma = 1
gamma.0 = 2
beta = c(rep(1, times = q), rep(0, times = p - q))


#Smallest and largest gamma candidates
gamma.min = gamma.0 - 1
gamma.max = gamma.0 + 1
#Step size for gamma candidates
gamma.step = 0.01
Gammas = seq(gamma.min, gamma.max, gamma.step)
len.G = length(Gammas)


### Generate data outside loop
X = matrix(rnorm(n * p, 10, 10), nrow = n, ncol = p)
mu.Y.raw = X %*% beta
mu.Y = mu.Y.raw + 10 * abs(min(mu.Y.raw))
Y.lin = mu.Y + rnorm(n, 0, sigma)
e = Y.lin - mu.Y

fit.lin = cv.glmnet(X, Y.lin, keep=T)
folds = fit.lin$foldid

Y = inv.BC(Y.lin, gamma.0)

### Start with lambda.min

info.lambda.min = pbsapply(seq_along(Gammas), function(j){

  this.gamma = Gammas[j]
  Z = BC(Y, this.gamma)
  
  fit.cv.glmnet = cv.glmnet(X, Z, foldid = folds)
  lambda.min = fit.cv.glmnet$lambda.min


  path = glmpath(X, Z, family = "gaussian", max.norm = 10000*ncol(X),
                 min.lambda = 0, lambda2=0)
  lambda.seq = path$lambda
  ind.lambda.min = which(diff(lambda.min > lambda.seq) != 0)
  int.lambda.min = lambda.seq[c(ind.lambda.min, ind.lambda.min + 1)]
  lower.lambda.min = int.lambda.min[2]
  upper.lambda.min = int.lambda.min[1]
  
  output = c(lambda.min, lower.lambda.min, upper.lambda.min)
  
  # all.lambda.mins[j] = lambda.min
  # all.lower.lambda.mins[j] = lower.lambda.min
  # all.upper.lambda.mins[j] = upper.lambda.min
  
  return(output)
})

ranges.lambda.min = data.frame(gamma = Gammas, 
                               lambda = log(info.lambda.min[1,]),
                               lower = log(info.lambda.min[2,]),
                               upper = log(info.lambda.min[3,]))


plot.lambda.min = ggplot(data=ranges.lambda.min, aes(x = gamma)) +
  geom_line(aes(y = lambda), size = 1.5) +
  geom_line(aes(y = lower), colour = "red") +
  geom_line(aes(y = upper), colour = "red") +
  geom_vline(xintercept = gamma.0, colour = "blue") +
  ylab("log(lambda)") + ggtitle("lambda-min")
plot(plot.lambda.min)

ranges.lambda.min = mutate(ranges.lambda.min,
                           diff = upper - lower)

plot.diff.min = ggplot(data = ranges.lambda.min, aes(x = gamma)) +
  geom_line(aes(y = diff)) + ggtitle("lambda - min") +
  ylab("Length of Interval (log-scale)") +
  geom_vline(xintercept = gamma.0, colour = "blue")
plot(plot.diff.min)

### Repeat the process for lambda.1se
info.lambda.1se = pbsapply(seq_along(Gammas), function(j){

  this.gamma = Gammas[j]
  Z = BC(Y, this.gamma)

  fit.cv.glmnet = cv.glmnet(X, Z, foldid = folds)
  folds = fit.cv.glmnet$foldid
  lambda.1se = fit.cv.glmnet$lambda.1se

  path = glmpath(X, Z, family = "gaussian", max.norm = 10000*ncol(X),
                 min.lambda = 0)
  lambda.seq = path$lambda
  ind.lambda.1se = which(diff(lambda.1se > lambda.seq) != 0)
  int.lambda.1se = lambda.seq[c(ind.lambda.1se, ind.lambda.1se + 1)]
  lower.lambda.1se = int.lambda.1se[2]
  upper.lambda.1se = int.lambda.1se[1]

  output = c(lambda.1se, lower.lambda.1se, upper.lambda.1se)


  return(output)
})

ranges.lambda.1se = data.frame(gamma = Gammas,
                               lambda = log(info.lambda.1se[1,]),
                               lower = log(info.lambda.1se[2,]),
                               upper = log(info.lambda.1se[3,]))

plot.lambda.1se = ggplot(data=ranges.lambda.1se, aes(x = gamma)) +
  geom_line(aes(y = lambda), size = 1.5) +
  geom_line(aes(y = lower), colour = "red") +
  geom_line(aes(y = upper), colour = "red") +
  geom_vline(xintercept = gamma.0, colour = "blue") +
  ylab("log(lambda)") + ggtitle("lambda-1se")
plot(plot.lambda.1se)

ranges.lambda.1se = mutate(ranges.lambda.1se,
                           diff = upper - lower)

plot.diff.1se = ggplot(data = ranges.lambda.1se, aes(x = gamma)) +
  geom_line(aes(y = diff)) + ggtitle("lambda - 1se") +
  ylab("Length of Interval (log-scale)") +
  geom_vline(xintercept = gamma.0, colour = "blue")
plot(plot.diff.1se)
