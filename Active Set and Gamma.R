###############################################################
### Find the neighbourhood of ganna values around the truth ###
###    within which the fitted active set doesn't change    ###
###       Note: we are doing test-set tuning, not CV        ###
###############################################################


library(MASS)
library(glmnet)
library(pbapply)
library(tidyverse)
library(reshape2)
library(egg)
library(glmpath)


set.seed(58716566)

source("LASSO_Likelihood_Helper_Functions.R")


n = 100 
p = 200
q = 10
sigma = 3
gamma.0 = 0.5
beta = c(rep(1, times = q), rep(0, times = p - q))


#Smallest and largest gamma candidates
gamma.min = 0.499
gamma.max = 0.501
#Step size for gamma candidates
gamma.step = 0.00001
Gammas = seq(gamma.min, gamma.max, gamma.step)
len.G = length(Gammas)


### Generate data outside loop
X = matrix(rnorm(n * p, 10, 10), nrow = n, ncol = p)
mu.Y.raw = X %*% beta
mu.Y = mu.Y.raw + 10 * abs(min(mu.Y.raw))
Y.lin = mu.Y + rnorm(n, 0, sigma)
e = Y.lin - mu.Y

Y = inv.BC(Y.lin, gamma.0)

fit.cv = cv.glmnet(X, Y)
lambda.1se = fit.cv$lambda.1se
active.cv = predict(fit.cv, s=lambda.1se, type="nonzero")

matching = pbsapply(seq_along(Gammas), function(i){
  this.gamma = Gammas[i]
  Z = BC(Y, this.gamma)
  
  fit = glmnet(X, Z)
  active = predict(fit, s=lambda.1se, type="nonzero")
  return(isTRUE(all.equal(active, active.cv)))
})

data.matching = data.frame(gamma = Gammas, match = matching)

plot.matching = ggplot(data.matching, aes(x = gamma, y = match)) +
  geom_point() + geom_vline(xintercept = gamma.0, colour = "blue")
plot(plot.matching)
