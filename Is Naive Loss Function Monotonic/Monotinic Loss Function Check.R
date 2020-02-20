
library(MASS)
library(glmnet)
library(pbapply)
library(tidyverse)
library(doParallel)
library(glmpath)


time = Sys.time()

set.seed(32249631)

source("LASSO_Likelihood_Helper_Functions.R")

K = 12 #Number of datasets to simulate

n = 100     #Sample Size
p = 200     #Number of predictors
q = 10      #Number of active predictors
sigma = 1   #Residual SD
gamma.0 = 2 #Correct value for Box-Cox transformation

#Smallest and largest gamma candidates
gamma.min = 1
gamma.max = 3
#Step size for gamma candidates
gamma.step = 0.1
#Candidate gamma values
Gammas = seq(gamma.min, gamma.max, gamma.step)
len.G = length(Gammas)

n.lambda.raw = 100



### Generate data
X = matrix(rnorm(n * p, 0, 1), nrow = n, ncol = p)
beta = c(rep(1, times = q), rep(0, times = p - q))
mu.Y.raw = X %*% beta
mu.Y = mu.Y.raw + 10 * abs(min(mu.Y.raw))
Y.lin = mu.Y + rnorm(n, 0, sigma)

### Transform Y so that gamma.0 is correct BC parameter
Z = inv.BC(Y.lin, gamma.0)

### Find lambda values to use for all datasets
### Note: lambda is chosen as a proportion of max(beta.hat.ls)
X.std = scale(X)


loss.LASSO = function(Y, Y.hat, b, lambda){
  L2 = sum((Y - Y.hat)^2)/(2*length(Y))
  L1 = lambda * sum(abs(b))
  return(L2 + L1)
}

lambda.type = "lambda.1se"

all.losses.raw = pblapply(Gammas, function(this.gamma){
  this.Z = BC(Z, this.gamma)
  this.fit = cv.glmnet(X.std, this.Z)
  this.Z.hat = predict(this.fit, X.std, s=lambda.type)
  
  #Error on the transformed data scale
  err.BC = mean((this.Z - this.Z.hat)^2)
  
  #Error on the observed data scale
  Z.hat = inv.BC(this.Z.hat, this.gamma)
  err.original = mean((Z - Z.hat)^2)
  
  #Return both errors
  this.err = data.frame(transformed = err.BC,
                        observed = err.original)
  return(this.err)
})

all.losses = data.frame(t(all.losses.raw))
results = data.frame(all.losses, gamma = Gammas)

plot.trans = ggplot(data = results, aes(x = Gammas,
                                        y = transformed)) +
  geom_line()
plot(plot.trans)


qplot(Gammas, log(all.losses$transformed), geom = "line")
