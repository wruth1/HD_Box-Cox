library(foreach)
library(MASS)
library(pbapply)
library(tidyverse)
library(doParallel)



time = Sys.time()

source("LASSO_Likelihood_Helper_Functions.R")


n = 100     #Sample Size
################################################## Uncomment these
# p = 200     #Number of predictors
# q = 10      #Number of active predictors
# sigma = 0.1   #Residual SD
# gamma.0 = 3 #Correct value for Box-Cox transformation
# beta.size = (1 + q) / 2 #Size of each non-zero coefficient
#Equal to average of 1:q

#Smallest and largest gamma candidates
gamma.min = 1
gamma.max = 5
#Step size for gamma candidates
gamma.step = 0.02
#Candidate gamma values
Gammas = seq(gamma.min, gamma.max, gamma.step)
len.G = length(Gammas)


sigmas = c(0.1, 0.01)
gamma.0s = c(2, 3, 4)
ps = c(10, 50)


all.pars = expand.grid(
  p = ps,
  sigma = sigmas,
  gamma.0 = gamma.0s
)



foreach(j = seq_len(nrow(all.pars))) %do% {
  set.seed(58322811)
  
# foreach(j = c(1, 10)) %do% {
  print(paste0(j, " of ", nrow(all.pars)))
  pars = all.pars[j, ]
  attach(pars)

  q = floor(sqrt(p))

  
  ### Generate data
  X = matrix(rnorm(n * p, 0, 1), nrow = n, ncol = p)
  beta = rep(0, times = p)
  mu.Y.raw = X %*% beta
  mu.Y = mu.Y.raw + 10 * abs(min(mu.Y.raw))
  Y.lin = mu.Y + rnorm(n, 0, sigma)
  
  ### Transform Y so that gamma.0 is correct BC parameter
  Z = inv.BC(Y.lin, gamma.0)
  
  ### Find lambda values to use for all datasets
  ### Note: lambda is chosen as a proportion of max(beta.hat.ls)
  # X.std = scale(X)
  X.std = X
  
  
  
  #Run simulation
  source("Profile Likelihood/Make One Prof Lik - LS.R")
  
  #Create name for plots and files
  plot.title = paste0("p = ", p, ", ",
                      "Sigma = ", sigma, ", ",
                      "gamma.0 = ", gamma.0)
  plot.title.full = paste0("Profile Likelihood for LS with ",
                           plot.title)
  
  #Make plots
  source("Profile Likelihood/Plot One Prof Lik - LS.R")
  
  
  detach(pars)
  
}


print(Sys.time() - time)