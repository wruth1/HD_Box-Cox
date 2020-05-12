library(glmnet)
library(doParallel)
library(pbapply)
library(stringr)
library(optimx)

set.seed(51798431)

time = Sys.time()

source("LASSO_Likelihood_Helper_Functions.R")


n = 100     #Sample Size

### Sizes of steps down from optimizer for profile likelihood CIs
step.sizes = 1:12

### Target number of Ys to generate for each X
### Actual value is slightly smaller to optimize parallelization
M.target = 50

#Smallest and largest gamma candidates
gamma.min = 1
gamma.max = 5
#Step size for gamma candidates
gamma.step = 0.01
#Candidate gamma values
Gammas = seq(gamma.min, gamma.max, gamma.step)
len.G = length(Gammas)


sigmas = c(0.1, 1)
gamma.0s = c(2, 3, 4)
ps = c(10, 50)#, 200)
deltas = c(1, 10) #SD of X %*% beta
q.strs = c("sqrt", "full")

#Construct list of parameter combinations
all.pars = expand.grid(
  p = ps,
  sigma = sigmas,
  gamma.0 = gamma.0s,
  delta = deltas,
  q.str = q.strs
)




# 
# #Initialize parallelization
# nclust = as.numeric(Sys.getenv("SLURM_NTASKS"))
# nclust = ifelse(is.na(nclust), detectCores(), nclust)
# nclust = 6
# cl = makeCluster(nclust)
# registerDoParallel(cl)
# 
# clusterSetRNGStream(cl=cl, iseed = 53567459)

### Number of Ys to generate for each X
### Note: Actual value is the greatest multiple 
### of nclust that is less than target
# M = M.target - (M.target %% nclust)
M = M.target


# #Pass info to cluster
# clusterExport(cl, c("n", "beta.size", "all.pars", "M",
#                     "Gammas", "folds", "step.sizes"))
# # "mu.Y", "X.std", "lambda.type", "n.folds", "folds",
# #                   "gamma.0", "n", "Gammas", "sigma", "step.sizes"))
# clusterEvalQ(cl, {
#   library(glmnet)
#   library(pbapply)
#   source("LASSO_Likelihood_Helper_Functions.R")
# })



all.cover.probs = pbsapply(seq_len(nrow(all.pars)), function(j){
  # print(j)
  # all.cover.probs = pbsapply(1:2, function(j){
  pars = all.pars[j,]
  attach(pars)
  
  
  ### Construct coefficient vector
  q = ifelse(q.str == "sqrt", sqrt(p), p)
  q = floor(q)
  # q = floor(sqrt(p))
  beta.size = delta / sqrt(q)
  beta = c(rep(beta.size, q),
           rep(0, p-q))
  
  
  
  ### Generate data
  X = matrix(rnorm(n * p, 0, 1), nrow = n, ncol = p)
  beta = c(rep(beta.size, times = q), rep(0, times = p - q))
  mu.Z.raw = X %*% beta
  mu.Z = mu.Z.raw + get.int(n, delta, sigma)
  Z = mu.Z + rnorm(n, 0, sigma)
  
  ### Transform Y so that gamma.0 is correct BC parameter
  Y = inv.BC(Z, gamma.0)
  
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
  
  # output = c(pars, cover.probs)
  # return(output)
})#, cl=cl)

# stopCluster(cl)

print(Sys.time() - time)












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
  beta.size = (1 + q) / 2 #Size of each non-zero coefficient
  
  
  ### Generate data
  X = matrix(rnorm(n * p, 0, 1), nrow = n, ncol = p)
  beta = c(rep(beta.size, times = q), rep(0, times = p - q))
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