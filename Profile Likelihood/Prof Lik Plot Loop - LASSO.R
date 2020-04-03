library(foreach)
library(MASS)
library(glmnet)
library(pbapply)
library(tidyverse)
library(doParallel)
library(glmpath)



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

n.lambda.raw = 100
n.folds = 10

sigmas = c(0.1, 0.01)
gamma.0s = c(2, 3, 4)
lambda.types = c("lambda.1se", "lambda.min")
ps = c(10, 50)

### Construct folds ###
fold.size = n %/% n.folds
n.leftover = n %% n.folds
folds.raw = rep(1:n.folds, times = fold.size)
leftovers = seq_len(n.leftover)
folds.raw = c(folds.raw, leftovers)
folds.select = sample(folds.raw)

fold.types = list(NULL, folds.select)

all.pars = expand.grid(
  p = ps,
  sigma = sigmas,
  gamma.0 = gamma.0s,
  lambda.type.list = lambda.types,
  fold.list = fold.types
)



foreach(j = seq_len(nrow(all.pars))) %do% {
  set.seed(32249631)
  
# foreach(j = c(1, 24)) %do% {
  cat(paste0(j, " of ", nrow(all.pars), ": "))
  pars = all.pars[j, ]
  attach(pars)
  if (is.null(fold.list)) {
    folds = fold.list
  } else{
    folds = fold.list[[1]]
  }
  lambda.type = as.character(lambda.type.list[[1]])
  
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
  cat("A")
  source("Profile Likelihood/Make One Prof Lik - LASSO.R")
  
  #Create name for plots and files
  CV.type = str_extract(lambda.type, "\\w+$")
  folds.type = ifelse(is.null(folds),
                      "Random", "Fixed")
  plot.title = paste0("p = ", p, ", ",
                      "Sigma = ", sigma, ", ",
                      "gamma.0 = ", gamma.0, ", ",
                      "CV Type = ", CV.type, ", ",
                      "Folds =", folds.type)
  plot.title.full = paste0("Profile Likelihood for CV Lasso with ",
                           plot.title)
  
  #Make plots
  cat("B")
  source("Profile Likelihood/Plot One CV Lambda.R")
  cat("C")
  source("Profile Likelihood/Plot One Prof Lik - LASSO.R")
  
  cat("\n")
  
  detach(pars)
  
}


print(Sys.time() - time)