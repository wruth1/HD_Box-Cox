library(foreach)
library(MASS)
library(glmnet)
library(pbapply)
library(tidyverse)
library(doParallel)
library(glmpath)



time = Sys.time()

source("LASSO_Likelihood_Helper_Functions.R")

K = 12 #Number of datasets to simulate

n = 100     #Sample Size
p = 200     #Number of predictors
q = 10      #Number of active predictors
################################################## Uncomment these
# sigma = 0.1   #Residual SD
# gamma.0 = 3 #Correct value for Box-Cox transformation
beta.size = (1 + q) / 2 #Size of each non-zero coefficient
#Equal to average of 1:q

#Smallest and largest gamma candidates
gamma.min = 1
gamma.max = 5
#Step size for gamma candidates
gamma.step = 0.1
#Candidate gamma values
Gammas = seq(gamma.min, gamma.max, gamma.step)
len.G = length(Gammas)

n.lambda.raw = 100
n.folds = 10

sigmas = c(0.1, 0.01)
gamma.0s = c(2, 3, 4)
lambda.types = c("lambda.1se", "lambda.min")

### Construct folds ###
fold.size = n %/% n.folds
n.leftover = n %% n.folds
folds.raw = rep(1:n.folds, times = fold.size)
leftovers = seq_len(n.leftover)
folds.raw = c(folds.raw, leftovers)
folds.select = sample(folds.raw)

fold.types = list(NULL, folds.select)


all.pars = expand.grid(
  sigma = sigmas,
  gamma.0 = gamma.0s,
  lambda.type.list = lambda.types,
  fold.list = fold.types
)



foreach(j = seq_len(nrow(all.pars))) %do% {
  set.seed(32249631)
  
# foreach(j = c(1, 24)) %do% {
  print(paste0(j, " of ", nrow(all.pars)))
  pars = all.pars[j, ]
  attach(pars)
  if (is.null(fold.list)) {
    folds = fold.list
  } else{
    folds = fold.list[[1]]
  }
  lambda.type = as.character(lambda.type.list[[1]])
  
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
  
  
  all.lambdas.raw = lapply(Gammas, function(this.gamma) {
    this.Z = BC(Z, this.gamma)
    this.fit = glmnet(X.std, this.Z)
    this.lambdas = this.fit$lambda
    return(this.lambdas)
  })
  
  all.lambdas = sort(unlist(all.lambdas.raw))
  
  #Run simulation
  source("Profile Likelihood/Make_One_Prof_Lik.R")
  
  #Create name for saved workspace
  CV.type = str_extract(lambda.type, "\\w+$")
  folds.type = ifelse(is.null(folds),
                      "Random", "Fixed")
  title = paste0("Sigma = ", sigma, ", ",
                      "gamma.0 = ", gamma.0, ", ",
                      "CV Type = ", CV.type, ", ",
                      "Folds =", folds.type)
  
  save.image(paste0("~/Output/", title, ".RData"))
  
  detach(pars)
  
}


print(Sys.time() - time)