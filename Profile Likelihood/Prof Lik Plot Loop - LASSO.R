library(glmnet)
library(doParallel)
library(pbapply)
library(stringr)
library(optimx)
library(ggplot2)

set.seed(51798431)

time = Sys.time()

source("Helper Functions/LASSO_Likelihood_Helper_Functions.R")
source("Helper Functions/Profile Likelihood Helper Functions.R")


n = 100     #Sample Size

### Amount to decrease likelihood by to construct CI
CI.step.size = qchisq(0.95, 1)/2

#Smallest and largest gamma candidates
gamma.min = -1
gamma.max = 3
#Step size for gamma candidates
gamma.step = 0.01
#Candidate gamma values
Gammas = seq(gamma.min, gamma.max, gamma.step)
len.G = length(Gammas)


sigmas = c(0.1, 1)
gamma.0s = c(0, 1, 2)
ps = c(10, 50)#, 200)
deltas = c(1, 3) #SD of X %*% beta
q.strs = c("sqrt", "full")




n.lambda.raw = 100
n.folds = 10


# lambda.types = c("lambda.1se", "lambda.min")
lambda.type = "lambda.1se"

### Construct folds ###
fold.size = n %/% n.folds
n.leftover = n %% n.folds
folds.raw = rep(1:n.folds, times = fold.size)
leftovers = seq_len(n.leftover)
folds.raw = c(folds.raw, leftovers)
folds = sample(folds.raw)


#Construct list of parameter combinations
all.pars = expand.grid(
  p = ps,
  sigma = sigmas,
  gamma.0 = gamma.0s,
  delta = deltas,
  q.str = q.strs
  # lambda.type.factor = lambda.types
)



pbsapply(seq_len(nrow(all.pars)), function(j){
  # print(paste0(j, " of ", nrow(all.pars)))
  set.seed(32249631)
  
  pars = all.pars[j,]
  attach(pars)
  # lambda.type.factor = as.character(lambda.type.fact)
  
  ### Construct coefficient vector
  q = ifelse(q.str == "sqrt", sqrt(p), p)
  q = floor(q)
  # q = floor(sqrt(p))
  beta.size = delta / sqrt(q)
  beta = c(rep(beta.size, q),
           rep(0, p-q))
  
  
  
  ### Generate data
  X = matrix(rnorm(n * p, 0, 1), nrow = n, ncol = p)
  mu.Z.raw = X %*% beta
  mu.Z = mu.Z.raw + get.int(n, delta, sigma)
  Z = mu.Z + rnorm(n, 0, sigma)
  
  ### Transform Y so that gamma.0 is correct BC parameter
  Y = inv.BC(Z, gamma.0)
  
  ### Find lambda values to use for all datasets
  ### Note: lambda is chosen as a proportion of max(beta.hat.ls)
  # X.std = scale(X)
  X.std = X
  
  
  ### Compute profile likelihood sequence
  prof.lik = sapply(Gammas, function(gamma){
    this.lik = prof.lik.CV.lasso(gamma=gamma, X=X, Y=Y, folds=folds)
    return(this.lik)
  })
  
  
  #Create name for plots and files
  pars[5] = as.character(q.str)
  var.vals = paste0(names(pars),"=", pars)
  plot.title = paste0(var.vals, collapse = ", ")
  plot.title = paste0("j = ", j, ", ", plot.title)
  
  source("Profile Likelihood/Plot One Prof Lik - LASSO.R", local=T)
  
  detach(pars)
  
  
  
})
  

print(Sys.time() - time)