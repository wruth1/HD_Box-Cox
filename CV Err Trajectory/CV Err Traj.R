library(glmnet)
library(doParallel)
library(pbapply)
library(stringr)
library(optimx)
library(ggplot2)
library(dplyr)
library(lars)

set.seed(37280337)

time = Sys.time()

source("Helper Functions/All Helper Scripts.R")


n = 100        #Sample Size
n.lambda = 400 # Target number of lambda candidates


#Smallest and largest gamma candidates
gamma.min = -1
gamma.max = 3
#Step size for gamma candidates
gamma.step = 0.01
#Candidate gamma values
Gammas = seq(gamma.min, gamma.max, gamma.step)
len.G = length(Gammas)


sigmas = c(0.1, 1)
# sigmas = 1
gamma.0s = c(1, 2)
# gamma.0s = 1
gamma.hats = c(1, 2)
ps = c(10, 50, 200)
# ps = 50
deltas = c(1, 3) #SD of X %*% beta
# deltas = 3
q.strs = c("sqrt", "full")
# q.strs = "sqrt"

#CV lambda strategy
CV.types = c("min", "1se")


n.folds = 10


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
  gamma.hat = gamma.hats,
  delta = deltas,
  q.str = q.strs,
  CV.type = CV.types
)



pbsapply(seq_len(nrow(all.pars)), function(j){
  # j = 5
  # print(paste0(j, " of ", nrow(all.pars)))
  set.seed(22949690)
  
  pars = all.pars[j,]
  attach(pars)

  lambda.type = paste0("lambda.", CV.type)
  
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
  
  
  # ### Find all candidate lambda values
  # all.lambdas.raw = sapply(Gammas, function(gamma){
  #   this.Z = BC(Y, gamma)
  #   this.fit = glmnet(X, this.Z)
  #   this.lambdas = this.fit$lambda
  #   return(this.lambdas)
  # })
  # all.lambdas.fine = sort(unlist(all.lambdas.raw))
  # 
  # ### Make the grid of lambda candidates coarser
  # all.lambdas = coarsen.grid(n.lambda, all.lambdas.fine)
  # 
  # 
  # 
  # ### Compute CV error for each lambda value
  # all.CV.errs = sapply(all.lambdas, function(l){
  #   
  # })
  # 
  # 
  # 
  # ### Compute profile likelihood sequence
  # prof.lik = pbsapply(Gammas, function(gamma){
  #   this.lik = pen.lik.CV.lasso(gamma=gamma, X=X, Y=Y, folds=folds,
  #                               all.lambdas = all.lambdas)
  #   return(this.lik)
  # })
  
  
  #Create name for plots and files
  pars["q.str"] = as.character(q.str)
  pars["CV.type"] = as.character(CV.type)
  var.vals = paste0(names(pars),"=", pars)
  plot.title = paste0(var.vals, collapse = ", ")
  plot.title = paste0("j = ", j, ", ", plot.title)
  
  file.name = paste0("Plots/CV Lambda Traj/", CV.type, "/", plot.title, ".jpg")
  jpeg(file.name)
  cv.lars(X, Y, se=F)
  dev.off()
  
  
  
  detach(pars)
  
  
  
})


print(Sys.time() - time)