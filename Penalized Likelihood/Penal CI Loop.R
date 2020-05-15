library(glmnet)
library(doParallel)
library(pbapply)
library(stringr)

set.seed(41460569)

time = Sys.time()
print(time)

source("Helper Functions/All Helper Scripts.R")


n = 100     #Sample Size
n.lambda = 100 #Minimum number of candidate lambdas

### Sizes of steps down from optimizer for profile likelihood CIs
step.sizes = 2

### Number of times to repeat the entire process
### Allows for assessment of uncertainty in coverage probs
B = 5

### Number of (X,Y) pairs to generate for each parameter setting
M = 50

#Smallest and largest gamma candidates
gamma.min = -1
gamma.max = 3
#Step size for gamma candidates
gamma.step = 0.05
#Candidate gamma values
Gammas = seq(gamma.min, gamma.max, gamma.step)
len.G = length(Gammas)

n.folds = 10

sigmas = c(1)
gamma.0s = c(0)
lambda.types = c("lambda.1se")
ps = c(10, 50, 200)
deltas = c(1, 3) # 2-norm of X %*% beta

### Construct folds ###
fold.size = n %/% n.folds
n.leftover = n %% n.folds
folds.raw = rep(1:n.folds, times = fold.size)
leftovers = seq_len(n.leftover)
folds.raw = c(folds.raw, leftovers)
folds = sample(folds.raw)

all.pars = expand.grid(
  p = ps,
  sigma = sigmas,
  gamma.0 = gamma.0s,
  lambda.type.list = lambda.types,
  delta = deltas)

# 
# #Initialize parallelization
# nclust = as.numeric(Sys.getenv("SLURM_NTASKS"))
# nclust = ifelse(is.na(nclust), detectCores(), nclust)
nclust = min(c(6, B))
cl = makeCluster(nclust)
registerDoParallel(cl)

clusterSetRNGStream(cl=cl, iseed = 43655740)

### Number of (X, Y) pairs to generate at each parameter combination
### Note: Actual value is the greatest multiple 
### of nclust that is less than target


# #Pass info to cluster
clusterExport(cl, c("all.pars", "M", "n", "gamma.min", "gamma.max",
                    "Gammas", "folds", "step.sizes", "gamma.step", 
                    "n.lambda"))


clusterEvalQ(cl, {
  library(glmnet)
  library(pbapply)
  library(optimx)
  source("Helper Functions/All Helper Scripts.R")
})

  
var.names = c("n", "p", "q", "sigma", "gamma.0", "delta",
              "M", "gamma.step", "step.sizes")
var.names = c(var.names,
              paste0("Minus ", step.sizes))
var.names = t(var.names)
write.table(var.names, "Penalized Likelihood/Coverages - LASSO.csv", 
            append = F, row.names = F, quote = F, sep = ",",
            col.names = F)


parSapply(X = seq_len(B), FUN = function(k){
sapply(X = seq_len(nrow(all.pars)), FUN = function(j){
  print(paste0(j, " of ", nrow(all.pars)))
# all.cover.probs = pbsapply(1:2, function(j){
  pars = all.pars[j,]
  attach(pars)
  
  lambda.type = as.character(lambda.type.list[[1]])
  
  ### Construct coefficient vector
  q = floor(sqrt(p))
  beta.size = delta / sqrt(q) # Value of each non-zero term in beta
  beta = c(rep(beta.size, times = q),
           rep(0, times = p-q))
  
  #Run simulation
  source("Penalized Likelihood/Make One Penal CI.R",
         local = T)
  
  detach(pars)
  
  # output = c(pars, cover.probs)
  # return(output)
  # return(j)
})

}, cl=cl)


stopCluster(cl)

print(Sys.time() - time)