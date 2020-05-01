library(glmnet)
library(doParallel)
library(pbapply)
library(stringr)

set.seed(41460569)

time = Sys.time()

source("LASSO_Likelihood_Helper_Functions.R")


n = 100     #Sample Size

### Sizes of steps down from optimizer for profile likelihood CIs
step.sizes = 1:6

### Target number of Ys to generate for each X
### Actual value is slightly smaller to optimize parallelization
M.target = 6

#Smallest and largest gamma candidates
gamma.min = 1
gamma.max = 5
#Step size for gamma candidates
gamma.step = 0.2
#Candidate gamma values
Gammas = seq(gamma.min, gamma.max, gamma.step)
len.G = length(Gammas)

n.folds = 10

sigmas = c(1, 0.1)
gamma.0s = c(2, 3, 4)
lambda.types = c("lambda.1se", "lambda.min")
ps = c(10, 50)#, 200)
deltas = c(1, 10) # 2-norm of X %*% beta

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
nclust = 6
cl = makeCluster(nclust)
registerDoParallel(cl)

clusterSetRNGStream(cl=cl, iseed = 53567459)

### Number of (X, Y) pairs to generate at each parameter combination
### Note: Actual value is the greatest multiple 
### of nclust that is less than target
M = M.target - (M.target %% nclust)


# #Pass info to cluster
clusterExport(cl, c("all.pars", "M", "n", "gamma.min", "gamma.max",
                    "Gammas", "folds", "step.sizes", "gamma.step"))


clusterEvalQ(cl, {
  library(glmnet)
  library(pbapply)
  library(optimx)
  source("LASSO_Likelihood_Helper_Functions.R")
  source("Helper Functions/Profile Likelihood Helper Functions.R")
})

  
var.names = c("n", "p", "q", "sigma", "gamma.0", "delta",
              "M", "gamma.step", "step.sizes")
var.names = c(var.names,
              paste0("Minus ", step.sizes))
var.names = t(var.names)
write.table(var.names, "LASSO CIs/Coverages - LASSO.csv", 
            append = F, row.names = F, quote = F, sep = ",",
            col.names = F)


all.cover.probs = pbsapply(seq_len(nrow(all.pars)), function(j){
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
  source("LASSO CIs/Make One LASSO CI.R",
         local = T)
  
  detach(pars)
  
  # output = c(pars, cover.probs)
  # return(output)
  # return(j)
}, cl=cl)

stopCluster(cl)

print(Sys.time() - time)