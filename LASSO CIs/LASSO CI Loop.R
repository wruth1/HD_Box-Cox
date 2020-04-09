library(glmnet)
library(doParallel)
library(pbapply)
library(stringr)

set.seed(41460569)

time = Sys.time()

source("LASSO_Likelihood_Helper_Functions.R")


n = 100     #Sample Size
beta.size = 20

### Sizes of steps down from optimizer for profile likelihood CIs
step.sizes = 1:12

### Target number of Ys to generate for each X
### Actual value is slightly smaller to optimize parallelization
M.target = 100

#Smallest and largest gamma candidates
gamma.min = 1
gamma.max = 5
#Step size for gamma candidates
gamma.step = 0.1
#Candidate gamma values
Gammas = seq(gamma.min, gamma.max, gamma.step)
len.G = length(Gammas)

n.folds = 10

sigmas = c(0.1, 0.01)
gamma.0s = c(2, 3, 4)
lambda.types = c("lambda.1se", "lambda.min")
ps = c(10)#, 50, 200)
beta.types = c("norm", "val")

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
  beta.type = beta.types
)





#Initialize parallelization
nclust = as.numeric(Sys.getenv("SLURM_NTASKS"))
nclust = ifelse(is.na(nclust), detectCores(), nclust)
nclust = 6
cl = makeCluster(nclust)
registerDoParallel(cl)

clusterSetRNGStream(cl=cl, iseed = 53567459)

### Number of Ys to generate for each X
### Note: Actual value is the greatest multiple 
### of nclust that is less than target
M = M.target - (M.target %% nclust)


#Pass info to cluster
clusterExport(cl, c("n", "beta.size", "all.pars", "M",
                    "Gammas", "folds", "step.sizes"))
  # "mu.Y", "X.std", "lambda.type", "n.folds", "folds",
  #                   "gamma.0", "n", "Gammas", "sigma", "step.sizes"))
clusterEvalQ(cl, {
  library(glmnet)
  library(pbapply)
  source("LASSO_Likelihood_Helper_Functions.R")
})

  
all.cover.probs = pbsapply(seq_len(nrow(all.pars)), function(j){
# all.cover.probs = pbsapply(1:2, function(j){
  pars = all.pars[j,]
  attach(pars)
  
  lambda.type <<- as.character(lambda.type.list[[1]])
  
  ### Construct coefficient vector
  q = floor(sqrt(p))
  beta = make.beta(p, q, beta.size, type = beta.type)
  
  ### Generate data
  X = matrix(rnorm(n * p, 0, 1), nrow = n, ncol = p)
  mu.Y.raw = X %*% beta
  mu.Y = mu.Y.raw + 10 * abs(min(mu.Y.raw))
  Y.lin = mu.Y + rnorm(n, 0, sigma)
  
  ### Transform Y so that gamma.0 is correct BC parameter
  Z <<- inv.BC(Y.lin, gamma.0)
  
  ### Find lambda values to use for all datasets
  ### Note: lambda is chosen as a proportion of max(beta.hat.ls)
  # X.std = scale(X)
  X.std <<- X
  
  #Run simulation
  source("LASSO CIs/Make One LASSO CI.R",
         local = T)
  
  detach(pars)
  
  output = c(pars, cover.probs)
  return(output)
}, cl=cl)

stopCluster(cl)

print(Sys.time() - time)