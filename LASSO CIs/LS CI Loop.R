library(glmnet)
library(doParallel)
library(pbapply)
library(stringr)

set.seed(51798431)

time = Sys.time()

source("LASSO_Likelihood_Helper_Functions.R")


n = 100     #Sample Size
beta.size = 20

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

n.folds = 10

sigmas = c(0.1, 0.01)
gamma.0s = c(2, 3, 4)
ps = c(10, 50)#, 200)
deltas = c(1, 10) #SD of a term in X %*% beta

#Construct list of parameter combinations
all.pars = expand.grid(
  p = ps,
  sigma = sigmas,
  gamma.0 = gamma.0s,
  delta = deltas
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

var.names = c("n", "p", "q", "sigma", "gamma.0", "beta.str",
              "M", "gamma.step")
var.names = c(var.names,
              paste0("Minus ", step.sizes))
var.names = t(var.names)
write.table(var.names, "LASSO CIs/Coverages - LS.csv", append = F,
            row.names = F, quote = F, sep = ",",
            col.names = F)

all.cover.probs = pbsapply(seq_len(nrow(all.pars)), function(j){
  # all.cover.probs = pbsapply(1:2, function(j){
  pars = all.pars[j,]
  attach(pars)
  

  ### Construct coefficient vector
  q = floor(sqrt(p))

  
  #Run simulation
  source("LASSO CIs/Make One LS CI.R",
         local = T)
  
  detach(pars)
  
  output = c(pars, cover.probs)
  return(output)
})#, cl=cl)

# stopCluster(cl)

print(Sys.time() - time)