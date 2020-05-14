library(glmnet)
library(doParallel)
library(pbapply)
library(stringr)
library(optimx)

set.seed(41457884)

time = Sys.time()

source("Helper Functions/All Helper Scripts.R")


n = 100     #Sample Size

### Sizes of steps down from optimizer for profile likelihood CIs
step.sizes = 2

### Number of times to repeat the entire process
### These repetitions allow assessment of uncertainty in coverage probs
B = 5

### Target number of Ys to generate for each X
### Actual value is slightly smaller to optimize parallelization
M.target = 200

#Smallest and largest gamma candidates
gamma.min = -2
gamma.max = 2
#Step size for gamma candidates
gamma.step = 0.01
#Candidate gamma values
Gammas = seq(gamma.min, gamma.max, gamma.step)
len.G = length(Gammas)


# sigmas = c(0.1, 1)
sigmas = 1
gamma.0s = c(0)
ps = c(10, 50)#, 200)
deltas = c(1, 3) #SD of X %*% beta
# q.strs = c("sqrt", "full")
q.strs = "sqrt"

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

var.names = c("n", "p", "q", "sigma", "gamma.0", "delta",
              "M", "gamma.step")
var.names = c(var.names,
              paste0("Minus ", step.sizes))
var.names = t(var.names)
write.table(var.names, "LASSO CIs/Coverages - LS.csv", append = F,
            row.names = F, quote = F, sep = ",",
            col.names = F)


pbsapply(seq_len(B), function(k) {
  sapply(seq_len(nrow(all.pars)), function(j) {
    # print(j)
    # all.cover.probs = pbsapply(1:2, function(j){
    pars = all.pars[j, ]
    attach(pars)
    
    
    ### Construct coefficient vector
    q = ifelse(q.str == "sqrt", sqrt(p), p)
    q = floor(q)
    beta.size = delta / sqrt(q)
    beta = c(rep(beta.size, q),
             rep(0, p - q))
    
    
    
    #Run simulation
    source("LASSO CIs/Make One LS CI.R",
           local = T)
    
    detach(pars)
    
    # output = c(pars, cover.probs)
    # return(output)
  })#, cl=cl)
  
})

# stopCluster(cl)

print(Sys.time() - time)