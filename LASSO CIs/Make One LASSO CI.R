# prof.file = "RProfile.txt"
# Rprof(prof.file)

time = Sys.time()

library(glmnet)
library(doParallel)
library(pbapply)
library(stringr)

set.seed(10782897)


source("LASSO_Likelihood_Helper_Functions.R")


n = 100     #Sample Size
p = 200     #Number of predictors
q = floor(sqrt(p)) #Number of active predictors
sigma = 0.1   #Residual SD
gamma.0 = 4 #Correct value for Box-Cox transformation
# beta.size = (1 + q) / 2 #Size of each non-zero coefficient
beta.size = 20
beta.type = "norm"
beta = make.beta(p, q, beta.size, type = beta.type)

### Sizes of steps down from optimizer for profile likelihood CIs
step.sizes = 2*(1:12)

### Target number of Ys to generate for each X
### Actual value is slightly smaller to optimize parallelization
M.target = 100


### lambda min or lambda 1se?
lambda.type = "lambda.1se"

#Smallest and largest gamma candidates
gamma.min = gamma.0 - 1
gamma.max = gamma.0 + 1
#Step size for gamma candidates
gamma.step = 0.1
#Candidate gamma values
Gammas = seq(gamma.min, gamma.max, gamma.step)
len.G = length(Gammas)

n.folds = 10

### Construct folds ###
fold.size = n %/% n.folds
n.leftover = n %% n.folds
folds.raw = rep(1:n.folds, times = fold.size)
leftovers = seq_len(n.leftover)
folds.raw = c(folds.raw, leftovers)
folds = sample(folds.raw)


# cat(paste0(j, " of ", nrow(all.pars), ": "))

### Generate design matrix and compute mean of Y
X = matrix(rnorm(n * p, 0, 1), nrow = n, ncol = p)
mu.Y.raw = X %*% beta
mu.Y = mu.Y.raw + 10 * abs(min(mu.Y.raw))
mu.Y = unlist(mu.Y)
### Standardize columns of X
X.std = scale(X)


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
clusterExport(cl, c("n", "mu.Y", "X.std", "Gammas", "sigma",
                    "n.folds", "lambda.type", "folds", 
                    "step.sizes", "gamma.0"))
clusterEvalQ(cl, {
  library(glmnet)
  library(pbapply)
  source("LASSO_Likelihood_Helper_Functions.R")
})




all.coverages = pbsapply(seq_len(M), function(i) {
  ### Generate Y
  Y.lin = mu.Y + rnorm(n, 0, sigma)
  
  ### Transform Y so that gamma.0 is correct BC parameter
  Z <<- inv.BC(Y.lin, gamma.0)
  
  #Run simulation
  source("LASSO CIs/Make One Prof Lik.R",
         local = T)


  ### Extract profile likelihood and find maximizer
  prof.lik = sapply(sim.output, function(info)
    info[[1]])
  opt.lik = max(prof.lik)

  ### Compute CIs for several likelihood drop sizes
  threshs = opt.lik - step.sizes
  ints = lapply(threshs, function(l) {
    this.int.inds = which(prof.lik > l)
    this.int = Gammas[this.int.inds]
  })

  ### Check coverage
  cover = sapply(ints, function(int) {
    a = min(int)
    b = max(int)
    output = (a <= gamma.0) && (gamma.0 <= b)
  })
  
}, cl=cl)

stopCluster(cl)

### Compute coverage rates
cover.probs = apply(all.coverages, 1, mean)
names(cover.probs) = paste0("M", step.sizes)
print(paste0("Estimated coverage prob. over ", M, " datasets"))
print(cover.probs)


pars = c(n = n, p = p, q = q, sigma = sigma,
         gamma.0 = gamma.0, 
         beta.str = paste0(beta.type, " = ", beta.size),
         M = M, gamma.step = gamma.step)
results.raw = c(pars, signif(cover.probs, digits = 3))
results = data.frame(t(results.raw))

write.table(results, "LASSO CIs/Coverages.csv", append = T,
          row.names = F, quote = F, sep = ",",
          col.names = F)

# print(Sys.time() - time)

# Rprof(NULL)
# 
# summaryRprof(prof.file)
