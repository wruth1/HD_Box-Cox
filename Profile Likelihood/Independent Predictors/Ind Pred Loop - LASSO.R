library(glmnet)
library(doParallel)
library(pbapply)
library(stringr)
library(dplyr)
library(ggplot2)



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
folds = sample(folds.raw)

all.pars = expand.grid(
  p = ps,
  sigma = sigmas,
  gamma.0 = gamma.0s,
  lambda.type.list = lambda.types
)


#Initialize parallelization
nclust = as.numeric(Sys.getenv("SLURM_NTASKS"))
nclust = ifelse(is.na(nclust), detectCores(), nclust)
nclust = 6
cl = makeCluster(nclust)
registerDoParallel(cl)

clusterSetRNGStream(cl=cl, iseed = 53567459)


#Pass info to cluster
clusterExport(cl, c("n", "Gammas",
                    "n.folds", "folds",
                    "all.pars"))
clusterEvalQ(cl, {
  library(glmnet)
  library(stringr)
  library(pbapply)
  library(dplyr)
  library(ggplot2)
  source("LASSO_Likelihood_Helper_Functions.R")
})

pblapply(seq_len(nrow(all.pars)), function(k){
  set.seed(32249631)
  print(paste0("Iteration - ", k))
  
  pars = all.pars[k, ]
  attach(pars)
  lambda.type = as.character(lambda.type.list[[1]])
  
  q = floor(sqrt(p))

  
  ### Generate data
  X = matrix(rnorm(n * p, 0, 1), nrow = n, ncol = p)
  beta = rep(0, times=p)
  mu.Y.raw = X %*% beta + 1
  mu.Y = mu.Y.raw + 10 * abs(min(mu.Y.raw))
  Y.lin = mu.Y + rnorm(n, 0, sigma)
  
  ### Transform Y so that gamma.0 is correct BC parameter
  Z = inv.BC(Y.lin, gamma.0)
  
  ### Find lambda values to use for all datasets
  ### Note: lambda is chosen as a proportion of max(beta.hat.ls)
  # X.std = scale(X)
  X.std = X
  
  
  
  #Run simulation
  # cat("A")
  source("Profile Likelihood/Make One Prof Lik - LASSO.R",
         local = T)
  if(is.null(sim.output)) return(NULL)
  
  #Create name for plots and files
  CV.type = str_extract(lambda.type, "\\w+$")
  plot.title = paste0("p = ", p, ", ",
                      "Sigma = ", sigma, ", ",
                      "gamma.0 = ", gamma.0, ", ",
                      "CV Type = ", CV.type)
  plot.title.full = paste0("Profile Likelihood for CV Lasso with ",
                           plot.title)
  
  # save.image(paste0("Profile Likelihood/Workspaces/LASSO/",
  #                   plot.title, ".RData"))
  
  #Make plots
  # cat("B")
  source("Profile Likelihood/Plot One CV Lambda.R",
         local = T)
  # cat("C")
  source("Profile Likelihood/Plot One Prof Lik - LASSO.R",
         local = T)
  
  # cat("\n")
  
  detach(pars)
  
}, cl=cl)

stopCluster(cl)

print(Sys.time() - time)