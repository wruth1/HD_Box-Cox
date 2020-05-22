library(glmnet)
library(doParallel)
library(pbapply)
library(stringr)
library(optimx)
library(ggplot2)
library(dplyr)
library(gtable)

set.seed(43655740)

time = Sys.time()
print(time)

source("Helper Functions/All Helper Scripts.R")


n = 100     #Sample Size
n.lambda = 200 # Target number of lambda candidates

lambda.type = "lambda.min"

### Amount to decrease likelihood by to construct CI
CI.step.size = qchisq(0.95, 1)/2

#Smallest and largest gamma candidates
gamma.min = -1
gamma.max = 3
#Step size for gamma candidates
gamma.step = 0.02
#Candidate gamma values
Gammas = seq(gamma.min, gamma.max, gamma.step)
len.G = length(Gammas)


sigmas = c(0.1, 1)
# sigmas = 1
gamma.0s = c(0, 1, 2)
# gamma.0s = 1
ps = c(10, 50, 200)
deltas = c(1, 3) #SD of X %*% beta
# deltas = 3
q.strs = c("sqrt", "full")
# q.strs = "sqrt"

#CV lambda strategy
CV.types = c("min", "1se")






### Construct folds ###
n.folds = 10
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
  q.str = q.strs,
  CV.type = CV.types
)


#Initialize parallelization
nclust = as.numeric(Sys.getenv("SLURM_NTASKS"))
nclust = ifelse(is.na(nclust), detectCores(), nclust)
cl = makeCluster(nclust)
registerDoParallel(cl)


#Pass info to cluster
clusterExport(cl, c("all.pars", "lambda.type", "n", "n.lambda", "folds",
                     "CI.step.size", "Gammas", "gamma.min", "gamma.max"))
clusterEvalQ(cl, {
  source("Helper Functions/All Helper Scripts.R")
  library(glmnet)
  library(ggplot2)
  library(dplyr)
  library(gtable)
})



pbsapply(seq_len(nrow(all.pars)), function(j){
# j = 5
  # print(paste0(j, " of ", nrow(all.pars)))
  set.seed(32249631)
  
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
  
  ### Find all candidate lambda values
  all.lambdas.raw = sapply(Gammas, function(gamma){
    this.Z = BC(Y, gamma)
    this.fit = glmnet(X, this.Z)
    this.lambdas = this.fit$lambda
    return(this.lambdas)
  })
  all.lambdas.fine = sort(unlist(all.lambdas.raw))
  
  ### Make the grid of lambda candidates coarser
  all.lambdas = coarsen.grid(n.lambda, all.lambdas.fine)
  
  
  
  ### Compute discretized profile penalized likelihood and active sets
  pen.analysis = lapply(Gammas, function(gamma){
    this.lik = pen.lik.CV.lasso(gamma=gamma, X=X, Y=Y, folds=folds,
                                 all.lambdas = all.lambdas, details=T,
                                lambda.type = lambda.type)
    return(this.lik)
  })
  
  ### Extract profile penalized likelihood and other details
  pen.prof.lik = sapply(pen.analysis, function(info) info[[1]])
  pen.active.sets = lapply(pen.analysis, function(info) info[[2]])
  pen.active.sizes = sapply(pen.analysis, function(info) {
    length(unlist(info[[2]]))
  })
  pen.lambda.hat = unlist(sapply(pen.analysis, function(info) info[[3]]))
  pen.b.norm = sapply(pen.analysis, function(info) info[[4]])
  
  ### Compute unpenalized profile likelihood
  prof.lik = sapply(Gammas, function(gamma){
    this.lik = prof.lik.CV.lasso(gamma=gamma, X=X, Y=Y, folds=folds,
                                all.lambdas = all.lambdas,
                                lambda.type = lambda.type)
    return(this.lik)
  })
  
  
  ### Store all simulation output in a data frame
  results = data.frame(lik.raw = prof.lik, lik = pen.prof.lik,
                       act = pen.active.sizes, gamma = Gammas,
                       l.hat = pen.lambda.hat, b.norm = pen.b.norm)
  
  #Create name for plots and files
  pars["q.str"] = as.character(q.str)
  pars["CV.type"] = as.character(CV.type)
  var.vals = paste0(names(pars),"=", pars)
  plot.title = paste0(var.vals, collapse = ", ")
  plot.title = paste0("j = ", j, ", ", plot.title)
  
  source("Penalized Likelihood/Plot One Pen Lik Decomp.R", local=T)
  
  detach(pars)
}, cl=cl)

stopCluster(cl)  

print(Sys.time() - time)