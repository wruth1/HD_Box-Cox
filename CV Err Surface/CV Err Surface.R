library(glmnet)
library(doParallel)
library(pbapply)
library(stringr)
# library(optimx)
library(ggplot2)
library(dplyr)
library(gtable)
library(plot3D)
library(reshape2)


set.seed(72524312)

time = Sys.time()
print(time)

source("Helper Functions/All Helper Scripts.R")


n = 100     #Sample Size
n.lambda = 200 # Target number of lambda candidates

CV.type = "lambda.min"

### Amount to decrease likelihood by to construct CI
CI.step.size = qchisq(0.95, 1) / 2


#Step size for gamma candidates
gamma.grid.step = 0.01
#Maximum increase and decrease from gamma.0
gamma.step.up = 1
gamma.step.down = 1



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
# CV.types = c("min", "1se")






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
  q.str = q.strs
  # CV.type = CV.types
)

#
# #Initialize parallelization
# nclust = as.numeric(Sys.getenv("SLURM_NTASKS"))
# nclust = ifelse(is.na(nclust), detectCores(), nclust)
# cl = makeCluster(nclust)
# registerDoParallel(cl)
#
#
# #Pass info to cluster
# clusterExport(cl, c("all.pars", "lambda.type", "n", "n.lambda", "folds",
#                     "CI.step.size", "gamma.grid.step",
#                     "gamma.step.down", "gamma.step.up"))
# clusterEvalQ(cl, {
#   source("Helper Functions/All Helper Scripts.R")
#   library(glmnet)
#   library(ggplot2)
#   library(dplyr)
#   library(gtable)
# })



# pbsapply(seq_len(nrow(all.pars)), function(j){
j = 5
# print(paste0(j, " of ", nrow(all.pars)))
set.seed(24812765)

pars = all.pars[j,]
attach(pars)

#Smallest and largest gamma candidates
gamma.min = gamma.0 - gamma.step.down
gamma.max = gamma.0 + gamma.step.up
#Candidate gamma values
Gammas = seq(gamma.min, gamma.max, gamma.grid.step)

lambda.type = paste0("lambda.", CV.type)

### Construct coefficient vector
q = ifelse(q.str == "sqrt", sqrt(p), p)
q = floor(q)
# q = floor(sqrt(p))
beta.size = delta / sqrt(q)
beta = c(rep(beta.size, q),
         rep(0, p - q))



### Generate data
X = matrix(rnorm(n * p, 0, 1), nrow = n, ncol = p)
mu.Z.raw = X %*% beta
mu.Z = mu.Z.raw + get.int(n, delta, sigma)
Z = mu.Z + rnorm(n, 0, sigma)

### Transform Z so that gamma.0 is correct BC parameter
Y = inv.BC(Z, gamma.0)

### Find all candidate lambda values
all.lambdas.raw = sapply(Gammas, function(gamma) {
  this.Z = BC(Y, gamma)
  this.fit = glmnet(X, this.Z)
  this.lambdas = this.fit$lambda
  return(this.lambdas)
})
all.lambdas.fine = sort(unlist(all.lambdas.raw))

### Make the grid of lambda candidates coarser
all.lambdas = coarsen.grid(n.lambda, all.lambdas.fine)


all.CV.errors.raw = pbsapply(Gammas, function(gamma) {
  this.Z = BC(Y, gamma)
  this.fit = cv.glmnet(X, this.Z, lambda = all.lambdas,
                       foldid = folds)
  this.errs = this.fit$cvm
  return(this.errs)
})


### Process array of CV errors
all.CV.errors = melt(all.CV.errors.raw)
colnames(all.CV.errors) = c("lambda", "gamma", "CV_Error")
# all.CV.errors = log(all.CV.errors.raw)
# colnames(all.CV.errors) = c("x", "y", "z")
# info = list(x = all.CV.errors$lambda,
#             y = all.CV.errors$gamma,
#             z = all.CV.errors$CV_Err)
all.CV.errors = all.CV.errors %>% mutate(log_err = log(CV_Error))

plot_ly(
  all.CV.errors,
  x = ~ lambda,
  y = ~ gamma,
  z = ~ log_err
) %>%
  add_surface()

plot_ly(x = Gammas, y = all.lambdas, z = log(all.CV.errors.raw)) %>%
  add_surface() %>%
  layout(scene = list(
    xaxis = list(title = "gamma"),
    yaxis = list(title = "lambda"),
    zaxis = list(title = "CV-Error")
  ))

### Plot the CV errors
x11()
persp3D(
  x = all.lambdas,
  y = Gammas,
  z = all.CV.errors,
  xlab = "lambda",
  ylab = "gamma",
  zlab = "CV Error"
)
dev.off()

kd <- with(MASS::geyser, MASS::kde2d(duration, waiting, n = 50))
fig <- plot_ly(x = kd$x, y = kd$y, z = kd$z) %>% add_surface()

detach(pars)
# }, cl=cl)

# stopCluster(cl)

print(Sys.time() - time)
