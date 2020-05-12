###########################################
### Box-Cox Analysis of Riboflavin Data ###
###########################################

library(hdi)
library(MASS)
library(glmnet)
library(pbapply)
library(ggplot2)
library(doParallel)

source("Helper Functions/Profile Likelihood Helper Functions.R")
source("Helper Functions/LASSO_Likelihood_Helper_Functions.R")

set.seed(3885486)

CI.step.size = qchisq(0.95, 1)/2


### Construct candidate gamma values
#Smallest and largest gamma candidates
gamma.min = -2
gamma.max = 2
#Step size for gamma candidates
gamma.step = 0.01
#Candidate gamma values
Gammas = seq(gamma.min, gamma.max, gamma.step)
len.G = length(Gammas)

### Import data
data("riboflavin")
data = riboflavin
# dim(data)
# str(data[])
X = data$x
Y = 2^data$y


############################################
###         Primitive Analysis           ###
### Variable selection followed by LS BC ###
############################################


fit = cv.glmnet(X, Y)
pred = predict(fit, X, s = "lambda.1se")
resid = Y - pred

# plot(pred, resid)
# qqnorm(resid)
# qqline(resid)


active = unlist(predict(fit, s = "lambda.1se", type = "nonzero"))


X.hat = X[,active]

fit.ls = lm(Y ~ X.hat)
# plot(fit.ls)

# w = boxcox(fit.ls, lambda = Gammas)
# liks = w$y
# ls.gammas = w$x
# max.lik = max(liks)
# inds.int = which(liks >= max.lik - CI.step.size)
# conf.set = ls.gammas[inds.int]
# 
# CI = c(min(conf.set), max(conf.set))

################################################
###        Naive likelihood analysis         ###
### LASSO residuals in LS profile likelihood ###
################################################

n.folds = 10
# Target length for all.lambdas
len.lambda = 200

### Construct folds
n = nrow(data)
fold.size = n %/% n.folds
n.leftover = n %% n.folds
folds.raw = rep(1:n.folds, times = fold.size)
leftovers = seq_len(n.leftover)
folds.raw = c(folds.raw, leftovers)
folds = sample(folds.raw)


#Initialize parallelization
nclust = as.numeric(Sys.getenv("SLURM_NTASKS"))
nclust = ifelse(is.na(nclust), detectCores(logical=F), nclust)
nclust = nclust - 1 ### Save one core for computing while code runs
# nclust = 6
cl = makeCluster(nclust)
registerDoParallel(cl)

clusterSetRNGStream(cl=cl, iseed = 53567459)


# #Pass info to cluster
clusterExport(cl, c("X", "Y", "folds"))
clusterEvalQ(cl, {
  library(glmnet)
  source("Helper Functions/Profile Likelihood Helper Functions.R")
  source("Helper Functions/LASSO_Likelihood_Helper_Functions.R")
})

### Extract all candidate lambda values
all.lambdas.raw = pbsapply(Gammas, function(gamma){
  this.Z = BC(Y, gamma)
  this.fit = glmnet(X, this.Z)
  this.lambdas = this.fit$lambda
  return(this.lambdas)
})
all.lambdas.raw = sort(all.lambdas.raw)

### Coarsify candidate lambda values
coarseness = floor(log2(length(all.lambdas.raw)/len.lambda)) - 3
all.lambdas = all.lambdas.raw
for(j in 1:coarseness){
  all.lambdas = coarser.grid(all.lambdas)
}


clusterExport(cl, "all.lambdas")

prof.lik = pbsapply(Gammas, function(gamma){
  this.lik = prof.lik.CV.lasso(gamma=gamma, X=X, Y=Y, folds=folds,
                               all.lambdas = all.lambdas)
  return(this.lik)
}, cl=cl)

stopCluster(cl)

### Plot LASSO profile likelihood
data.lik = data.frame(gamma = Gammas, lik = prof.lik)
lik.plot = ggplot(data = data.lik, mapping = aes(x=gamma, y=lik)) +
  # geom_point()
  geom_line()
# plot(lik.plot)


### Maximize LASSO profile likelihood and construct CI
max.lik = max(prof.lik)
lasso.mle = Gammas[which.max(prof.lik)]
inds.int = which(prof.lik >= max.lik - CI.step.size)
CI = Gammas[inds.int]
CI.bounds = c(min(CI), max(CI))


### Add CI info to LASSO profile likelihood plot
lik.CI.plot = lik.plot + geom_vline(xintercept = lasso.mle) +
  geom_hline(yintercept = max.lik - CI.step.size, lty = 2) +
  geom_vline(xintercept = CI.bounds, lty = 2)
plot(lik.CI.plot)
