
library(MASS)
library(glmnet)
library(pbapply)
library(tidyverse)
library(doParallel)
library(glmpath)


time = Sys.time()

set.seed(32249631)

source("LASSO_Likelihood_Helper_Functions.R")

K = 12 #Number of datasets to simulate

n = 100     #Sample Size
p = 200     #Number of predictors
q = 10      #Number of active predictors
sigma = 0.1   #Residual SD
gamma.0 = 2 #Correct value for Box-Cox transformation
beta.size = (1+q)/2 #Size of each non-zero coefficient
                    #Equal to average of 1:q

#Smallest and largest gamma candidates
gamma.min = 1
gamma.max = 3
#Step size for gamma candidates
gamma.step = 0.01
#Candidate gamma values
Gammas = seq(gamma.min, gamma.max, gamma.step)
len.G = length(Gammas)

n.lambda.raw = 100



### Generate data
X = matrix(rnorm(n * p, 0, 1), nrow = n, ncol = p)
beta = c(rep(beta.size, times = q), rep(0, times = p - q))
mu.Y.raw = X %*% beta
mu.Y = mu.Y.raw + 10 * abs(min(mu.Y.raw))
Y.lin = mu.Y + rnorm(n, 0, sigma)

### Transform Y so that gamma.0 is correct BC parameter
Z = inv.BC(Y.lin, gamma.0)

### Find lambda values to use for all datasets
### Note: lambda is chosen as a proportion of max(beta.hat.ls)
# X.std = scale(X)
X.std = X


loss.LASSO = function(Y, Y.hat, b, lambda){
  L2 = sum((Y - Y.hat)^2)/(2*length(Y))
  L1 = lambda * sum(abs(b))
  return(L2 + L1)
}

lambda.type = "lambda.1se"

sim.output = pblapply(Gammas, function(this.gamma){
  this.Z = BC(Z, this.gamma)
  this.fit = cv.glmnet(X.std, this.Z)
  this.Z.hat = predict(this.fit, X.std, s=lambda.type)
  # this.data = data.frame(this.Z, X.std)
  # this.fit = lm(this.Z ~ ., data=this.data)
  # this.Z.hat = predict(this.fit, this.data)
  
  #Error on the transformed data scale
  err.BC = mean((this.Z - this.Z.hat)^2)
  
  #Error on the observed data scale
  Z.hat = inv.BC(this.Z.hat, this.gamma)
  err.original = mean((Z - Z.hat)^2)
  
  #Profile likelihood
  prof.lik = profile.lik.formula(Z, this.Z, this.Z.hat, this.gamma)

  #Jacobian
  Jacob = (this.gamma - 1) * sum(log(Z))
  
  #Active set
  A.hat = predict(this.fit, s=lambda.type, type = "nonzero")
  
  #Return all computations
  this.err = data.frame(transformed = err.BC,
                        observed = err.original,
                        lik = prof.lik,
                        Jacob = Jacob)
  output = list(this.err, A.hat)
  return(output)
})



#Separate losses and active set
all.losses.raw = sapply(seq_along(sim.output), function(i){
  return(sim.output[[i]][[1]])
})
all.A.hats = lapply(sim.output, function(out){
  return(out[[2]])
})


####################################
### Investigate different losses ###
####################################
all.losses.list = data.frame(t(all.losses.raw))
all.losses = apply(all.losses.list, 2, unlist)
results = data.frame(all.losses, gamma = Gammas)


plot.trans = ggplot(data = results, aes(x = Gammas, y = transformed)) +
  geom_line() + xlab("Gamma") + ylab("MSE") +
  ggtitle("MSE Computed on the BC-Transformed Data") +
  geom_vline(xintercept = gamma.0)
plot(plot.trans)

plot.obs = ggplot(data = results, aes(x = Gammas, y = observed)) +
  geom_point() + xlab("Gamma") + ylab("MSE") +
  ggtitle(paste0("MSE Computed on the Observed Data Scale", ". ",
          "I.e. After inverting the BC-Transformation")) +
  geom_vline(xintercept = gamma.0)
plot(plot.obs)

plot.lik = ggplot(data = results,
                  aes(x = Gammas, y = lik)) +
  geom_point() + xlab("Gamma") + ylab("Profile Likelihood") +
  # ggtitle(paste0("MSE Computed on the Observed Data Scale", ". ",
  #                "I.e. After inverting the BC-Transformation")) +
  geom_vline(xintercept = gamma.0)
plot(plot.lik)


#Decompose likelihood into MSE component and Jacobian
data.lik = results %>%
  mutate(MSE = -n*log(transformed)/2) %>%
  gather(key = "Type", value = "Score", MSE, Jacob, lik) %>%
  select(-observed, -transformed)
gamma.hat = Gammas[which.min(results$lik)]
plot.decomp = ggplot(data.lik,
                     aes(x = gamma, y = Score, color = Type)) +
  geom_line() + geom_vline(xintercept = gamma.0) +
  geom_vline(xintercept = gamma.hat, linetype=2)
plot(plot.decomp)


###############################
### Investigate active sets ###
###############################

all.lengths = sapply(all.A.hats, nrow)
print(all.lengths)

# all.same = sapply(all.A.hats, function(q){
#   all(q == w)
# })
# all(all.lengths)

plot.lik = ggplot(data = results,
                  aes(x = Gammas, y = lik)) +
  geom_point() + xlab("Gamma") + ylab("Profile Likelihood") +
  # ggtitle(paste0("MSE Computed on the Observed Data Scale", ". ",
  #                "I.e. After inverting the BC-Transformation")) +
  geom_vline(xintercept = gamma.0) 
plot(plot.lik)


##############################################
### Investigate banding in likelihood plot ###
##############################################
data.clust = select(results, gamma, lik)
dist.clust = dist(data.clust)
clust = hclust(dist.clust, method = "single")


m = 1  #Minimim number of clusters to plot
M = 10 #Maximum number of clusters to plot
all.groupings = cutree(clust, k = m:M)
for(j in rev(seq_len(M - m + 1))){
  this.groups = factor(all.groupings[,j])
  this.plot = ggplot(data = results,
                     aes(x = gamma, y = lik, 
                         colour = !!this.groups)) +
    geom_point()
  plot(this.plot)
}

