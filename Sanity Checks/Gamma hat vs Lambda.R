#########################################################
### Which gamma values are chosen at each lambda, and ###
###    which lambda values are chosen at each gamma   ###
###      Also, find gamma hat at CV lambda values     ###
#########################################################


library(MASS)
library(glmnet)
library(pbapply)
library(tidyverse)
library(reshape2)
library(egg)
library(doParallel)


time = Sys.time()

set.seed(74799272)

source("LASSO_Likelihood_Helper_Functions.R")

K = 12 #Number of CV folds

n = 100
p = 100
q = 10
sigma = 1
gamma.0 = 2

#Smallest and largest gamma candidates
gamma.min = 1
gamma.max = 5
#Step size for gamma candidates
gamma.step = 0.1
Gammas = seq(gamma.min, gamma.max, gamma.step)
len.G = length(Gammas)

n.lambda.raw = 1000


################################################
### Choose gamma at various values of lambda ###
###        across multiple datasets          ###
################################################

#Initialize parallelization
nclust = as.numeric(Sys.getenv("SLURM_NTASKS"))
nclust = ifelse(is.na(nclust), detectCores(), nclust)
cl = makeCluster(nclust)
registerDoParallel(cl)

#Pass info to cluster
clusterExport(cl, c("n", "p", "q", "sigma", "gamma.0", "len.G",
                    "Gammas"))

clusterEvalQ(cl, {
  library(glmnet)
  
  source("LASSO_Likelihood_Helper_Functions.R")
})


sim.output = pblapply(seq_len(K), function(i) {
  set.seed(30795084 + 10000 * i)
  
  ### Generate data inside loop
  X = matrix(runif(n * p, 0, 10), nrow = n, ncol = p)
  beta = c(rep(1, times = q), rep(0, times = p - q))
  mu.Y.raw = X %*% beta
  mu.Y = mu.Y.raw + 10 * abs(min(mu.Y.raw))
  Y = mu.Y + rnorm(n, 0, sigma)
  
  Z = inv.BC(Y, gamma.0)
  
  fit.raw = glmnet(X, Z)
  all.lambdas = fit.raw$lambda
  
  all.gamma.hats = rep(0, times = length(all.lambdas))
  
  for (j in seq_along(all.lambdas)) {
    # print(paste0(i, ":", j, " of ",
    #              K, ":", len.L))
    this.lambda = all.lambdas[j]
    
    this.likelihoods = rep(0, times = len.G)
    
    
    for (k in seq_along(Gammas)) {
      this.gamma = Gammas[k]
      this.Z = BC(Z, this.gamma) #Y.test is not used
      
      this.fit = glmnet(X, this.Z)
      
      
      
      this.lik = get.profile.lik(this.gamma, X, this.Z,
                                 this.fit, this.lambda)
      this.likelihoods[k] = this.lik
    }
    
    ind.gamma = which.max(this.likelihoods)
    gamma.hat = Gammas[ind.gamma]
    all.gamma.hats[j] = gamma.hat
  }
  
  
  all.lambda.mins = rep(0, times = len.G)
  all.lambda.1ses = rep(0, times = len.G)
  
  for (k in seq_along(Gammas)) {
    this.gamma = Gammas[k]
    this.Z = BC(Z, this.gamma) #Y.test is not used
    
    this.fit = cv.glmnet(X, this.Z)
    
    this.lambda.min = this.fit$lambda.min
    all.lambda.mins[k] = this.lambda.min
    this.lambda.1se = this.fit$lambda.1se
    all.lambda.1ses[k] = this.lambda.1se
  }
  
  
  output = list(
    gammas = all.gamma.hats,
    lambdas = all.lambdas,
    lambda.mins = all.lambda.mins,
    lambda.1ses = all.lambda.1ses
  )
  return(output)
}, cl = cl)


#Terminate parallelization
stopCluster(cl)


data.gamma.hat = c()

for (i in seq_along(sim.output)) {
  this.info = sim.output[[i]]
  this.gamma.hats = this.info$gammas
  this.lambdas = this.info$lambdas
  this.iterations = rep(i, times = length(this.lambdas))
  this.data = data.frame(iteration = this.iterations,
                         lambda = this.lambdas,
                         gamma = this.gamma.hats)
  data.gamma.hat = rbind(data.gamma.hat, this.data)
}

plot.gamma.hat = ggplot(data.gamma.hat,
                        aes(
                          x = lambda,
                          y = gamma,
                          group = iteration,
                          colour = iteration
                        )) +
  geom_line() +
  geom_hline(yintercept = gamma.0) +
  ggtitle(
    paste0(
      "Fitted BC parameters for various lambda values ",
      "across multiple simulated datasets"
    )
  )

plot(plot.gamma.hat)


data.lambda.hat = c()
for (i in seq_along(sim.output)) {
  this.info = sim.output[[i]]
  this.len = length(this.info$lambda.mins)
  lambda.type = rep(c("min", "1se"), each = this.len)
  this.lambda.hats = c(this.info$lambda.mins,
                       this.info$lambda.1ses)
  this.iterations = rep(i, times = this.len)
  this.data = data.frame(
    iteration = rep(this.iterations, times = 2),
    gamma = rep(Gammas, times = 2),
    type = lambda.type,
    lambda = this.lambda.hats
  )
  data.lambda.hat = rbind(data.lambda.hat, this.data)
}

mean.lambda.hats = data.lambda.hat %>%
  group_by(gamma, type) %>%
  summarise(ave = mean(lambda))
mean.lambda.hats = data.frame(mean.lambda.hats)

plot.lambda.hat = ggplot() +
  geom_line(aes(
    x = gamma,
    y = lambda,
    group = interaction(iteration, type),
    colour = iteration,
    linetype = type
  ),
  data = data.lambda.hat) +
  geom_line(
    aes(
      x = gamma,
      y = ave,
      group = type,
      linetype = type
    ),
    colour = "red",
    size = 1.1,
    data = mean.lambda.hats
  ) +
  ggtitle(
    paste0(
      "CV selected lambdas for various gamma values ",
      "across multiple simulated datasets"
    )
  ) +
  geom_vline(xintercept = gamma.0)
plot(plot.lambda.hat)


### Find which gamma values are chosen at the CV lambda values
lambda.gamma.0 = filter(data.lambda.hat, gamma == gamma.0)
lambda.gamma.0$gamma.hat = 0
for (i in seq_len(K)) {
  print(paste0(i, " of ", K))
  this.lambdas = filter(lambda.gamma.0, iteration == i)
  
  set.seed(30795084 + 10000 * i)
  
  ### Generate data inside loop
  X = matrix(runif(n * p, 0, 10), nrow = n, ncol = p)
  beta = c(rep(1, times = q), rep(0, times = p - q))
  mu.Y.raw = X %*% beta
  mu.Y = mu.Y.raw + 10 * abs(min(mu.Y.raw))
  Y = mu.Y + rnorm(n, 0, sigma)
  
  Z = inv.BC(Y, gamma.0)
  
  
  these.lambdas = this.lambdas$lambda
  
  these.gamma.hats = rep(0, times = length(these.lambdas))
  
  
  for (j in seq_along(these.lambdas)) {
    # print(paste0(i, ":", j, " of ",
    #              K, ":", len.L))
    this.lambda = these.lambdas[j]
    
    this.likelihoods = rep(0, times = len.G)
    
    
    for (k in seq_along(Gammas)) {
      this.gamma = Gammas[k]
      this.Z = BC(Z, this.gamma) #Y.test is not used
      
      this.fit = glmnet(X, this.Z)
      
      
      
      this.lik = get.profile.lik(this.gamma, X, this.Z,
                                 this.fit, this.lambda)
      this.likelihoods[k] = this.lik
    }
    
    ind.gamma = which.max(this.likelihoods)
    gamma.hat = Gammas[ind.gamma]
    ind.data = 2 * (i - 1) + j
    lambda.gamma.0$gamma.hat[ind.data] = gamma.hat
  }
}

plot.gamma.0 = ggplot(lambda.gamma.0, 
                      aes(x = lambda, y = gamma.hat,
                          colour = type)) +
  geom_line()
plot(plot.gamma.0)

print(Sys.time() - time)