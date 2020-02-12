#########################################################
### Which gamma values are chosen at each lambda, and ###
###   which lambda values are chosen at each gamma.   ###
###      Also, find gamma hat at CV lambda values     ###
#########################################################


library(MASS)
library(glmnet)
library(pbapply)
library(tidyverse)
library(reshape2)
library(egg)
library(doParallel)
library(glmpath)


time = Sys.time()

set.seed(74799272)

source("LASSO_Likelihood_Helper_Functions.R")

K = 12 #Number of datasets to simulate

n = 100     #Sample Size
p = 200     #Number of predictors
q = 10      #Number of active predictors
sigma = 1   #Residual SD
gamma.0 = 2 #Correct value for Box-Cox transformation

#Smallest and largest gamma candidates
gamma.min = 1
gamma.max = 3
#Step size for gamma candidates
gamma.step = 0.1
#Candidate gamma values
Gammas = seq(gamma.min, gamma.max, gamma.step)
len.G = length(Gammas)

n.lambda.raw = 100


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
  library(glmpath)
  
  source("LASSO_Likelihood_Helper_Functions.R")
})


sim.output = pblapply(seq_len(K), function(i) {
  set.seed(30795084 + 10000 * i)
  
  #Generate linear model data
  X = matrix(rnorm(n * p, 0, 1), nrow = n, ncol = p)
  beta = c(rep(1, times = q), rep(0, times = p - q))
  mu.Y.raw = X %*% beta
  mu.Y = mu.Y.raw + 10 * abs(min(mu.Y.raw))
  Y = mu.Y + rnorm(n, 0, sigma)
  
  ### Transform Y so that gamma.0 is correct BC parameter
  Z = inv.BC(Y, gamma.0)
  
  ### Find lambda values to use for all datasets
  ### Note: lambda is chosen as a proportion of max(beta.hat.ls)
  X.std = scale(X)
  fit.raw = glmnet(X.std, Z)
  # fit.raw = glmpath(X.std, Z, family = "gaussian",
  #                   max.norm = 10000*ncol(X),
  #                   min.lambda = 0)
  ### Extract selected lambda values and convert to proportions
  all.lambda.vals = fit.raw$lambda
  max.lambda = max(all.lambda.vals)
  all.lambda.props = all.lambda.vals / max.lambda
  
  len.L = length(all.lambda.props)
  
  #Contrainer for fitted BC parameters
  all.gamma.hats = rep(0, times = len.L)
  
  ##############################################
  ### Estimate gamma at each value of lambda ###
  ##############################################
  for (j in seq_len(len.L)) {
    # print(paste0(i, ":", j, " of ",
    #              K, ":", len.L))
    this.lambda.prop = all.lambda.props[j]
    
    #Container for profile likelihoods
    this.likelihoods = rep(0, times = len.G)
    
    #Compute profile likelihoods
    for (k in seq_along(Gammas)) {
      this.gamma = Gammas[k]
      
      ### Transform observed data using this gamma value
      this.Z = BC(Z, this.gamma)
      
      ### Fit LASSO to transformed response
      this.fit = glmnet(X.std, this.Z, family = "gaussian")
      # Extra arguments for glmpath
      # max.norm = 10000*ncol(X),
      # min.lambda = 0)
      
      ### Get predictions at this value of lambda
      ### (Use lambda.prop times largest lambda selected by glmnet)
      this.max.lambda = max(abs(this.fit$lambda))
      this.lambda = this.max.lambda * this.lambda.prop
      
      ### Get and store profile likelihood value
      this.lik = get.profile.lik(this.gamma, X.std, this.Z,
                                 this.fit, this.lambda)
      this.likelihoods[k] = this.lik
    }
    
    ind.gamma = which.max(this.likelihoods)
    gamma.hat = Gammas[ind.gamma]
    all.gamma.hats[j] = gamma.hat
  }
  
  
  all.lambda.min.props = rep(0, times = len.G)
  all.lambda.1se.props = rep(0, times = len.G)
  
  for (k in seq_along(Gammas)) {
    this.gamma = Gammas[k]
    this.Z = BC(Z, this.gamma)
    
    this.fit = cv.glmnet(X.std, this.Z, family = "gaussian")
    
    this.max.lambda = max(abs(this.fit$lambda))
    
    
    
    this.lambda.min.prop = this.fit$lambda.min / this.max.lambda
    all.lambda.min.props[k] = this.lambda.min.prop
    this.lambda.1se.prop = this.fit$lambda.1se / this.max.lambda
    all.lambda.1se.props[k] = this.lambda.1se.prop
  }
  
  
  output = list(
    gamma.hats = all.gamma.hats,
    lambdas = all.lambda.props,
    lambda.mins = all.lambda.min.props,
    lambda.1ses = all.lambda.1se.props
  )
  return(output)
}, cl = cl)


#Terminate parallelization
stopCluster(cl)


#Container for lambda and fitted gamma values
data.sim = c()
#Container for gamma and CV-fitted lambda values
data.sim.cv.raw = c()


#Convert simulation output to a data frame
for (i in seq_along(sim.output)) {
  this.info = sim.output[[i]]
  
  #Extract and store gamma.hat vs lambda
  this.gamma.hats = this.info$gamma.hats
  this.lambda.props = this.info$lambdas
  this.iterations = rep(i, times = length(this.lambda.props))
  this.data = data.frame(iteration = this.iterations,
                         lambda.prop = this.lambda.props,
                         gamma = this.gamma.hats)
  data.sim = rbind(data.sim, this.data)
  
  #Extract and store lambda.hat vs gamma
  this.lambda.hats.min = this.info$lambda.mins
  this.lambda.hats.1se = this.info$lambda.1ses
  this.iterations.cv = rep(i, times = length(this.lambda.hats.min))
  this.data.cv = data.frame(
    iteration = this.iterations.cv,
    gamma = Gammas,
    lambda.min = this.lambda.hats.min,
    lambda.1se = this.lambda.hats.1se
  )
  data.sim.cv.raw = rbind(data.sim.cv.raw, this.data.cv)
}

### Convert CV data from wide to tall
data.sim.cv = data.sim.cv.raw %>%
  #Convert two lambda columns into one lambda column
  #and a label column
  melt(
    measure.vars = c("lambda.min", "lambda.1se"),
    value.name = "lambda",
    variable.name = "type.raw"
  ) %>%
  #Shorten lambda labels
  mutate(type = str_extract(type.raw, "...$")) %>%
  select(-type.raw)

##################################################
### Plot gamma.hat as a function of lambda     ###
### Note: lambda is expressed as a prop of the ###
###       largest value selected by glmnet     ###                              ###
##################################################
plot.gamma.hat = ggplot(data.sim,
                        aes(
                          x = log(lambda.prop),
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

mean.lambda.hats = data.sim.cv %>%
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
  data = data.sim.cv) +
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
      "across multiple simulated datasets",
      "\n",
      "(Dotted for lambda.min, solid for lambda.1se; ",
      "Red lines are averages across iterations)"
    )
  ) +
  geom_vline(xintercept = gamma.0)
plot(plot.lambda.hat)


############################################################
### Find which gammas are chosen at the CV lambda values ###
############################################################

lambda.gamma.0 = filter(data.sim.cv, gamma == gamma.0)
lambda.gamma.0$gamma.hat = 0
for (i in seq_len(K)) {
  print(paste0(i, " of ", K))
  this.lambdas = filter(lambda.gamma.0, iteration == i)
  
  #Generate the same data used to fit lambda
  set.seed(30795084 + 10000 * i)
  X = matrix(runif(n * p, 0, 10), nrow = n, ncol = p)
  beta = c(rep(1, times = q), rep(0, times = p - q))
  mu.Y.raw = X %*% beta
  mu.Y = mu.Y.raw + 10 * abs(min(mu.Y.raw))
  Y = mu.Y + rnorm(n, 0, sigma)
  
  Z = inv.BC(Y, gamma.0)
  
  X.std = scale(X)
  
  these.lambdas = this.lambdas$lambda
  
  these.gamma.hats = rep(0, times = length(these.lambdas))
  
  
  for (j in seq_along(these.lambdas)) {
    this.lambda = these.lambdas[j]
    
    #Container to store profile likelihoods for the current lambda
    this.likelihoods = rep(0, times = len.G)
    
    #Compute the profile likelihoods for the current lambda
    for (k in seq_along(Gammas)) {
      this.gamma = Gammas[k]
      this.Z = BC(Z, this.gamma)
      
      this.fit = glmnet(X.std, this.Z)
      
      
      this.lik = get.profile.lik(this.gamma, X.std, this.Z,
                                 this.fit, this.lambda)
      this.likelihoods[k] = this.lik
    }
    
    ind.gamma = which.max(this.likelihoods)
    gamma.hat = Gammas[ind.gamma]
    
    #Position in lambda.gamma.0 where gamma.hat goes
    ind.data = 2 * (i - 1) + j
    
    lambda.gamma.0$gamma.hat[ind.data] = gamma.hat
  }
}

plot.gamma.0 = ggplot(lambda.gamma.0,
                      aes(x = lambda, y = gamma.hat,
                          colour = type)) + geom_line() +
  geom_hline(yintercept = gamma.0)
plot(plot.gamma.0)

print(Sys.time() - time)