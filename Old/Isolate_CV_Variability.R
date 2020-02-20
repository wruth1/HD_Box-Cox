library(MASS)
library(glmnet)
library(pbapply)
library(tidyverse)
library(reshape2)
library(egg)

time = Sys.time()

set.seed(92794619)

source("LASSO_Likelihood_Helper_Functions.R")

M = 20 #Number of times to replicate simulation

n = 100
p = 100
q = 10
sigma = 1
gamma.0 = 1

#Smallest and largest gamma candidates
gamma.min = -3
gamma.max = 5
#Step size for gamma candidates
gamma.step = 0.05
Gammas = seq(gamma.min, gamma.max, gamma.step)
len.G = length(Gammas)

### Store all profile likelihoods
### Dimensions are:
### Type = Oracle, LASSO Min or LASSO 1se
### Iteration
### Gamma (i.e. BC parameter)
pr.lik.collection = array(0, dim = c(3, M, len.G))


### Store all lambda values chosen by CV
### Dimensions are:
### Type = LASSO Min or LASSO 1se
### Iteration
### Gamma
all.lambdas = array(0, dim = c(2, M, len.G))

### Count how many times each variable is selected by each LASSO
### Dimensions are:
### Type = LASSO Min or LASSO 1se
### Gamma
### Variable
all.vars = array(0, dim = c(2, len.G, p))

### Generate data outside loop
X = matrix(rnorm(n * p, 10, 10), nrow = n, ncol = p)
beta = c(rep(1, times=q), rep(0, times = p-q))
mu.Y.raw = X %*% beta
mu.Y = mu.Y.raw + 10 * abs(min(mu.Y.raw))
Y = mu.Y + rnorm(n, 0, sigma)

Z = inv.BC(Y, gamma.0)



for (j in seq_len(M)) {
  print(paste0(j, " of ", M))
  set.seed(22977417 + 10000 * j)
  
  all.lik.oracle = sapply(seq_along(Gammas), function(i) {
    this.gamma = Gammas[i]
    this.lik = profile.lik(this.gamma, X[,1:q], Z)
    return(this.lik)
  })
  
  pr.lik.collection[1,j,] = all.lik.oracle
  
  fit.cv.glmnet = cv.glmnet(X, Z)
  folds = fit.cv.glmnet$foldid
  
  all.lik.lasso = pbsapply(seq_along(Gammas), function(i) {
    this.gamma = Gammas[i]
    this.info = profile.lik.lasso(this.gamma, X, Z, folds)
    this.lik = this.info[[1]]
    this.lambdas = this.info[[2]]
    all.lambdas[,j,i] = this.lambdas
    this.vars = this.info[[3]]
    all.vars[1,i,] <<- increment.counts(all.vars[1,i,], this.vars[[1]])
    all.vars[2,i,] <<- increment.counts(all.vars[2,i,], this.vars[[2]])
    return(this.lik)
  })
  
  pr.lik.collection[2:3,j,] = all.lik.lasso
  
}

### Compute SD within each parameter setting
pr.lik.sd = apply(pr.lik.collection, c(1, 3), function(W) {
  return(sd(W) / sqrt(length(W)))
})
pr.lik.mean = apply(pr.lik.collection, c(1, 3), function(W) {
  return(mean(W))
})

#########################################################
### Process and plot variable selection probabilities ###
#########################################################

### Format probabilities for tidyverse
var.probs = all.vars/M
vars.tidy = melt(var.probs)
colnames(vars.tidy) = c("Type", "Gamma", "Variable", "Prob")
vars.tidy$Type = c("Min", "1SE")[vars.tidy$Type]
vars.tidy$Gamma = Gammas[vars.tidy$Gamma]

### Split probs by LASSO flavor
vars.min = vars.tidy %>% filter(Type == "Min")
vars.1se = vars.tidy %>% filter(Type == "1SE")
  
### Plot LASSO min probabilities
select.plot.min = ggplot(data= vars.min, aes(x=Gamma, y = Variable)) +
  geom_tile(aes(fill = Prob), colour = "white") + 
  scale_fill_gradient(low = "white", high = "blue") + 
  geom_hline(yintercept = q, size=1.5) + 
  geom_vline(xintercept = gamma.0) 
#plot(select.plot.min)

### Plot LASSO 1SE Probs
select.plot.1se = ggplot(data= vars.1se, aes(x=Gamma, y = Variable)) +
  geom_tile(aes(fill = Prob), colour = "white") + 
  scale_fill_gradient(low = "white", high = "blue")+ 
  geom_hline(yintercept = q, size=1.5) + 
  geom_vline(xintercept = gamma.0) 
#plot(select.plot.1se)


##############################################
### Process and plot "profile likelihoods" ###
##############################################

### Format "likelihoods" for tidyverse
liks.tidy = melt(pr.lik.collection)
colnames(liks.tidy) = c("Type", "Iteration", "Gamma", "Likelihood")
liks.tidy$Type = c("Oracle", "Min", "1SE")[liks.tidy$Type]
liks.tidy$Gamma = Gammas[liks.tidy$Gamma]

### Separate "likelihoods" by LASSO flavor. Keep only 1st iteration
###   Want to plot a single iteration. Mean of all iter's would
###   under-represent variability. 1se is chosen arbitrarily
liks.min = liks.tidy %>% filter(Type == "Min", Iteration == 1)
liks.1se = liks.tidy %>% filter(Type == "1SE", Iteration == 1)
liks.oracle = liks.tidy %>% filter(Type == "Oracle", Iteration == 1)

### Plot "profile log-likelihood" for LASSO MIN
lik.plot.min = ggplot(data = liks.min, aes(x=Gamma, y=Likelihood))+
  geom_line(size = 2) + geom_vline(xintercept = gamma.0) +
  geom_line(data = liks.oracle, aes(x=Gamma, y=Likelihood)) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank()) +
  ylab("\"Log-Likelihood\"") +
  ggtitle(paste0("\"Log-Profile-Likelihood\" and ",
  "Variable Selection Probabilities Under LASSO Min \n",
  "(thick: observed, thin: oracle LS)"))
#plot(lik.plot.min)

ggarrange(lik.plot.min, select.plot.min, heights = c(0.25, 0.75))


### Plot "profile log-likelihood" for LASSO 1SE
lik.plot.1se = ggplot(data = liks.1se, aes(x=Gamma, y=Likelihood))+
  geom_line(size = 2) + geom_vline(xintercept = gamma.0) +
  geom_line(data = liks.oracle, aes(x=Gamma, y=Likelihood)) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank()) +
  ylab("\"Log-Likelihood\"") +
  ggtitle(paste0("\"Log-Profile-Likelihood\" and ",
                 "Variable Selection Probabilities Under LASSO 1se \n",
                 " (thick: observed, thin: oracle LS)"))
#plot(lik.plot.1se)

ggarrange(lik.plot.1se, select.plot.1se, heights = c(0.25, 0.75))


### Zoom plots around the true gamma
zoom = xlim(0,2)
lik.plot.min.zoom = lik.plot.min + zoom
lik.plot.1se.zoom = lik.plot.1se + zoom
select.plot.min.zoom = select.plot.min + zoom
select.plot.1se.zoom = select.plot.1se + zoom

ggarrange(lik.plot.min.zoom, select.plot.min.zoom,
          heights = c(0.25, 0.75))
ggarrange(lik.plot.1se.zoom, select.plot.1se.zoom,
          heights = c(0.25, 0.75))

print(Sys.time() - time)