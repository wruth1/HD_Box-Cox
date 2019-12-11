library(MASS)
library(glmnet)
library(pbapply)
library(tidyverse)
library(reshape2)


set.seed(96533752)

source("LASSO_Likelihood_Helper_Functions.R")

M = 5 #Number of times to replicate simulation

n = 100
p = 100
q = 10
sigma = 1
gamma.0 = 1

#Smallest and largest gamma candidates
gamma.min = -3
gamma.max = 5
#Step size for gamma candidates
gamma.step = 0.5
Gammas = seq(gamma.min, gamma.max, gamma.step)
len.G = length(Gammas)

### Store all profile likelihoods
### Dimensions are:
### Type = LASSO Min or LASSO 1se
### Iteration
### Gamma (i.e. BC parameter)
pr.lik.collection = array(0, dim = c(2, M, len.G))


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
  set.seed(100 * 5 ^ j)
  
  all.lik.lasso = pbsapply(seq_along(Gammas), function(i) {
    this.gamma = Gammas[i]
    this.info = profile.lik.lasso(this.gamma, X, Z)
    this.lik = this.info[[1]]
    this.lambdas = this.info[[2]]
    all.lambdas[,j,i] = this.lambdas
    this.vars = this.info[[3]]
    all.vars[1,i,] <<- increment.counts(all.vars[1,i,], this.vars[[1]])
    all.vars[2,i,] <<- increment.counts(all.vars[2,i,], this.vars[[2]])
    return(this.lik)
  })
  
  pr.lik.collection[,j,] = all.lik.lasso
  
}

### Compute SD within each parameter setting
pr.lik.sd = apply(pr.lik.collection, c(1, 3), function(W) {
  return(sd(W) / sqrt(length(W)))
})


var.probs = all.vars/M
vars.tidy = melt(var.probs)
colnames(vars.tidy) = c("Type", "Gamma", "Variable", "Prob")
vars.tidy$Type = c("Min", "1SE")[vars.tidy$Type]
vars.tidy$Gamma = Gammas[vars.tidy$Gamma]

vars.min = vars.tidy %>% filter(Type == "Min")
vars.1se = vars.tidy %>% filter(Type == "1SE")
  
select.plot.min = ggplot(data= vars.min, aes(x=Gamma, y = Variable)) +
  geom_tile(aes(fill = Prob), colour = "white") + 
  scale_fill_gradient(low = "white", high = "steelblue") + 
  geom_hline(yintercept = q, size=2) + 
  geom_vline(xintercept = gamma.0)
plot(select.plot.min)

select.plot.1se = ggplot(data= vars.1se, aes(x=Gamma, y = Variable)) +
  geom_tile(aes(fill = Prob), colour = "white") + 
  scale_fill_gradient(low = "white", high = "steelblue")+ 
  geom_hline(yintercept = q, size=2) + 
  geom_vline(xintercept = gamma.0)
plot(select.plot.1se)
