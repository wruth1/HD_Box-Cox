######################################################################
### Compute the BC parameter estimator from Klaassen et al. (2017) ###
######################################################################

library(ggplot2)
library(glmnet)

source("LASSO_Likelihood_Helper_Functions.R")

set.seed(74380493)




time = Sys.time()

K = 12 #Number of CV folds

n = 100
p = 10
q = 10
sigma = 1
gamma.0 = 2

#Smallest and largest gamma candidates
gamma.min = gamma.0 - 1
gamma.max = gamma.0 + 1
#Step size for gamma candidates
gamma.step = 0.01
Gammas = seq(gamma.min, gamma.max, gamma.step)
len.G = length(Gammas)


### Generate data inside loop
X = matrix(rnorm(n * p, 10, 10), nrow = n, ncol = p)
# beta = c(rep(1, times = q), rep(0, times = p - q))
beta = c(1:q, rep(0, times = p - q))
mu.Y.raw = X %*% beta
mu.Y = mu.Y.raw + 10 * abs(min(mu.Y.raw))
Y = mu.Y + rnorm(n, 0, sigma)

Z = inv.BC(Y, gamma.0)

fit.raw = glmnet(X, Z)
all.lambdas = fit.raw$lambda
len.L = length(all.lambdas)

all.gamma.hats = rep(0, times = length(all.lambdas))
all.liks = array(0, dim = c(len.L, len.G))


for (j in seq_along(all.lambdas)) {
  print(paste0(j, " of ",
               len.L))
  this.lambda = all.lambdas[j]
  
  this.likelihoods = rep(0, times = len.G)
  
  
  for (k in seq_along(Gammas)) {
    this.gamma = Gammas[k]
    
    this.lik = profile.lik.lasso(this.gamma, X, Z, this.lambda)
    
    this.likelihoods[k] = this.lik
  }
  
  all.liks[j, ] = this.likelihoods
  
  ind.gamma = which.max(this.likelihoods)
  gamma.hat = Gammas[ind.gamma]
  all.gamma.hats[j] = gamma.hat
}

single.curve = ggplot(mapping = aes(y = all.liks[1, ],
                                    x = Gammas)) +
  geom_line() + ggtitle("A Single Profile Likelihood")
plot(single.curve)

mean.liks = apply(all.liks, 2, mean)
se.liks = apply(all.liks, 2, sd) / sqrt(len.L)

data.liks = data.frame(
  gamma = Gammas,
  mean = mean.liks,
  upper = mean.liks + se.liks,
  lower = mean.liks - se.liks
)



plot.liks = ggplot(data.liks, aes(x = gamma)) +
  geom_line(aes(y = mean)) + geom_vline(xintercept = gamma.0) +
  geom_line(aes(y = upper), colour = "red") +
  geom_line(aes(y = lower), colour = "red")
plot(plot.liks)


#
# plot.liks.realized = ggplot(mapping = aes(x = Gammas)) +
#   geom_vline(xintercept = gamma.0)
# num.curves = 5
# for(i in 1:num.curves){
#   plot.liks.realized = plot.liks.realized +
#     geom_line(aes(y = all.liks[!!i,]))
# }
# plot(plot.liks.realized)



##################
### LS version ###
##################

ls.liks = array(0, dim = c(len.G))

for (j in seq_along(Gammas)) {
  this.gamma = Gammas[j]
  
  this.lik = profile.lik.ls(this.gamma, X, Z)
  
  ls.liks[j] = this.lik
}

boxcox(Z ~ X, lambda = Gammas)

ind.gamma = which.max(ls.liks)
ls.gamma.hat = Gammas[ind.gamma]

plot.liks.ls = ggplot(mapping = aes(x = Gammas, y = ls.liks)) +
  geom_line() + geom_vline(xintercept = gamma.0) +
  ggtitle("Profile log-likelihood for LS") +
  ylab("Log Profile Likelihood")
plot(plot.liks.ls)

print(Sys.time() - time)