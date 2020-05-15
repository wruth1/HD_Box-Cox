library(MASS)
library(glmnet)
library(pbapply)
library(grDevices)
library(leaps)

set.seed(96533752)

M = 5 #Number of times to replicate simulation

n = 100
p = 1000
q = 10
sigma = 1
beta = c(rep(1, times = q), rep(0, times = p - q))
gamma.0 = 1

#Smallest and largest gamma candidates
gamma.min = -1
gamma.max = 3
#Step size for gamma candidates
gamma.step = 0.05
Gammas = seq(gamma.min, gamma.max, gamma.step)
len.G = length(Gammas)

### Store all profile likelihoods
### Dimensions are:
### Type = LS, AIC
### Iteration
### Gamma (i.e. BC parameter)
pr.lik.collection = array(0, dim = c(2, M, len.G))



#source(Roboflavin_Investigation.R)
profile.lik = function(gamma, X, Y) {
  n = length(Y)
  Y.new = BC(Y, gamma)
  fit = lm(Y.new ~ X[, 1:q])
  sse = sum((fit$residuals) ^ 2)
  lik = -n * log(sse) / 2
  Jacob = (gamma - 1) * sum(log(Y))
  lik = lik + Jacob
  return(lik)
}

profile.lik.AIC = function(gamma, X, Y) {
  n = length(Y)
  Y.new = BC(Y, gamma)
  data = data.frame(X, Y.new)
  fit.LS = lm(Y.new ~ ., data = data)
  fit.AIC = stepAIC(fit.LS,
                    scope = list(upper = ~ ., lower = ~ 1),
                    trace = 0)
  Y.hat = predict(fit.AIC, data.frame(X))
  sse = sum((fit.AIC$residuals) ^ 2)
  lik = -n * log(sse) / 2
  Jacob = (gamma - 1) * sum(log(Y))
  lik = lik + Jacob
  return(lik)
}

#Returns a vector that the BC transform (with par. gamma) maps to Y
#Note: Does not account for negative values
inv.BC = function(Y, gamma) {
  if (gamma == 0) {
    return(exp(Y))
  } else{
    Z = 1 + gamma * Y
    Z = Z ^ (1 / gamma)
    return(Z)
  }
}

#Returns the BC transformation of Y with par. gamma
#Note: Does not account for negative values
BC = function(Y, gamma) {
  if (gamma == 0) {
    return(log(Y))
  } else{
    Z = Y ^ gamma
    Z = (Z - 1) / gamma
    return(Z)
  }
}


lm.largest = function(X, Y){
  n = nrow(X)
  p = ncol(X)
  
  Y.new = BC(Z, gamma)
  data = data.frame(X, Y.new)
  
  fit.step = regsubsets(Z ~ ., data = data, nbest = 1,
                    method = "forward",
                    nvmax = n/2)
  BICs.step = summary(fit.step)$bic
  errs.step = summary(fit.step)$rss
  AICs.step = n*log(errs.step) + 2*(1:length(errs.step)+1)
}

for (j in 1:M) {
  set.seed(100 * 5 ^ j)
  
  all.lambdas = c()
  
  X = matrix(rnorm(n * p, 10, 10), nrow = n, ncol = p)
  mu.Y.raw = X %*% beta
  mu.Y = mu.Y.raw + 10 * abs(min(mu.Y.raw))
  Y = mu.Y + rnorm(n, 0, sigma)
  
  Z = inv.BC(Y, gamma.0)
  
  # plot(X, Y)
  # plot(X, Z)
  
  
  
  all.likelihoods = sapply(seq_along(Gammas), function(i) {
    this.gamma = Gammas[i]
    this.lik = profile.lik(this.gamma, X, Z)
    return(this.lik)
  })
  
  pr.lik.collection[1, j,] = all.likelihoods
  
  all.lik.AIC = pbsapply(seq_along(Gammas), function(i) {
    this.gamma = Gammas[i]
    this.lik = profile.lik.AIC(this.gamma, X, Z)
    return(this.lik)
  })
  
  pr.lik.collection[2, j,] = all.lik.AIC
  
}

### Compute mean and SE within each parameter setting
pr.lik.sd = apply(pr.lik.collection, c(1, 3), function(W) {
  return(sd(W) / sqrt(length(W)))
})


for (j in 1:M) {
  plot(
    Gammas,
    pr.lik.collection[1, j,],
    type = "l",
    xlim = c(gamma.0 - 1, gamma.0 + 1),
    xlab = "BC Parameter",
    ylab = "\"log-likelihood\"",
    main = paste0("Profile log-likelihood for BC parameter under Oracle - ",
                  j)
  )
  lines(Gammas,
        pr.lik.collection[1, j,] + pr.lik.sd[1,],
        col = "Red",
        lty = 2)
  lines(Gammas,
        pr.lik.collection[1, j,] - pr.lik.sd[1,],
        col = "Red",
        lty = 2)
  abline(v = 1, lwd = 1, col = "blue")
}

for (j in 1:M) {
  plot(
    Gammas,
    pr.lik.collection[2, j,],
    type = "l",
    xlim = c(gamma.0 - 1, gamma.0 + 1),
    #ylim = c(-800,-300),
    xlab = "BC Parameter",
    ylab = "\"log-likelihood\"",
    main = paste0(
      "\"Profile log-likelihood\" for BC parameter under AIC Selection - ",
      j
    )
  )
  lines(Gammas,
        pr.lik.collection[2, j,] + pr.lik.sd[2,],
        col = "Red",
        lty = 2)
  lines(Gammas,
        pr.lik.collection[2, j,] - pr.lik.sd[2,],
        col = "Red",
        lty = 2)
  abline(v = 1, lwd = 1, col = "blue")
  
}
