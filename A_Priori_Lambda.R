library(MASS)
library(glmnet)
library(pbapply)
library(grDevices)

set.seed(64501945)

M = 5 #Number of times to replicate simulation

n = 100
p = 10
q = 5
beta = c(rep(1, times = q), rep(0, times = p - q))
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
### Type = LS, LASSO
### Iteration
### Gamma (i.e. BC parameter)
pr.lik.collection = array(0, dim = c(2, M, len.G))

for (j in seq_len(M)) {
  set.seed(72045690 + 10000 * 5 ^ j)
  
  
  #source(Roboflavin_Investigation.R)
  profile.lik = function(gamma, X, Y) {
    n = length(Y)
    Y.new = BC(Y, gamma)
    fit = lm(Y.new ~ X)
    sse = sum((fit$residuals) ^ 2)/n
    lik = -n * log(sse) / 2
    Jacob = (gamma - 1) * sum(log(Y))
    lik = lik + Jacob
    return(lik)
  }
  
  profile.lik.lasso = function(gamma, X, Y) {
    Y.new = BC(Y, gamma)
    fit = glmnet(X, Y.new)
    this.lambda = 4 * sigma * sqrt(log(p) / n)
    Y.hat = predict(fit, X, s = this.lambda)
    resid = Y.new - Y.hat
    sse = sum((resid) ^ 2)/n
    lik.raw = -n * log(sse) / 2
    Jacob = (gamma - 1) * sum(log(Y))
    lik = lik.raw + Jacob
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
  
  
  
  
  
  X = matrix(rnorm(n * p, 10, 10), nrow = n, ncol = p)
  beta = 1:p
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
  
  pr.lik.collection[1, j, ] = all.likelihoods
  
  
  all.lik.lasso = pbsapply(seq_along(Gammas), function(i) {
    this.gamma = Gammas[i]
    this.lik = profile.lik.lasso(this.gamma, X, Z)
    return(this.lik)
  })
  
  pr.lik.collection[2, j, ] = all.lik.lasso
  
}

### Compute mean and SE within each parameter setting
pr.lik.mean = apply(pr.lik.collection, c(1, 3), function(W) {
  return(mean(W))
})
pr.lik.se = apply(pr.lik.collection, c(1, 3), function(W) {
  return(sd(W) / sqrt(length(W)))
})

for (j in 1:M) {
  ### Plot "log-likelihood function" using residuals from
  ### lm and from LASSO
  plot(
    Gammas,
    pr.lik.collection[1, j, ],
    type = "l",
    # ylim = c(-625, -200),
    xlab = "BC Parameter",
    ylab = "\"log-likelihood\"",
    main = "Profile log-likelihood for BC parameter under LS"
  )
  abline(v = 1, lwd = 1, col = "blue")
  
  plot(
    Gammas,
    pr.lik.collection[2, j, ],
    type = "l",
    # ylim = c(-625, -200),
    xlab = "BC Parameter",
    ylab = "\"log-likelihood\"",
    main = "\"Profile log-likelihood\" for BC parameter under LASSO"
  )
  abline(v = 1, lwd = 1, col = "blue")
  
  
  
  ### Plot "log-likelihood function" using residuals from
  ### lm and from LASSO, ZOOMED-IN around the maximum
  plot(
    Gammas,
    pr.lik.collection[1, j, ],
    type = "l",
    xlim = c(gamma.0 - 1, gamma.0 + 1),
    xlab = "BC Parameter",
    ylab = "\"log-likelihood\"",
    main = "Profile log-likelihood for BC parameter under LS"
  )
  abline(v = 1, lwd = 1, col = "blue")
  
  plot(
    Gammas,
    pr.lik.collection[2, j, ],
    type = "l",
    xlim = c(gamma.0 - 1, gamma.0 + 1),
    xlab = "BC Parameter",
    ylab = "\"log-likelihood\"",
    main = "\"Profile log-likelihood\" for BC parameter under LASSO"
  )
  abline(v = 1, lwd = 1, col = "blue")
  
}

### Plot average profile likelihood with +/- 1SE

plot(
  Gammas,
  pr.lik.mean[1, ],
  type = "l",
  xlim = c(gamma.0 - 1, gamma.0 + 1),
  xlab = "BC Parameter",
  ylab = "\"log-likelihood\"",
  main = "Mean Profile log-likelihood for BC parameter under LS"
)
lines(Gammas,
      pr.lik.mean[1,] + pr.lik.se[1, ],
      col = "Red",
      lty = 2)
lines(Gammas,
      pr.lik.mean[1,] - pr.lik.se[1, ],
      col = "Red",
      lty = 2)
abline(v = 1, lwd = 1, col = "blue")

plot(
  Gammas,
  pr.lik.mean[2,],
  type = "l",
  xlim = c(gamma.0 - 1, gamma.0 + 1),
  ylim = c(-800,-300),
  xlab = "BC Parameter",
  ylab = "\"log-likelihood\"",
  main = "\"Profile log-likelihood\" for BC parameter under LASSO Min"
)
lines(Gammas,
      pr.lik.mean[2,] + pr.lik.se[2, ],
      col = "Red",
      lty = 2)
lines(Gammas,
      pr.lik.mean[2,] - pr.lik.sd[2,],
      col = "Red",
      lty = 2)
abline(v = 1, lwd = 1, col = "blue")
