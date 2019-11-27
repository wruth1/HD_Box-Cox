library(MASS)
library(glmnet)
library(pbapply)

set.seed(96533752)


all.lambdas = c()

#source(Roboflavin_Investigation.R)
profile.lik = function(gamma, X, Y){
  n = length(Y)
  Y.new = BC(Y, gamma)
  fit = lm(Y.new ~ X)
  sse = sum((fit$residuals)^2)
  lik = -n*log(sse)/2
  Jacob = (gamma - 1) * sum(log(Y))
  lik = lik + Jacob
  return(lik)
}

profile.lik.lasso = function(gamma, X, Y, type = "1se"){
  n = length(Y)
  Y.new = BC(Y, gamma)
  fit.cv = cv.glmnet(X, Y.new)
  lambda = paste0("lambda.", type)
  all.lambdas <<- cbind(all.lambdas, c(fit.cv$lambda.1se,
                                       fit.cv$lambda.min))
  Y.hat = predict(fit.cv, X, s = lambda)
  resid = Y.new - Y.hat
  sse = sum((resid)^2)
  lik = -n*log(sse)/2
  Jacob = (gamma - 1) * sum(log(Y))
  lik = lik + Jacob
  return(lik)
}

#Returns a vector that the BC transform (with par. gamma) maps to Y
#Note: Does not account for negative values
inv.BC = function(Y, gamma){
  if(gamma == 0){
    return(exp(Y))
  } else{
  Z = 1 + gamma*Y
  Z = Z^(1/gamma)
  return(Z)
  }
}

#Returns the BC transformation of Y with par. gamma
#Note: Does not account for negative values
BC = function(Y, gamma){
  if(gamma ==0){
    return(log(Y))
  } else{
  Z = Y^gamma
  Z = (Z - 1)/gamma
  return(Z)
  }
}


n = 100
p = 10
sigma = 1
gamma.0 = 2

#Smallest and largest gamma candidates
gamma.min = -2
gamma.max = 4
#Step size for gamma candidates
gamma.step = 0.01


X = matrix(rnorm(n*p, 10, 10), nrow = n, ncol = p)
beta = 1:p
mu.Y.raw = X %*% beta
mu.Y = mu.Y.raw + 10*abs(min(mu.Y.raw))
Y = mu.Y + rnorm(n, 0, sigma)

Z = inv.BC(Y, gamma.0)

# plot(X, Y)
# plot(X, Z)

Gammas = seq(gamma.min, gamma.max, gamma.step)
m = length(Gammas)

all.likelihoods = sapply(seq_along(Gammas),function(i){
  this.gamma = Gammas[i]
  this.lik = profile.lik(this.gamma, X, Z)
  return(this.lik)
})

all.lik.lasso = pbsapply(seq_along(Gammas),function(i){
  this.gamma = Gammas[i]
  this.lik = profile.lik.lasso(this.gamma, X, Z, "min")
  return(this.lik)
})


plot(Gammas, all.likelihoods, type="l")
plot(Gammas, all.lik.lasso, type = "l")



fit = lm(Z ~ X)
test = boxcox(fit, lambda = seq(-2, 4, 1/10))



fit.1se = lm(log(all.lambdas[1,]) ~ Gammas)
summary(fit.1se)
plot(fit.1se, which=1)
plot(Gammas, log(all.lambdas[1,]))
abline(fit.1se)

fit.min = lm(log(all.lambdas[2,]) ~ Gammas)
summary(fit.min)
plot(fit.min, which=1)
plot(Gammas, log(all.lambdas[2,]))
abline(fit.min)


### Trim Gamma Values with weird log-lambdas
ind.drop = which(Gammas > 1 & Gammas < 3)
Gammas.trim = Gammas[-ind.drop]
all.lambdas.trim = all.lambdas[,-ind.drop]

### Re-fit regression models on trimmed data
fit.1se.trim = lm(log(all.lambdas.trim[1,]) ~ Gammas.trim)
summary(fit.1se.trim)
plot(fit.1se.trim, which=1)
plot(Gammas.trim, log(all.lambdas.trim[1,]))
abline(fit.1se.trim)
plot(Gammas, log(all.lambdas[1,]))
abline(fit.1se.trim)

fit.min.trim = lm(log(all.lambdas.trim[2,]) ~ Gammas.trim)
summary(fit.min.trim)
plot(fit.min.trim, which=1)
plot(Gammas.trim, log(all.lambdas.trim[2,]))
abline(fit.min.trim)
plot(Gammas, log(all.lambdas[2,]))
abline(fit.min.trim)
