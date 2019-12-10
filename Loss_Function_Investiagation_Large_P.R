library(MASS)
library(glmnet)
library(pbapply)
library(grDevices)

set.seed(96533752)

# all.ps = c(2, 

for(j in seq_len(5)){
  set.seed(100*5^j)

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

profile.lik.lasso = function(gamma, X, Y, type = "1se", folds=NULL){
  n = length(Y)
  Y.new = BC(Y, gamma)
  if(is.null(folds)){
    fit.cv = cv.glmnet(X, Y.new)
  } else{
    fit.cv = cv.glmnet(X, Y.new, foldid = folds)
  }
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

#Returns a vector that the BC transform (with par. = gamma) maps to Y
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
p = 1000
q = 10
sigma = 1
gamma.0 = 1

#Smallest and largest gamma candidates
gamma.min = -3
gamma.max = 5
#Step size for gamma candidates
gamma.step = 0.05


X = matrix(rnorm(n*p, 10, 10), nrow = n, ncol = p)
beta = c(rep(1, times=q), rep(0, times = p-q))
mu.Y.raw = X %*% beta
mu.Y = mu.Y.raw + 10*abs(min(mu.Y.raw))
Y = mu.Y + rnorm(n, 0, sigma)

Z = inv.BC(Y, gamma.0)

# plot(X, Y)
# plot(X, Z)

Gammas = seq(gamma.min, gamma.max, gamma.step)
m = length(Gammas)

### Get CV folds
fit.lasso = cv.glmnet(X, Y, keep=T)
folds = fit.lasso$foldid

all.lik.lasso = pbsapply(seq_along(Gammas),function(i){
  this.gamma = Gammas[i]
  this.lik = profile.lik.lasso(this.gamma, X, Z, "min", folds)
  return(this.lik)
})


### Plot "log-likelihood function" using residuals from LASSO
jpeg(paste0("Plots/Profile Likelihood - LASSO/",
            "Profile Likelihood LASSO (P=1000) - ",
            j, ".jpg"))
plot(Gammas, all.lik.lasso, type = "l",# ylim = c(-625, -200),
     xlab = "BC Parameter", ylab = "\"log-likelihood\"",
     main = "\"Profile log-likelihood\" for BC parameter under LASSO - P=1000")
abline(v = 1, lwd = 1, col = "blue")
dev.off()

### Plot "log-likelihood function" using residuals from
### lm and from LASSO focused around the maximum
jpeg(paste0("Plots/Zoomed Profile Likelihood - LASSO/",
            "Zoomed Profile Likelihood LASSO (P=1000) - ",
            j, ".jpg"))
plot(Gammas, all.lik.lasso, type = "l",
     xlim = c(gamma.0 - 1, gamma.0 + 1),
     xlab = "BC Parameter", ylab = "\"log-likelihood\"",
     main = "\"Profile log-likelihood\" for BC parameter under LASSO - P=1000")
abline(v = 1, lwd = 1, col = "blue")
dev.off()


### Model log(lambda) ~ BC par. using linear regression
### Plot results
fit.1se = lm(log(all.lambdas[1,]) ~ Gammas)
jpeg(paste0("Plots/Log-Lambda Min/",
            "Log Lambda Min (P=1000) - ",
            j, ".jpg"))
plot(Gammas, log(all.lambdas[1,]), type="l", lwd=2,
     xlab = "BC Parameter", ylab = "log(lambda)",
     main = "lambda.1se chosen by CV for various BC parameters - P=1000")
abline(v = 1, lwd = 1, col = "blue")
abline(fit.1se)
dev.off()

fit.min = lm(log(all.lambdas[2,]) ~ Gammas)
jpeg(paste0("Plots/Log-Lambda 1SE/",
            "Log Lambda 1SE (P=1000) - ",
            j, ".jpg"))
plot(Gammas, log(all.lambdas[2,]), type="l", lwd=2,
     xlab = "BC Parameter", ylab = "log(lambda)",
     main = "lambda.min chosen by CV for various BC parameters - P=1000")
# abline(v = c(0.55, 1.45), lwd = 2, col = "red")
abline(v = 1, lwd = 1, col = "blue")
abline(fit.min)
dev.off()

}