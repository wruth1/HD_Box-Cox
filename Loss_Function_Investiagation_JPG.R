library(MASS)
library(glmnet)
library(pbapply)
library(grDevices)

set.seed(96533752)

M = 5 #Number of times to replicate simulation

n = 100
p = 10
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
    ### Type = LS, LASSO Min or LASSO 1se
    ### Iteration
    ### Gamma (i.e. BC parameter)
pr.lik.collection = array(0, dim = c(3, M, len.G))

for(j in seq_len(M)){
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

profile.lik.lasso = function(gamma, X, Y, folds=NULL){
  n = length(Y)
  Y.new = BC(Y, gamma)
  if(is.null(folds)){
    fit.cv = cv.glmnet(X, Y.new)
  } else{
    fit.cv = cv.glmnet(X, Y.new, foldid = folds)
  }
  all.lambdas <<- cbind(all.lambdas, c(fit.cv$lambda.1se,
                                       fit.cv$lambda.min))
  Y.hat.min = predict(fit.cv, X, s = "lambda.min")
  Y.hat.1se = predict(fit.cv, X, s = "lambda.min")
  resid.min = Y.new - Y.hat.min
  resid.1se = Y.new - Y.hat.1se
  sse.min = sum((resid.min)^2)
  sse.1se = sum((resid.1se)^2)
  lik.min = -n*log(sse.min)/2
  lik.1se = -n*log(sse.1se)/2
  Jacob = (gamma - 1) * sum(log(Y))
  lik.min = lik.min + Jacob
  lik.1se = lik.1se + Jacob
  output = c(lik.min, lik.1se)
  return(output)
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





X = matrix(rnorm(n*p, 10, 10), nrow = n, ncol = p)
beta = 1:p
mu.Y.raw = X %*% beta
mu.Y = mu.Y.raw + 10*abs(min(mu.Y.raw))
Y = mu.Y + rnorm(n, 0, sigma)

Z = inv.BC(Y, gamma.0)

# plot(X, Y)
# plot(X, Z)



all.likelihoods = sapply(seq_along(Gammas),function(i){
  this.gamma = Gammas[i]
  this.lik = profile.lik(this.gamma, X, Z)
  return(this.lik)
})

pr.lik.collection[1,j,] = all.likelihoods

### Get CV folds
fit.lasso = cv.glmnet(X, Y, keep=T)
folds = fit.lasso$foldid

all.lik.lasso = pbsapply(seq_along(Gammas),function(i){
  this.gamma = Gammas[i]
  this.lik = profile.lik.lasso(this.gamma, X, Z, folds)
  return(this.lik)
})

pr.lik.collection[2:3,j,] = all.lik.lasso


### Compute mean and SE within each parameter setting
pr.lik.means = apply(pr.lik.collection, c(1,3), mean)
pr.lik.sd = apply(pr.lik.collection, c(1,3), function(W){
  return(sd(W)/sqrt(length(W)))
})


### Plot "log-likelihood function" using residuals from
### lm and from LASSO
jpeg(paste0("Plots/Profile Likelihood - LS/Profile Likelihood LS - ",
            j, ".jpg"))
plot(Gammas, pr.lik.means[1,], type="l",# ylim = c(-625, -200),
     xlab = "BC Parameter", ylab = "\"log-likelihood\"",
     main = "Profile log-likelihood for BC parameter under LS")
abline(v = 1, lwd = 1, col = "blue")
dev.off()
jpeg(paste0("Plots/Profile Likelihood - LASSO/",
            "Profile Likelihood LASSO Min - ",
            j, ".jpg"))
plot(Gammas, pr.lik.means[2,], type = "l",# ylim = c(-625, -200),
     xlab = "BC Parameter", ylab = "\"log-likelihood\"",
     main = "\"Profile log-likelihood\" for BC parameter under LASSO Min")
abline(v = 1, lwd = 1, col = "blue")
dev.off()
jpeg(paste0("Plots/Profile Likelihood - LASSO/",
            "Profile Likelihood LASSO 1se - ",
            j, ".jpg"))
plot(Gammas, pr.lik.means[3,], type = "l",# ylim = c(-625, -200),
     xlab = "BC Parameter", ylab = "\"log-likelihood\"",
     main = "\"Profile log-likelihood\" for BC parameter under LASSO 1se")
abline(v = 1, lwd = 1, col = "blue")
dev.off()

### Plot "log-likelihood function" using residuals from
### lm and from LASSO focused around the maximum
jpeg(paste0("Plots/Zoomed Profile Likelihood - LS/",
            "Zoomed Profile Likelihood LS - ",
            j, ".jpg"))
plot(Gammas, pr.lik.means[1,], type="l",
     xlim = c(gamma.0 - 1, gamma.0 + 1),
     xlab = "BC Parameter", ylab = "\"log-likelihood\"",
     main = "Profile log-likelihood for BC parameter under LS")
abline(v = 1, lwd = 1, col = "blue")
dev.off()
jpeg(paste0("Plots/Zoomed Profile Likelihood - LASSO/",
            "Zoomed Profile Likelihood LASSO min - ",
            j, ".jpg"))
plot(Gammas, pr.lik.means[2,], type = "l",
     xlim = c(gamma.0 - 1, gamma.0 + 1),
     xlab = "BC Parameter", ylab = "\"log-likelihood\"",
     main = "\"Profile log-likelihood\" for BC parameter under LASSO Min")
abline(v = 1, lwd = 1, col = "blue")
dev.off()
jpeg(paste0("Plots/Zoomed Profile Likelihood - LASSO/",
            "Zoomed Profile Likelihood LASSO 1se - ",
            j, ".jpg"))
plot(Gammas, pr.lik.means[3,], type = "l",
     xlim = c(gamma.0 - 1, gamma.0 + 1),
     xlab = "BC Parameter", ylab = "\"log-likelihood\"",
     main = "\"Profile log-likelihood\" for BC parameter under LASSO 1se")
abline(v = 1, lwd = 1, col = "blue")
dev.off()


### Model log(lambda) ~ BC par. using linear regression
### Plot results
fit.1se = lm(log(all.lambdas[1,]) ~ Gammas)
jpeg(paste0("Plots/Log-Lambda Min/",
            "Log Lambda Min - ",
            j, ".jpg"))
plot(Gammas, log(all.lambdas[1,]), type="l", lwd=2,
     xlab = "BC Parameter", ylab = "log(lambda)",
     main = "lambda.1se chosen by CV for various BC parameters")
abline(v = 1, lwd = 1, col = "blue")
abline(fit.1se)
dev.off()

fit.min = lm(log(all.lambdas[2,]) ~ Gammas)
jpeg(paste0("Plots/Log-Lambda 1SE/",
            "Log Lambda 1SE - ",
            j, ".jpg"))
plot(Gammas, log(all.lambdas[2,]), type="l", lwd=2,
     xlab = "BC Parameter", ylab = "log(lambda)",
     main = "lambda.min chosen by CV for various BC parameters")
# abline(v = c(0.55, 1.45), lwd = 2, col = "red")
abline(v = 1, lwd = 1, col = "blue")
abline(fit.min)
dev.off()

}