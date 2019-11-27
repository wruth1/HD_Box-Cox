library(MASS)
library(hdi)
library(glmnet)

set.seed(45945726)

data(riboflavin)
data.rib = riboflavin

Y = data.rib$y
X.raw = data.rib$x
X = scale(X.raw)

n = nrow(X)
p = ncol(X)
l = sqrt(2*log(p)/n)

fit.cv = cv.glmnet(X, Y)
l.min = fit.cv$lambda.min
l.1se = fit.cv$lambda.1se

sd.Y = sd(Y)


########################
### Helper Functions ###
########################

get.Y.lik = function(Y, X, b, sigma, tau){
  L = prod(Y)^tau
  L = L * prod################################################################################
}

#Computes the profile likelihood for a LASSO fit to the univariate
#Box-Cox transformed Y values with transformation parameter gamma
get.profile.likelihood = function(gamma){
  #Get Box-Cox transformed Y
  Y.new = get.new.Y(Y, gamma)
  
  #Fit LASSO model and extract predictions
  fit = glmnet(X, Y.new)
  Y.hat = predict(fit, X, s = l.min)
  
  #Compute log-likelihood (up to constants not depending on l)
  sd.Y.new = sd(Y.new)
  A = - sum((Y.new - Y.hat)^2) / 2*sd.Y.new^2
  B = - 0.5*n*log(2*pi*sd.Y.new^2)
  C = prod(Y^(gamma - 1))
  
  ### Finish constructing likelihood ##############################################################################
  
  
  return(log.lik)
}

get.new.Y = function(Y.raw, gamma){
  S = sign(Y.raw)
  Y = abs(Y.raw)
  if(gamma == 0){
    Z = log(Y)
  } else{
    Z = Y^gamma - 1
    Z = Z / gamma
  }
  Z = Z * S
  return(Z)
}



########################
### Perform Analysis ###
########################

#Do CV LASSO to choose candidate active set
fit.lasso = cv.glmnet(X, Y)
l.min = fit.lasso$lambda.min
l.1se = fit.lasso$lambda.1se
b.hat.min = predict(fit.lasso, X, s=l.min, type = "coefficients")
b.hat.1se = predict(fit.lasso, X, s=l.1se, type = "coefficients")
A.hat.min = which(b.hat.min != 0)
A.hat.1se = which(b.hat.1se != 0)

m = 21
# Box-Cox parameter candidates
gamma.seq = seq(-5, 2, length.out = m)

log.lik.seq = sapply(seq_along(gamma.seq), function(i){
  this.gamma = gamma.seq[i]
  this.log.lik = get.profile.likelihood(this.gamma, A.hat.1se)
  return(this.log.lik)
})

plot(gamma.seq, log.lik.seq)



#Plot residuals for various gamma values
for(i in seq_along(gamma.seq)){
  this.gamma = gamma.seq[i]
  this.Y = get.new.Y(Y, this.gamma)
  fit = lm(this.Y ~ X[,A.hat.1se])
  plot(fit, which=1, caption = paste0("gamma = ", this.gamma))
}
