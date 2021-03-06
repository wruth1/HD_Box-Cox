library(optimx)
library(ggplot2)
library(dplyr)
library(numDeriv)

#####################
### Least Squares ###
#####################

### Create functions to compute profile likelihood and its grad
my.prof.lik = function(gamma) prof.lik.ls(gamma, X, Y)

my.grad.prof.lik = function(gamma) {
   info = prof.lik.ls(gamma, X, Y, grad=T)
   this.grad = info$grad
   return(this.grad)
}



### Optimize the profile likelihood
test.ls = optimr(par = 2, fn = my.prof.lik, 
              gr = my.grad.prof.lik,
              method = "BFGS",
              control = list(maximize = T))


### Compute profile likelihood and its gradient over Gammas
liks = sapply(Gammas, my.prof.lik)
grads = sapply(Gammas, my.grad.prof.lik)
LS.data = data.frame(gamma = Gammas,
                     lik = liks,
                     grad = grads)


### Compute some intermediate quantities for optimizing prof lik
max.lik = max(liks)
thresh = max.lik - 2
CI = filter(LS.data, lik >= thresh)
a = min(CI$gamma)
b = max(CI$gamma)
crit = Gammas[which.min(abs(grads))]


### Plot profile likelihood with maximizer, gamma.0 and CI
LS.plot = ggplot(mapping = aes(x = Gammas,
                               y = liks))+
   geom_line() + geom_vline(xintercept = crit) +
   geom_vline(xintercept = gamma.0, color = "red") +
   geom_hline(yintercept = thresh, linetype = "dashed") +
   geom_vline(xintercept = c(a, b), linetype = "dashed")
plot(LS.plot)

### Plot gradients of profile likelihood
plot(Gammas, grads)



#############
### LASSO ###
#############

### Get grid of candidate lambdas, increase resolution,
### and get lambda.1se from CV
fit0 = glmnet(X, Y)
lambda.0 = fit0$lambda
lambda.1 = finer.grid(lambda.0)
fit = cv.glmnet(X, Y, lambda = lambda.1)
l = fit$lambda.1se


### Create functions to compute profile likelihood and its grad
lasso.prof.lik = function(gamma) prof.lik.lasso(gamma, X, Y, l)
lasso.penal.prof.lik = function(gamma) penal.prof.lik.lasso(gamma, X, Y, l)

### Optimize profile likelihood using numerical derivative
test.lasso = optimize(prof.lik.lasso, c(gamma.min, gamma.max),
                      maximum = T, X=X, Y=Y, l=l)
crit = test.lasso$maximum
test.penal.lasso = optimize(penal.prof.lik.lasso, c(gamma.min, gamma.max),
                      maximum = T, X=X, Y=Y, l=l)
crit.penal = test.penal.lasso$maximum

### Compute profile likelihood and its grad over Gammas
liks = sapply(Gammas, lasso.prof.lik)

### Compute penalized profile likelihood
pen.liks = sapply(Gammas, lasso.penal.prof.lik)

### Compute approximate grad using numerical method
num.grad = sapply(Gammas, function(g){
   return(grad(lasso.prof.lik, g, "simple"))
})




### Plot profile likelihood with critical value
plot(Gammas, liks)
abline(v = crit)
abline(v = gamma.0, col="red")

### Plot penalized profile likelihood with critical value
plot(Gammas, pen.liks)
abline(v = crit.penal)
abline(v = gamma.0, col="red")


### Plot analytic and numeric gradients of prof lik
plot(Gammas, grads)
plot(Gammas, num.grad)

