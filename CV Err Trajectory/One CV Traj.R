##############################################
### Split data into training and test sets ###
##############################################

### Construct training set indices
prop.train = 0.75
size.train = floor(n*prop.train)
inds.train = sample(1:n, size.train, replace = F)

### Split X and Y into training and test
X.train = X[inds.train,]
Y.train = Y[inds.train]
X.test = X[-inds.train,]
Y.test = Y[-inds.train]


################################################
### Fit LASSO to training set, extract knots ###
################################################

### Fit model
Z.train = BC(Y.train, gamma.hat)
fit.lars = lars(X.train, Z.train)
# cv.lars(X.train, Z.train)

### Get knots as proportions of largest lambda
knot.props = fit.lars$lambda




###############################################
### Fit LASSO to training set using glmnet, ### 
### get prediction error for various lambda ###
###############################################

### Fit model and extract lambdas
fit = glmnet(X.train, Z.train)
this.lambdas = fit$lambda
knots = max(this.lambdas)*knot.props

# ### Add more lambda candidates around optimizer and re-fit model
# new.lambdas = seq(0.005, 0.015, by = 0.0002)
# lambdas = c(lambdas.raw, new.lambdas)
# lambdas = sort(lambdas, decreasing=T)
# fit = glmnet(X.train, Z.train, lambda = lambdas)

### Compute prediction errors for all lambda values on the Y-scale
preds = predict(fit, X.test)
errs = apply(preds, 2, function(Z.hat){
  Y.hat = inv.BC(Z.hat, gamma.hat)
  err = mean((Y.test - Y.hat)^2)
  return(err)
})


###################################################
### Plot MSPE as function of lambda, with knots ###
###################################################

plot.MSPE = ggplot(mapping = aes(x=this.lambdas, y=errs)) + geom_line() +
  geom_vline(xintercept = knots, lty=2) +
  # xlim(0.005, 0.015) + ylim(1, 1.3) +
  geom_rug(sides = "b")
plot(plot.MSPE)


###############################################
### Investigate quadratic fits to MSPE path ###
###############################################

data.interp.raw = data.frame(lambda = lambdas, err = errs)
data.interp = filter(data.interp.raw, lambda <=0.015, lambda >= 0.005)

interp = lm(err ~ poly(lambda, 2), data=data.interp)
err.interp = predict(interp, data.frame(lambda = lambdas))

plot.MSPE.interp = plot.MSPE + geom_line(y = err.interp, col = "red")
plot(plot.MSPE.interp)




# 
# 
# 
# ### Computes the CV error along a grid of lambda values and return both in a DF
# get.CV.traj = function(X, Y, gamma){
#   this.Z = BC(Y, gamma)
#   fit = cv.glmnet(X, this.Z)
#   errs = fit$cvm
#   lambdas = fit$lambda
#   output = data.frame(lambda = lambdas, err = errs)
#   return(output)
# }
# 
# this.traj = get.CV.traj(X, Y, 0)
# 
# traj.plot = ggplot(this.traj, aes(x=lambda, y=err)) + geom_line() +
#   geom_point(shape = 20)
# plot(traj.plot)
# 
# 
# library(lars)
# this.Z = BC(Y, 0)
# fit = lars(X, this.Z)
# 
# cv.fit = cv.lars(X, this.Z, mode = "step")
# 
# cv.err = cv.fit$cv
# lambdas = fit$lambda
# plot(lambdas, cv.err)
