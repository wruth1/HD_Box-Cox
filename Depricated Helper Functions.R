profile.lik = function(gamma, X, Y) {
  n = length(Y)
  Y.new = BC(Y, gamma)
  fit = lm(Y.new ~ X)
  mse = sum((fit$residuals) ^ 2) / n
  lik = -n * log(mse) / 2
  Jacob = (gamma - 1) * sum(log(Y))
  lik = lik + Jacob
  return(lik)
}


### Compute `profile likelihood' at gamma based on LASSO
get.profile.lik = function(gamma, X, Y, fit, lambda) {
  Y.hat = predict(fit, X, s = lambda)
  Y.obs = inv.BC(Y, gamma)
  lik = profile.lik.formula(Y.obs, Y, Y.hat, gamma)
  return(lik)
}






### Old version
# profile.lik.lasso = function(gamma, X, Y, folds = NULL) {
#   n = length(Y)
#   Y.new = BC(Y, gamma)
#   if (is.null(folds)) {
#     fit.cv = cv.glmnet(X, Y.new)
#   } else{
#     fit.cv = cv.glmnet(X, Y.new, foldid = folds)
#   }
#   Y.hat.min = predict(fit.cv, X, s = "lambda.min")
#   Y.hat.1se = predict(fit.cv, X, s = "lambda.1se")
#   resid.min = Y.new - Y.hat.min
#   resid.1se = Y.new - Y.hat.1se
#   sse.min = sum((resid.min) ^ 2) / n
#   sse.1se = sum((resid.1se) ^ 2) / n
#   lik.min = -n * log(sse.min) / 2
#   lik.1se = -n * log(sse.1se) / 2
#   Jacob = (gamma - 1) * sum(log(Y))
#   lik.min = lik.min + Jacob
#   lik.1se = lik.1se + Jacob
#   this.liks = c(lik.min, lik.1se)
#   this.lambdas = c(fit.cv$lambda.min, fit.cv$lambda.1se)
#   this.vars = list(
#     predict(fit.cv, X, s = "lambda.min", type = "nonzero"),
#     predict(fit.cv, X, s = "lambda.1se", type = "nonzero")
#   )
#   output = list(this.liks, this.lambdas, this.vars)
#   return(output)
# }

increment.counts = function(counts, inds.list) {
  inds = unlist(inds.list)
  counts[inds] = counts[inds] + 1
  return(counts)
}
