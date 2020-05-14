# library(profvis)
library(optimx)

# 
# time = Sys.time()
# 
# library(glmnet)
# library(doParallel)
# library(pbapply)
# library(stringr)

# set.seed(10782897)


# source("LASSO_Likelihood_Helper_Functions.R")

# profvis({

# M = 5


all.intervals = pblapply(seq_len(M), function(i) {
  ### Generate data
  X = matrix(rnorm(n * p, 0, 1), nrow = n, ncol = p)
  mu.Z.raw = X %*% beta
  intercept = get.int(n, delta, sigma)
  mu.Z = mu.Z.raw + intercept
  Z = mu.Z + rnorm(n, 0, sigma)
  
  ### Transform Y so that gamma.0 is correct BC parameter
  Y = inv.BC(Z, gamma.0)
  
  ### Find lambda values to use for all datasets
  ### Note: lambda is chosen as a proportion of max(beta.hat.ls)
  # X.std = scale(X)
  X.std = X
  
  ### Find all candidate lambda values
  all.lambdas.raw = sapply(Gammas, function(gamma){
    this.Z = BC(Y, gamma)
    this.fit = glmnet(X, this.Z)
    this.lambdas = this.fit$lambda
    return(this.lambdas)
  })
  all.lambdas.fine = sort(unlist(all.lambdas.raw))
  
  ### Make the grid of lambda candidates coarser
  all.lambdas = coarsen.grid(n.lambda, all.lambdas.fine)
  
  
  ##############################################################################
  ### Use optimization software to maximize the profile likelihood for gamma ###
  ##############################################################################
  
  opt.lik = optimize(prof.lik.CV.lasso, c(gamma.min, gamma.max),
                     maximum = T, X=X, Y=Y, folds = folds, 
                     all.lambdas = all.lambdas)
  gamma.hat = opt.lik$maximum
  lik.hat = opt.lik$objective
  
  
  #####################################################
  ### Compute CIs for several likelihood drop sizes ###
  #####################################################
  
  ### Likelihood thresholds for inclusion in CIs
  threshs = lik.hat - step.sizes
  
  ### Likelihoods at endpoints of candidate range
  ### Used to check if interval extends outside candidate range for gamma
  lik.left = prof.lik.CV.lasso(gamma.min, X, Y, folds=folds, 
                               all.lambdas = all.lambdas)
  lik.right = prof.lik.CV.lasso(gamma.max, X, Y, folds=folds, 
                                all.lambdas = all.lambdas)
  
  ### Compute intervals
  ints = lapply(threshs, function(thresh) {
    ### Check if endpoints are within the range being considered
    ### If yes, return endpoint. If no, compute endpoint inside interval
    if(lik.left > thresh) {
      a = gamma.min
    } else{
      a = uniroot(lik.root.CV.lasso, c(gamma.min, gamma.hat),
                  X = X, Y = Y, val = thresh, folds=folds, 
                  all.lambdas = all.lambdas)$root
    }
    if(lik.right > thresh) {
      b = gamma.max
    } else{
      b = uniroot(lik.root.CV.lasso, c(gamma.hat, gamma.max),
                  X = X, Y = Y, val = thresh, folds=folds, 
                  all.lambdas = all.lambdas)$root
    }
    return(c(a,b))
  })
  
  ### Check coverage
  cover = sapply(ints, function(int) {
    a = min(int)
    b = max(int)
    check = (a <= gamma.0) && (gamma.0 <= b)
    length = b-a
    output = list(check, length)
  })
  
  
  
  
  
  
})


### Compute coverage rates
all.coverages = sapply(all.intervals, function(info){
  cover.check = info[1,]
})

# Handle the case of a single step.size separately
if(length(step.sizes == 1)){
  cover.probs = mean(unlist(all.coverages))
}else{
  cover.probs = apply(all.coverages, 1, function(checks) {
    return(mean(unlist(checks)))
  })
}
# print(paste0("Estimated coverage prob. over ", M, " datasets"))
# print(cover.probs)


### Extract interval lengths
all.lengths = sapply(all.intervals, function(info){
  lengths = info[2,]
})
# Handle the case of a single step.size separately
if(length(step.sizes == 1)){
  mean.lengths = mean(unlist(all.lengths))
} else{
  mean.lengths = apply(all.lengths, 1, function(lengths) {
    return(mean(unlist(lengths)))
  })
}


### Combine coverage probabilities with interval lengths
info = paste0(signif(cover.probs, digits = 3),
              "(", 
              signif(mean.lengths, digits = 3),
              ")")
names(info) = paste0("Minus ", step.sizes)


pars = c(n = n, p = p, q = q, sigma = sigma,
         gamma.0 = gamma.0, 
         delta = delta,
         M = M, gamma.step = gamma.step)
results.raw = c(pars, info)
results = data.frame(t(results.raw))

write.table(results, "LASSO CIs/Coverages - LASSO.csv",
            append = T, row.names = F, quote = F,
            sep = ",", col.names = F)

# })

# print(Sys.time() - time)


