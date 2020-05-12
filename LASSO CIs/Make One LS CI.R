# prof.file = "RProfile.txt"
# Rprof(prof.file)
# 
# time = Sys.time()
# 
# library(glmnet)
# library(doParallel)
# library(pbapply)
# library(stringr)

# set.seed(10782897)


# source("LASSO_Likelihood_Helper_Functions.R")


all.intervals = lapply(seq_len(M), function(i) {
  # if(i == 19) browser()
  # print(paste0(j, ":", i))
  ### Generate data
  X = matrix(rnorm(n * p, 0, 1), nrow = n, ncol = p)
  mu.Z.raw = X %*% beta
  mu.Z = mu.Z.raw + get.int(n, delta, sigma)
  Z = mu.Z + rnorm(n, 0, sigma)
  
  ### Transform Y so that gamma.0 is correct BC parameter
  Y = inv.BC(Z, gamma.0)
  
  ### Find lambda values to use for all datasets
  ### Note: lambda is chosen as a proportion of max(beta.hat.ls)
  # X.std = scale(X)
  X.std = X
  
  
  
  ##############################################
  ### Find likelihood CIs using optimization ###
  ### and root finding functions             ###
  ##############################################
  
  # opt.lik = optimr(par = 1, fn = prof.lik.ls, 
  #                  gr = grad.prof.lik.ls,
  #                  method = "BFGS",
  #                  control = list(maximize = T),
  #                  X = X, Y = Y)
  
  ### Find maximizer of profile likelihood
  opt.lik = optimize(prof.lik.ls, c(gamma.min, gamma.max),
                     X=X, Y=Y, maximum=T)

  gamma.hat = opt.lik$maximum
  lik.hat = opt.lik$objective
  
  ### Compute CIs for several likelihood drop sizes
  threshs = lik.hat - step.sizes
  ints = lapply(threshs, function(thresh) {
    ### Find left endpoint of interval
    # Check if endpoints are within the range being considered
    lik.left = prof.lik.ls(gamma.min, X, Y)
    lik.right = prof.lik.ls(gamma.max, X, Y)
    if(lik.left > thresh) {
      a = gamma.min
    } else{
      a = uniroot(lik.root.ls,
                  c(gamma.min, gamma.hat),
                  X = X, Y = Y, val = thresh)$root
    }
    if(lik.right > thresh) {
      b = gamma.max
    } else{
      b = uniroot(lik.root.ls,
                  c(gamma.hat, gamma.max),
                  X = X, Y = Y, val = thresh)$root
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

write.table(results, "LASSO CIs/Coverages - LS.csv", append = T,
            row.names = F, quote = F, sep = ",",
            col.names = F)

# print(Sys.time() - time)

# Rprof(NULL)
# 
# summaryRprof(prof.file)
