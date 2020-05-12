source("Helper Functions/Profile Likelihood Helper Functions.R")

geom_mean = function(x, na.rm = TRUE) {
  exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
}

#Extracts the specified variable from simulation output
get.from.sim = function(sim.output, var.name) {
  vals.raw = lapply(sim.output, function(G) {
    values = G[[1]]
    this.val = values[var.name]
    return(this.val)
  })
  vals = unlist(vals.raw)
}

### Construct a coefficient vector with the specified size
### Size can refer to the size of the non-zero terms (type=="val")
### of the 2-norm of beta (type == "norm")
make.beta = function(p, q, size = 0, type = "val") {
  if (type == "val") {
    b = rep(size, times = q)
    b = c(b, rep(0, times = p - q))
  } else if (type == "norm") {
    val = size / sqrt(q)
    b = rep(val, times = q)
    b = c(b, rep(0, times = p - q))
  } else{
    stop("In make.beta (Invalid type specified)")
  }
}





#Returns a vector that the BC transform (with par. gamma) maps to Y
#Note: Does not account for negative values
inv.BC = function(Z, gamma, type = "good") {
  if (type == "good") {
    if (gamma == 0) {
      return(exp(Z))
    } else{
      Y = 1 + gamma * Z
      Y = Y ^ (1 / gamma)
      return(Y)
    }
  } else if (type == "simple") {
    if (gamma == 0) {
      return(exp(Z))
    } else{
      Y = Z ^ (1 / gamma)
      return(Y)
    }
  }
}

#Returns the BC transformation of Y with par. gamma
#Note: Does not account for negative values
BC = function(Y, gamma, type = "good") {
  if (type == "good") {
    if (gamma == 0) {
      return(log(Y))
    } else{
      Z = Y ^ gamma
      Z = (Z - 1) / gamma
      return(Z)
    }
  } else if (type == "simple") {
    if (gamma == 0) {
      return(log(Y))
    } else{
      Z = Y ^ gamma
      return(Z)
    }
  }
}

# Computes a value with one-sided tail prob for Y
# less than tail.prob. See general notes in Overleaf
get.int = function(n, delta, sigma, tail.prob = 0.01) {
  v = delta ^ 2 + sigma ^ 2
  a = sqrt(2 * v * log(2 * n))
  b = sqrt(-2 * v * log(2*tail.prob))
  return(a + b)
}

### Generates M (X, Y) pairs with the specified parameters
### Optionally also returns the specified parameter settings
### Note: X is a matrix of iid standard normals
### Note: delta is the SD of an entry in X %*% beta
### Note: tail.prob upper bounds P(Y < 0)
make.data = function(n,
                     p,
                     q,
                     delta,
                     sigma,
                     M = 1,
                     tail.prob = 0.01,
                     return.pars = F) {
  all.data = lapply(seq_len(M), function(j) {
    X = matrix(rnorm(n * p), nrow = n, ncol = p)
    # Size of each non-zero element of beta
    beta.term = delta / sqrt(q)
    beta = c(rep(beta.term, times = q),
             rep(0, times = p - q))
    intercept = get.int(n, delta, sigma, tail.prob)
    Z.mean = X %*% beta + intercept
    epsilon = rnorm(n, sd = sigma)
    Z = Z.mean + epsilon
    data = list(X = X, Z = Z)
    return(data)
  })
  if (!return.pars) {
    return(all.data)
  } else{
    pars = list(
      n = n,
      p = p,
      q = q,
      delta = delta,
      sigma = sigma,
      beta = beta
    )
    output = list(data = all.data, pars = pars)
    return(output)
  }
}

### Computes the gradient of the BC transformation wrt gamma
grad.BC.gamma = function(Y, gamma) {
  if (gamma == 0) {
    c = log(Y)
    return(c ^ 2 / 2)
  } else{
    a = Y ^ gamma * log(Y)
    b = BC(Y, gamma, type = "good")
    output = (a - b) / gamma
    return(output)
  }
}


### Creates a new grid based on X by replacing every interior pair
### with their midpoint
### If length(X) is odd, keep the (more precisely, a) median
coarser.grid = function(X){
  n = length(X)
  Y = X[1]
  if(n %% 2 != 0){
    med = median(X)
    ind.med = min(which(X == med))
    X = X[-ind.med]
    Y = c(Y, med)
  }
  n = length(X) - 2
  m = n/2
  Y = X[1]
  for(j in 1:m){
    a = X[2*j]
    b = X[2*j+1]
    c = (a+b)/2
    Y = c(Y, c)
  }
  Y = c(Y, X[length(X)])
  Y = sort(Y)
  return(Y)
}

### Apply coarser.grid to fine.grid as many times as possible while keeping
### n.target elements
coarsen.grid = function(n.target, fine.grid){
  N = length(fine.grid)
  passes = floor(log2(N/n.target))
  coarse.grid = fine.grid
  for(j in 1:passes){
    coarse.grid = coarser.grid(coarse.grid)
  }
  return(coarse.grid)
}

### Creates a new grid based on X by adding all the midpoints
finer.grid = function(X){
  n = length(X)
  m = 2*n - 1
  Y = rep(0, times = m)
  inds.X = 1:n
  inds.X = 2*inds.X - 1
  Y[inds.X] = X
  for(i in 1:(n-1)){
    j = 2*i
    a = Y[j-1]
    b = Y[j+1]
    Y[j] = (a+b)/2
  }
  return(Y)
}
