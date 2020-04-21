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

# Computes a value with two-sided tail prob for Y
# less than tail.prob. See general notes in Overleaf
get.int = function(n, delta, sigma, tail.prob = 0.01) {
  v = delta ^ 2 + sigma ^ 2
  a = sqrt(2 * v * log(2 * n))
  b = sqrt(-2 * v * log(tail.prob))
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
      sigma = sigma
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



