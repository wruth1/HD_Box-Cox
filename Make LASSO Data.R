
library(glmnet)
library(doParallel)
library(pbapply)
library(stringr)

set.seed(51798431)

time = Sys.time()

source("LASSO_Likelihood_Helper_Functions.R")


n = 100     #Sample Size
beta.size = 20

### Sizes of steps down from optimizer for profile likelihood CIs
step.sizes = 1:12

### Target number of Ys to generate for each X
### Actual value is slightly smaller to optimize parallelization
M.target = 50

#Smallest and largest gamma candidates
gamma.min = 1
gamma.max = 5
#Step size for gamma candidates
gamma.step = 0.01
#Candidate gamma values
Gammas = seq(gamma.min, gamma.max, gamma.step)
len.G = length(Gammas)

n.folds = 10

sigmas = c(0.1, 0.01)
gamma.0s = c(2, 3, 4)
ps = c(10, 50)#, 200)
deltas = c(1, 10) #SD of a term in X %*% beta

#Construct list of parameter combinations
all.pars = expand.grid(
  p = ps,
  sigma = sigmas,
  gamma.0 = gamma.0s,
  delta = deltas
)

j = 1
# j = 5

pars = all.pars[j, ]
attach(pars)


### Construct coefficient vector
q = floor(sqrt(p))


M = 1

all.data = make.data(n, p, q, delta, sigma, M)#, return.pars=T)


beta.size = delta/sqrt(q)
beta = c(rep(beta.size, q), rep(0, p-q))
X = all.data[[1]]$X
Z = all.data[[1]]$Z
Y = inv.BC(Z, gamma.0)
# beta = all.data[[2]]


detach(pars)





H = X %*% solve(t(X) %*% X) %*% t(X)
det(H)
G = diag(n) - H
test = eigen(G, symmetric = T)
norm(G, "2")