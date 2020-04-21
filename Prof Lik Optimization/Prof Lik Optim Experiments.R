library(optimx)


my.prof.lik = function(gamma) prof.lik.ls(gamma, X, Y)

my.grad.prof.lik = function(gamma) {
   info = prof.lik.ls(gamma, X, Y, grad=T)
   this.grad = info$grad
   return(this.grad)
}



test = optimr(par = 2, fn = my.prof.lik, 
              gr = my.grad.prof.lik,
              method = "BFGS",
              control = list(maximize = T))


liks = my.prof.lik(Gammas, Y, this.Z, this.Z.hat)

