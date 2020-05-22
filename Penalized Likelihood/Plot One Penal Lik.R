##################################################
### Plots the (discretized) profile likelihood ###
##################################################



### Find optimizer of profile likelihood
opt.lik = optimize(pen.lik.CV.lasso, c(gamma.min, gamma.max), X=X, Y=Y,
                   all.lambdas = all.lambdas, folds = folds, 
                   lambda.type = lambda.type, maximum = T)
print(opt.lik)
gamma.hat = opt.lik$maximum



### Find stepdown CI for gamma ###
max.lik = opt.lik$objective
thresh = max.lik - CI.step.size
lik.left = pen.lik.CV.lasso(gamma.min, X, Y, folds = folds, 
                             lambda.type = lambda.type, 
                            all.lambdas = all.lambdas)
lik.right = pen.lik.CV.lasso(gamma.max, X, Y, folds = folds, 
                              lambda.type = lambda.type, 
                             all.lambdas = all.lambdas)
if(lik.left > thresh) {
  # Interval extends below candidate region
  warning(paste0("Confidence interval extends below candidate region for gamma",
                 " at iteration: ", j))
  a = gamma.min
} else{
  # Find left endpoint of interval inside candidate region
  a = uniroot(pen.lik.root.CV.lasso,
              c(gamma.min, gamma.hat), X = X, Y = Y, folds=folds,
               val = thresh, lambda.type = lambda.type, 
              all.lambdas = all.lambdas)$root
}
if(lik.right > thresh) {
  # Interval extends above candidate region
  warning(paste0("Confidence interval extends above candidate region for gamma",
                 " at iteration: ", j))
  b = gamma.max
} else{
  # Find right endpoint of interval inside candidate region
  b = uniroot(pen.lik.root.CV.lasso,
              c(gamma.hat, gamma.max), X = X, Y = Y, folds=folds,
              val = thresh, lambda.type = lambda.type, 
              all.lambdas = all.lambdas)$root
}
CI = c(a,b)

pen.lik.CV.lasso(b, X, Y, folds, all.lambdas = all.lambdas)

### Annotate profile likelihood plot with:
### -gamma.0
### -gamma.hat
### -CI for gamma.0

annot.plot.lik = plot.lik + 
  # geom_vline(xintercept = gamma.0, colour = "red", lwd=1.5) +
  geom_vline(xintercept = gamma.hat) + geom_vline(xintercept = CI, lty = 2) +
  ggtitle(plot.title)
plot(annot.plot.lik)


g = Gammas[which.max(prof.lik)]

plot(plot.lik)
pres.plot1 = plot.lik + geom_vline(xintercept = gamma.hat)
plot(pres.plot1)
pres.plot11 = pres.plot1 + geom_vline(xintercept = g, lty = 2)
plot(pres.plot11)
pres.plot2 = pres.plot1 + geom_hline(yintercept = max.lik, lty = 2)
plot(pres.plot2)
pres.plot3 = pres.plot1 + geom_hline(yintercept = thresh, lty = 2)
plot(pres.plot3)
pres.plot4 = pres.plot3 + geom_vline(xintercept = CI, lty=2)
plot(pres.plot4)
pres.plot5 = pres.plot1 + geom_vline(xintercept = CI, lty=2)
plot(pres.plot5)
pres.plot6 = pres.plot5 + geom_vline(xintercept = gamma.0, lwd=2, colour = "red")
plot(pres.plot6)


