##################################################
### Run this after "Make One Prof Lik - LS.R"  ###
### Plots the (discretized) profile likelihood ###
##################################################


this.prof.lik = unlist(sim.output)

results = data.frame(lik = this.prof.lik,
                     gamma = Gammas)

plot.lik = ggplot(data = results, mapping = aes(x = gamma, y = lik)) +
  geom_line() + xlab("gamma") + ylab("Profile Likelihood")
plot(plot.lik)



### Find optimizer of profile likelihood
opt.lik = optimize(prof.lik.ls, c(gamma.min, gamma.max), X=X, Y=Y,
                   maximum = T)
gamma.hat = opt.lik$maximum


### Find stepdown CI for gamma ###
max.lik = opt.lik$objective
thresh = max.lik - CI.step.size
lik.left = prof.lik.ls(gamma.min, X, Y)
lik.right = prof.lik.ls(gamma.max, X, Y)
if(lik.left > thresh) {
  # Interval extends below candidate region
  warning(paste0("Confidence interval extends below candidate region for gamma",
                 " at iteration: ", j))
  a = gamma.min
} else{
  # Find left endpoint of interval inside candidate region
  a = uniroot(lik.root.ls,
              c(gamma.min, gamma.hat),
              X = X, Y = Y, val = thresh)$root
}
if(lik.right > thresh) {
  # Interval extends above candidate region
  warning(paste0("Confidence interval extends above candidate region for gamma",
                 " at iteration: ", j))
  b = gamma.max
} else{
  # Find right endpoint of interval inside candidate region
  b = uniroot(lik.root.ls,
              c(gamma.hat, gamma.max),
              X = X, Y = Y, val = thresh)$root
}
CI = c(a,b)


### Annotate profile likelihood plot with:
### -gamma.0
### -gamma.hat
### -CI for gamma.0

# annot.plot.lik = plot.lik + geom_vline(xintercept = gamma.0, colour = "red", lwd=2) +
#   geom_vline(xintercept = gamma.hat) + geom_vline(xintercept = CI, lty = 2) +
#   ggtitle(plot.title)
# plot(annot.plot.lik)

pres.plot1 = plot.lik + geom_vline(xintercept = gamma.hat)
plot(pres.plot1)
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
