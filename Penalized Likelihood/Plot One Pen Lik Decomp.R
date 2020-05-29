#########################################################################
### Plots the profile penalized likelihood, along with the following: ###
###   -Unpenalized profile likelihood                                 ###
###   -log lambda hat                                                 ###
###   -log L1 norm of beta hat                                        ###
###   -Number of active predictors                                    ###
#########################################################################


### Plot profile penalized likelihood
plot.lik = ggplot(data = results, 
                  mapping = aes(x = gamma, y = lik)) +
  geom_point(size = 0.5) + xlab("gamma") + ylab("Prof Pen Lik") +
  geom_vline(xintercept = gamma.0, lty = 2) +
  ggtitle(plot.title)
# plot(plot.lik)

### Plot unpenalized profile likelihood
plot.lik.raw = ggplot(data = results, 
                      mapping = aes(x = gamma, y = lik.raw)) +
  geom_point(size = 0.5) + xlab("gamma") + ylab("Prof Lik") +
  geom_vline(xintercept = gamma.0, lty = 2)
# plot(plot.lik.raw)

### Plot log fitted lambda value
plot.l.hat = ggplot(data = results, aes(x=gamma, y=log(l.hat))) +
  geom_point(size = 0.5) + xlab("gamma") + ylab("log(Fitted Lambda)") +
  geom_vline(xintercept = gamma.0, lty = 2)
# plot(plot.l.hat)


### Plot log L1 norm of beta hat
plot.norm.b = ggplot(data = results, aes(x=gamma, y=b.norm)) +
  geom_point(size = 0.5) + xlab("gamma") + ylab("log(L1 Norm)") +
  geom_vline(xintercept = gamma.0, lty = 2)
# plot(plot.norm.b)


### Plot active set size
plot.active = ggplot(data = results, aes(x = gamma, y = act)) +
  geom_point(size = 0.5) + xlab("gamma") + ylab("Active Set Size") +
  geom_vline(xintercept = gamma.0, lty = 2)

### Stack above plots into one figure and plot it
plot.lik.facet = ggplotGrob(plot.lik)
plot.lik.raw.facet = ggplotGrob(plot.lik.raw)
plot.l.hat.facet = ggplotGrob(plot.l.hat)
plot.norm.b.facet = ggplotGrob(plot.norm.b)
plot.active.facet = ggplotGrob(plot.active)
plot.lik.active = rbind(plot.lik.facet, plot.lik.raw.facet, plot.l.hat.facet,
                        plot.norm.b.facet, plot.active.facet)

file.name = paste0("Plots/Pen Prof Lik Decomp/", CV.type, "/", plot.title, ".jpg")
jpeg(file.name)
plot(plot.lik.active)
dev.off()
