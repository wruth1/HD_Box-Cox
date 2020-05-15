#Separate losses and active set
all.losses.raw = sapply(seq_along(sim.output), function(i){
  return(sim.output[[i]][[1]])
})


###############################
### Extract relevant output ###
###############################
all.losses.list = data.frame(t(all.losses.raw))
all.losses = apply(all.losses.list, 2, unlist)
results = data.frame(all.losses, gamma = Gammas)




plot.lik = ggplot(data = results,
                  aes(x = gamma, y = lik)) +
  geom_point() + xlab("Gamma") + ylab("Profile Likelihood") +
  ggtitle(plot.title) +
  geom_vline(xintercept = gamma.0, color = "red",
             size = 1.5)


### Add likelihood CI to plot
opt = results[which.max(results$lik),]
gamma.hat = opt$gamma
lik.max = opt$lik
thresh = lik.max - 2
range = filter(results, lik > thresh)

plot.lik = plot.lik +
  geom_vline(xintercept = gamma.hat, linetype = "dashed") +
  geom_rug(aes(x = gamma, y = NULL), data=range) +
  geom_hline(yintercept = thresh)


jpeg(paste0("Plots/Profile Likelihood/LS/JPEG/",
            plot.title, ".jpeg"))
plot(plot.lik)
dev.off()

png(paste0("Plots/Profile Likelihood/LS/PNG/",
           plot.title, ".png"))
plot(plot.lik)
dev.off()

