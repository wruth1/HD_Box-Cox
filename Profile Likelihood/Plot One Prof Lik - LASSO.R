#Separate losses and active set
all.losses.raw = sapply(seq_along(sim.output), function(i){
  return(sim.output[[i]][[1]])
})
all.A.hats.raw = lapply(sim.output, function(out){
  return(out[[2]])
})


####################################
### Investigate different losses ###
####################################
cat(1)

all.losses.list = data.frame(t(all.losses.raw))
all.losses = apply(all.losses.list, 2, unlist)
results = data.frame(all.losses, gamma = Gammas)

# 
# plot.trans = ggplot(data = results, aes(x = Gammas, y = transformed)) +
#   geom_line() + xlab("Gamma") + ylab("MSE") +
#   ggtitle("MSE Computed on the BC-Transformed Data") +
#   geom_vline(xintercept = gamma.0)
# plot(plot.trans)
# 
# plot.obs = ggplot(data = results, aes(x = Gammas, y = observed)) +
#   geom_point() + xlab("Gamma") + ylab("MSE") +
#   ggtitle(paste0("MSE Computed on the Observed Data Scale", ". ",
#                  "I.e. After inverting the BC-Transformation")) +
#   geom_vline(xintercept = gamma.0)
# plot(plot.obs)



plot.lik = ggplot(data = results,
                  aes(x = Gammas, y = lik)) +
  geom_point() + xlab("Gamma") + ylab("Profile Likelihood") +
  ggtitle(plot.title) + 
  geom_vline(xintercept = gamma.0)

jpeg(paste0("Plots/Profile Likelihood/LASSO/Plain/JPEG/",
            plot.title, ".jpeg"))
plot(plot.lik)
dev.off()

png(paste0("Plots/Profile Likelihood/LASSO/Plain/PNG/",
           plot.title, ".png"))
plot(plot.lik)
dev.off()

# 
# #Decompose likelihood into MSE component and Jacobian
# data.lik = results %>%
#   mutate(MSE = -n*log(transformed)/2) %>%
#   gather(key = "Type", value = "Score", MSE, Jacob, lik) %>%
#   select(-observed, -transformed)
# gamma.hat = Gammas[which.min(results$lik)]
# plot.decomp = ggplot(data.lik,
#                      aes(x = gamma, y = Score, color = Type)) +
#   geom_line() + geom_vline(xintercept = gamma.0) +
#   geom_vline(xintercept = gamma.hat, linetype=2)
# plot(plot.decomp)


###############################
### Investigate active sets ###
###############################
cat(2)

all.A.hats = lapply(all.A.hats.raw, paste)
A.hats.factor = factor(unlist(all.A.hats))

plot.lik = ggplot(data = results,
                  aes(x = Gammas, y = lik, colour = A.hats.factor)) +
  geom_point() + xlab("Gamma") + ylab("Profile Likelihood") +
  ggtitle(plot.title) + 
  geom_vline(xintercept = gamma.0) +
  theme(legend.position = "none")


jpeg(paste0("Plots/Profile Likelihood/LASSO/Active Sets/JPEG/",
            plot.title, ".jpeg"))
plot(plot.lik)
dev.off()

png(paste0("Plots/Profile Likelihood/LASSO/Active Sets/PNG/",
           plot.title, ".png"))
plot(plot.lik)
dev.off()

# plot(Gammas, a)