
lambda.hats = get.from.sim(sim.output, lambda.type)
sd.Zs = get.from.sim(sim.output, "sd.Z")
gm.Zs = get.from.sim(sim.output, "gm.Z")


data.lambda = data.frame(gamma = Gammas, lambda = lambda.hats,
                         sd = sd.Zs, gm = gm.Zs)
data.lambda = data.lambda %>% mutate(log.lambda = log(lambda),
                                     log.gm = log(gm))

lambda.plot = ggplot(data.lambda, aes(x=gamma, y = lambda)) +
  geom_line()
# plot(lambda.plot)



### Plot log(lambda.hat)
cat(0)
lambda.plot = ggplot(data.lambda, 
                         aes(x=gamma, y = lambda)) +
  xlab("Gamma") + ylab("CV Lambda") +
  ggtitle(plot.title) + 
  geom_vline(xintercept = gamma.0) +
  geom_point()

jpeg(paste0("Plots/CV lambda/lambda/JPEG/",
            plot.title, ".jpeg"))
plot(lambda.plot)
dev.off()

png(paste0("Plots/CV lambda/lambda/PNG/",
           plot.title, ".png"))
plot(lambda.plot)
dev.off()



### Plot log(lambda.hat)
cat(1)
log.lambda.plot = ggplot(data.lambda, 
                         aes(x=gamma, y = log.lambda)) +
  xlab("Gamma") + ylab("log(CV Lambda)") +
  ggtitle(plot.title) + 
  geom_vline(xintercept = gamma.0) +
  geom_point()

jpeg(paste0("Plots/CV lambda/log-lambda/JPEG/",
            plot.title, ".jpeg"))
plot(log.lambda.plot)
dev.off()

png(paste0("Plots/CV lambda/log-lambda/PNG/",
           plot.title, ".png"))
plot(log.lambda.plot)
dev.off()


### Plot lambda/SE(Y)
cat(2)
lambda.over.sd = ggplot(data.lambda, 
                        aes(x=gamma, y = lambda/sd)) +
  xlab("Gamma") + ylab("CV Lambda / SE(Y)") +
  ggtitle(plot.title) + 
  geom_vline(xintercept = gamma.0) +
  geom_point()


jpeg(paste0("Plots/CV lambda/lambda over sd/JPEG/",
            plot.title, ".jpeg"))
plot(lambda.over.sd)
dev.off()

png(paste0("Plots/CV lambda/lambda over sd/PNG/",
           plot.title, ".png"))
plot(lambda.over.sd)
dev.off()



### Plot log(lambda)/SE(Y)
cat(3)
log.lambda.over.sd = ggplot(data.lambda, 
                        aes(x=gamma, y = log.lambda/sd)) +
  xlab("Gamma") + ylab("log(CV Lambda) / SE(Y)") +
  ggtitle(plot.title) + 
  geom_vline(xintercept = gamma.0) +
  geom_point()


jpeg(paste0("Plots/CV lambda/log-lambda over sd/JPEG/",
            plot.title, ".jpeg"))
plot(log.lambda.over.sd)
dev.off()

png(paste0("Plots/CV lambda/log-lambda over sd/PNG/",
           plot.title, ".png"))
plot(log.lambda.over.sd)
dev.off()



### Plot lambda/geom_mean(Y)
cat(4)
lambda.over.gm = ggplot(data.lambda, 
                            aes(x=gamma, y = lambda/gm)) +
  xlab("Gamma") + ylab("CV Lambda / geom_mean(Y)") +
  ggtitle(plot.title) + 
  geom_vline(xintercept = gamma.0) +
  geom_point()


jpeg(paste0("Plots/CV lambda/lambda over GM(Y)/JPEG/",
            plot.title, ".jpeg"))
plot(lambda.over.gm)
dev.off()

png(paste0("Plots/CV lambda/lambda over GM(Y)/PNG/",
           plot.title, ".png"))
plot(lambda.over.gm)
dev.off()



### Plot log(lambda)/geom_mean(Y)
cat(5)
log.lambda.over.gm = ggplot(data.lambda, 
                        aes(x=gamma, y = log.lambda/gm)) +
  xlab("Gamma") + ylab("log-Lambda / geom_mean(Y)") +
  ggtitle(plot.title) + 
  geom_vline(xintercept = gamma.0) +
  geom_point()


jpeg(paste0("Plots/CV lambda/log-lambda over GM(Y)/JPEG/",
            plot.title, ".jpeg"))
plot(log.lambda.over.gm)
dev.off()

png(paste0("Plots/CV lambda/log-lambda over GM(Y)/PNG/",
           plot.title, ".png"))
plot(log.lambda.over.gm)
dev.off()



### Plot log(lambda)/log(geom_mean(Y))
cat(6)
log.lambda.over.log.gm = ggplot(data.lambda, 
                            aes(x=gamma, y = log.lambda/log.gm)) +
  xlab("Gamma") + ylab("log-Lambda / log-geom_mean(Y)") +
  ggtitle(plot.title) + 
  geom_vline(xintercept = gamma.0) +
  geom_point()


jpeg(paste0("Plots/CV lambda/log-lambda over log-GM(Y)/JPEG/",
            plot.title, ".jpeg"))
plot(log.lambda.over.log.gm)
dev.off()

png(paste0("Plots/CV lambda/log-lambda over log-GM(Y)/PNG/",
           plot.title, ".png"))
plot(log.lambda.over.log.gm)
dev.off()



### Plot log(lambda) - log(geom_mean(Y))
cat(7)
log.lambda.minus.log.gm = ggplot(data.lambda, 
                                aes(x=gamma, y = log.lambda-log.gm)) +
  xlab("Gamma") + ylab("log-Lambda - log-geom_mean(Y)") +
  ggtitle(plot.title) + 
  geom_vline(xintercept = gamma.0) +
  geom_point()


jpeg(paste0("Plots/CV lambda/log-lambda minus log-GM(Y)/JPEG/",
            plot.title, ".jpeg"))
plot(log.lambda.minus.log.gm)
dev.off()

png(paste0("Plots/CV lambda/log-lambda minus log-GM(Y)/PNG/",
           plot.title, ".png"))
plot(log.lambda.minus.log.gm)
dev.off()