library(doParallel)

n = 10

#Initialize parallelization
nclust = as.numeric(Sys.getenv("SLURM_NTASKS"))
nclust = ifelse(is.na(nclust), detectCores(), nclust)
cl = makeCluster(nclust)
registerDoParallel(cl)

### This isn't really necessary
### You can just set the seed in the function passed to parSapply
clusterSetRNGStream(cl=cl, iseed = 123456)

#Pass info to cluster
clusterExport(cl, c("n"))
clusterEvalQ(cl, {
             library(dplyr)
             library(ggplot2)
  })


M = 100

num = parSapply(1:M, function(j) {
  # set.seed(1)
  rnorm(n)
}, cl = cl)

constants = foreach(i = 1:M, .combine="rbind") %dopar% {
                          return(m)

                        }
                        
#Terminate parallelization
stopCluster(cl)





