library(doParallel)


#Initialize parallelization
nclust = as.numeric(Sys.getenv("SLURM_NTASKS"))
nclust = ifelse(is.na(nclust), detectCores(), nclust)
cl = makeCluster(nclust)
registerDoParallel(cl)

clusterSetRNGStream(cl=cl, iseed = 123456)

#Pass info to cluster
clusterExport(cl, c("n"))
clusterEvalQ(cl, {
             library(dplyr)
             library(ggplot2)
  })


m = 10
constants = foreach(i = 1:m, .combine="rbind") %dopar% {
                          return(m)

                        }
                        
#Terminate parallelization
stopCluster(cl)





