
library(BART3)
library(parallel)
library(doRNG)
library(doParallel)
source("~/BART/a3/codes/wrapper.R")
rep <- 20
bart2 <- vector(mode="list", length=rep)
for (i in 1:1){
    tdata <- read.csv(file=paste0("~/BART/a3/data/training/data_lin_",i,".csv"))
    res <- dtr3(tdata, 16)
    bart2[[i]] <- res
}
saveRDS(bart2, file="~/BART/a3/results/bart2.RDS")
