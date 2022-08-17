library(BART3)
library(parallel)
library(doRNG)
library(doParallel)
source("~/BART/a3/codes/wrapper.R")
rep <- 200
newdata <- read.csv(file="~/BART/a3/test/data_nln.csv")
newdata$eta <- 1
for (i in 1:rep){
    data <- read.csv(file=paste0("~/BART/a3/train/data_nln_",i,".csv"))
    res <- dtr1(x1=c("x1","b1","z1"), a1="a1", time1="y1",
                 x2=c("x1","b1","z1","x2","b2","z2"), a2="a2", stg2ind="eta", time2="y2", delta="delta",
                 data, newdata=newdata, mc.cores=16)
    saveRDS(res, file=paste0("~/BART/a3/results/bart/bart2_",i,".RDS"))
}

