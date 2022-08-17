
library(survival)
rep <- 200
q1t2t <- q1f2t <- q1t2f <- q1f2f <- qmis2 <- qmis3 <- qmis4 <- qmis5 <- vector(mode="list", length=rep)
#################################
##  Test set import            ##
#################################
newdata <- read.csv(file="~/BART/a3/test/data_lin.csv")
newn <- nrow(newdata)
#newdata$x2_3 <- newdata$x2^3
#newdata$x1_4 <- newdata$x1^4
#newdata$x2_4 <- newdata$x2^4
#newdata$x1_2 <- newdata$x1^2
for (i in 1:rep){
    print(i)
    data <- read.csv(file=paste0("~/BART/a3/train/data_lin_",i,".csv"))
    n <- nrow(data)
#    data$x2_3 <- data$x2^3
#    data$x1_4 <- data$x1^4
#    data$x2_4 <- data$x2^4
#    data$x1_2 <- data$x1^2
###############################
## Stage 2 fitting           ##
###############################
    stg2 <- which(data$eta==1)
    data.stg2 <- data[stg2, ]  #stg2 data
    stg2n <- nrow(data.stg2)
    a2option <- unique(data$a2)  #all possible stg2 actions
    a2n <- length(a2option)  #total number of stg2 actions
    stg2t <- survreg(Surv(y2, delta) ~ x2 + x2*b2 + b2 + x1 + b1 + x1*b1 + a2 + a2*x2 + a2*b2, data=data.stg2, dist="lognormal")
#    q1stg2 <- survreg(Surv(y2, delta) ~ x2+b2+x2*b2+x1+b1+x1*b1+a2+a2*x2+a2*b2, data=data.stg2, dist="lognormal")
    stg2f <- survreg(Surv(y2, delta) ~ x2 + b2 + z2 + x1 + b1 + a2 + a2*x2 + a2*z2, data=data.stg2, dist="lognormal")
#    qmis2stg2 <- survreg(Surv(y2, delta) ~ x1 + a2 + a2*x1, data=data.stg2, dist="lognormal")
#    qmis3stg2 <- survreg(Surv(y2, delta) ~ x2+a2, data=data.stg2, dist="lognormal")
#    qmis4stg2 <- survreg(Surv(y2, delta) ~ x2_4+x1+a2+a2*x1, data=data.stg2, dist="lognormal")
#################################
## Predict opt stg2 act & time ##
#################################
    x.mat <- rbind(data.stg2[,c("x1","x2","b1","b2","z1","z2")],newdata[,c("x1","x2","b1","b2","z1","z2")])  #x matrix for stg2
    x.test.stg2 <- data.frame(x.mat[rep(seq_len(nrow(x.mat)), each=a2n),], a2=rep(a2option,stg2n+newn))
    pred <- function(model, newdata, nobs, ntest, ntrt, trt){
        mu <- predict(model, newdata=newdata)
        aopt <- sapply(1:(nobs+ntest)*ntrt-1, function(column) {which.max(mu[column:(column+ntrt-1)])})  #opt col at stg2
        a.opt <- sapply(aopt, function(x) trt[x])  #opt action at stg2
        yhatoptmean <- sapply(1:(nobs+ntest)*ntrt-1, function(column) {max(mu[column:(column+ntrt-1)])})  #opt y mean at stg2
        return(list(a.opt=a.opt[1:nobs], yhatoptmean=yhatoptmean[1:nobs], newa.opt=a.opt[-(1:nobs)], newyhatoptmean=yhatoptmean[-(1:nobs)]))
    }
    resstg2t <- pred(stg2t, x.test.stg2, stg2n, newn, a2n, a2option)
    resstg2f <- pred(stg2f, x.test.stg2, stg2n, newn, a2n, a2option)
#    resmis2stg2 <- pred(qmis2stg2, x.test.stg2, stg2n, newn, a2n, a2option)
#    resmis3stg2 <- pred(qmis3stg2, x.test.stg2, stg2n, newn, a2n, a2option)
#    resmis4stg2 <- pred(qmis4stg2, x.test.stg2, stg2n, newn, a2n, a2option)
############################
## Stage 1 fitting        ##
############################
    data$delta1 <- 1-(1-data$delta)*(1-data$eta)
    a1option <- unique(data$a1)  #all possible stage1 actions
    a1n <- length(a1option)  #total number of stage1 actions
    augdata <- function(data, stg2, yhat2){
        data$y2opt <- 0
        data[stg2,]$y2opt <- yhat2
        data$y1new <- data$y1+data$y2opt
        return(data)
    }
    data2t <- augdata(data, stg2, resstg2t$yhatoptmean)
    data2f <- augdata(data, stg2, resstg2f$yhatoptmean)
#    data2mis2 <- augdata(data, stg2, resmis2stg2$yhatoptmean)
#    data2mis3 <- augdata(data, stg2, resmis3stg2$yhatoptmean)
#    data2mis4 <- augdata(data, stg2, resmis4stg2$yhatoptmean)
    stg1t2t <- survreg(Surv(y1new, delta1) ~ x1 + b1 + x1*b1 + a1 + a1*x1 + a1*b1, data=data2t, dist="lognormal")  #Q1tQ2t
    stg1f2t <- survreg(Surv(y1new, delta1) ~ x1 + b1 + z1 + a1 + a1*x1 + a1*z1, data=data2t, dist="lognormal")  #Q1fQ2t
    stg1t2f <- survreg(Surv(y1new, delta1) ~ x1 + b1 + x1*b1 + a1 + a1*x1 + a1*b1, data=data2f, dist="lognormal")  #Q1tQ2f
    stg1f2f <- survreg(Surv(y1new, delta1) ~ x1 + b1 + z1 + a1 + a1*x1 + a1*z1, data=data2f, dist="lognormal")  #Q1fQ2f
#    mis2stg1 <- survreg(Surv(y1new, delta1) ~ x1_2 + a1 + a1*x1, data=data2mis2, dist="lognormal")
#    mis3stg1 <- survreg(Surv(y1new, delta1) ~ x1 + x1_2 + a1 + a1*x1, data=data2mis3, dist="lognormal")
#    mis4stg1 <- survreg(Surv(y1new, delta1) ~ x1 + a1 + a1*x1, data=data2mis4, dist="lognormal")
#    mis5stg1 <- survreg(Surv(y1new, delta1) ~ x1 + x1_4 + a1, data=data2t, dist="lognormal")
#################################
##  Test set rbind            ##
#################################
    x.mat <- as.matrix(rbind(data[,c("x1","b1","z1")],newdata[,c("x1","b1","z1")]))
    x.test.stg1 <- data.frame(x.mat[rep(seq_len(nrow(x.mat)), each=a1n),], a1=rep(a1option,n+newn))  #stg1 actions for everyone
#################################
## Predict opt stg1 act & time ##
#################################
    res1t2t <- pred(stg1t2t, x.test.stg1, n, newn, a1n, a1option)
    res1f2t <- pred(stg1f2t, x.test.stg1, n, newn, a1n, a1option)
    res1t2f <- pred(stg1t2f, x.test.stg1, n, newn, a1n, a1option)
    res1f2f <- pred(stg1f2f, x.test.stg1, n, newn, a1n, a1option)
#    resmis2stg1 <- pred(mis2stg1, x.test.stg1, n, newn, a1n, a1option)
#    resmis3stg1 <- pred(mis3stg1, x.test.stg1, n, newn, a1n, a1option)
#    resmis4stg1 <- pred(mis4stg1, x.test.stg1, n, newn, a1n, a1option)
#    resmis5stg1 <- pred(mis5stg1, x.test.stg1, n, newn, a1n, a1option)
#################################
##    Save the results         ##
#################################
    q1t2t[[i]] <- list(a2.opt=resstg2t$a.opt, yhat2optmean=resstg2t$yhatoptmean, newa2.opt=resstg2t$newa.opt, newyhat2optmean=resstg2t$newyhatoptmean, a1.opt=res1t2t$a.opt, yhat1optmean=res1t2t$yhatoptmean, newa1.opt=res1t2t$newa.opt, newyhat1optmean=res1t2t$newyhatoptmean)
    q1f2t[[i]] <- list(a2.opt=resstg2t$a.opt, yhat2optmean=resstg2t$yhatoptmean, newa2.opt=resstg2t$newa.opt, newyhat2optmean=resstg2t$newyhatoptmean, a1.opt=res1f2t$a.opt, yhat1optmean=res1f2t$yhatoptmean, newa1.opt=res1f2t$newa.opt, newyhat1optmean=res1f2t$newyhatoptmean)
    q1t2f[[i]] <- list(a2.opt=resstg2f$a.opt, yhat2optmean=resstg2f$yhatoptmean, newa2.opt=resstg2f$newa.opt, newyhat2optmean=resstg2f$newyhatoptmean, a1.opt=res1t2f$a.opt, yhat1optmean=res1t2f$yhatoptmean, newa1.opt=res1t2f$newa.opt, newyhat1optmean=res1t2f$newyhatoptmean)
    q1f2f[[i]] <- list(a2.opt=resstg2f$a.opt, yhat2optmean=resstg2f$yhatoptmean, newa2.opt=resstg2f$newa.opt, newyhat2optmean=resstg2f$newyhatoptmean, a1.opt=res1f2f$a.opt, yhat1optmean=res1f2f$yhatoptmean, newa1.opt=res1f2f$newa.opt, newyhat1optmean=res1f2f$newyhatoptmean)
#    qmis2[[i]] <- list(a2.opt=resmis2stg2$a.opt, yhat2optmean=resmis2stg2$yhatoptmean, newa2.opt=resmis2stg2$newa.opt, newyhat2optmean=resmis2stg2$newyhatoptmean, a1.opt=resmis2stg1$a.opt, yhat1optmean=resmis2stg1$yhatoptmean, newa1.opt=resmis2stg1$newa.opt, newyhat1optmean=resmis2stg1$newyhatoptmean)
#    qmis3[[i]] <- list(a2.opt=resmis3stg2$a.opt, yhat2optmean=resmis3stg2$yhatoptmean, newa2.opt=resmis3stg2$newa.opt, newyhat2optmean=resmis3stg2$newyhatoptmean, a1.opt=resmis3stg1$a.opt, yhat1optmean=resmis3stg1$yhatoptmean, newa1.opt=resmis3stg1$newa.opt, newyhat1optmean=resmis3stg1$newyhatoptmean)
#    qmis4[[i]] <- list(a2.opt=resmis4stg2$a.opt, yhat2optmean=resmis4stg2$yhatoptmean, newa2.opt=resmis4stg2$newa.opt, newyhat2optmean=resmis4stg2$newyhatoptmean, a1.opt=resmis4stg1$a.opt, yhat1optmean=resmis4stg1$yhatoptmean, newa1.opt=resmis4stg1$newa.opt, newyhat1optmean=resmis4stg1$newyhatoptmean)
#    qmis5[[i]] <- list(a2.opt=resstg2t$a.opt, yhat2optmean=resstg2t$yhatoptmean, newa2.opt=resstg2t$newa.opt, newyhat2optmean=resstg2t$newyhatoptmean, a1.opt=resmis5stg1$a.opt, yhat1optmean=resmis5stg1$yhatoptmean, newa1.opt=resmis5stg1$newa.opt, newyhat1optmean=resmis5stg1$newyhatoptmean)
}

saveRDS(q1t2t, file="~/BART/a3/results/q1t2t.RDS")
saveRDS(q1f2t, file="~/BART/a3/results/q1f2t.RDS")
saveRDS(q1t2f, file="~/BART/a3/results/q1t2f.RDS")
saveRDS(q1f2f, file="~/BART/a3/results/q1f2f.RDS")
#saveRDS(qmis2, file="~/BART/a3/results/qmis2.RDS")
#saveRDS(qmis3, file="~/BART/a3/results/qmis3.RDS")
#saveRDS(qmis4, file="~/BART/a3/results/qmis4.RDS")
#saveRDS(qmis5, file="~/BART/a3/results/qmis5.RDS")
