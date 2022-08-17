
library(survival)
rep <- 200
qoracle <- qmislin <- qmisint <- vector(mode="list", length=rep)
#################################
##  Test set import            ##
#################################
newdata <- read.csv(file="~/BART/a3/test/data_nln.csv")
newn <- nrow(newdata)
newdata$cosx2_3 <- cos(newdata$x2^3)
newdata$x2b2_2 <- (newdata$x2*newdata$b2+0.5)^2
newdata$sinx1b1 <- sin(pi*newdata$x1*newdata$b1)
newdata$x2_2 <- newdata$x2^2
newdata$sinx1_2 <- sin(newdata$x1^2)
newdata$x1_4 <- newdata$x1^4
newdata$x1b1 <- newdata$x1*newdata$b1
newdata$x1_3 <- newdata$x1^3
for (i in 1:rep){
    print(i)
    data <- read.csv(file=paste0("~/BART/a3/train/data_nln_",i,".csv"))
    n <- nrow(data)
    data$cosx2_3 <- cos(data$x2^3)
    data$x2b2_2 <- (data$x2*data$b2+0.5)^2
    data$sinx1b1 <- sin(pi*data$x1*data$b1)
    data$x2_2 <- data$x2^2
    data$sinx1_2 <- sin(data$x1^2)
    data$x1_4 <- data$x1^4
    data$x1b1 <- data$x1*data$b1
    data$x1_3 <- data$x1^3
###############################
## Stage 2 fitting           ##
###############################
    stg2 <- which(data$eta==1)
    data.stg2 <- data[stg2, ]  #stg2 data
    stg2n <- nrow(data.stg2)
    a2option <- unique(data$a2)  #all possible stg2 actions
    a2n <- length(a2option)  #total number of stg2 actions
    stg2t <- survreg(Surv(y2, delta) ~ cosx2_3 + x2b2_2 + x1 + sinx1b1 + a2 + a2*x2_2, data=data.stg2, dist="lognormal")
    stg2lin <- survreg(Surv(y2, delta) ~ x2 + b2 + z2 + x1 + b1 + z1 + a2, data=data.stg2, dist="lognormal")
    stg2int <- survreg(Surv(y2, delta) ~ x2 + b2 + z2 + x2*b2 + x2*z2 + b2*z2 + x1 + b1 + z1 + x1*b1 + x1*z1 + b1*z1 + a2 + a2*x2 + a2*b2 + a2*z2, data=data.stg2, dist="lognormal")
#################################
## Predict opt stg2 act & time ##
#################################
    x.mat <- rbind(data.stg2,newdata)[,c("x1","x2","b1","b2","z1","z2","cosx2_3","x2b2_2","sinx1b1","x2_2")]  #x matrix for stg2
    x.test.stg2 <- data.frame(x.mat[rep(seq_len(nrow(x.mat)), each=a2n),], a2=rep(a2option,stg2n+newn))
    pred <- function(model, newdata, nobs, ntest, ntrt, trt){
        mu <- predict(model, newdata=newdata)
        aopt <- sapply(1:(nobs+ntest)*ntrt-1, function(column) {which.max(mu[column:(column+ntrt-1)])})  #opt col at stg2
        a.opt <- sapply(aopt, function(x) trt[x])  #opt action at stg2
        yhatoptmean <- sapply(1:(nobs+ntest)*ntrt-1, function(column) {max(mu[column:(column+ntrt-1)])})  #opt y mean at stg2
        return(list(a.opt=a.opt[1:nobs], yhatoptmean=yhatoptmean[1:nobs], newa.opt=a.opt[-(1:nobs)], newyhatoptmean=yhatoptmean[-(1:nobs)]))
    }
    resstg2t <- pred(stg2t, x.test.stg2, stg2n, newn, a2n, a2option)
    resstg2lin <- pred(stg2lin, x.test.stg2, stg2n, newn, a2n, a2option)
    resstg2int <- pred(stg2int, x.test.stg2, stg2n, newn, a2n, a2option)
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
    data2lin <- augdata(data, stg2, resstg2lin$yhatoptmean)
    data2int <- augdata(data, stg2, resstg2int$yhatoptmean)
    stg1t <- survreg(Surv(y1new, delta1) ~ sinx1_2 + x1_4 + x1b1 + a1 + a1*x1_3, data=data2t, dist="lognormal")  #Q1tQ2t
    stg1lin <- survreg(Surv(y1new, delta1) ~ x1 + b1 + z1 + a1, data=data2lin, dist="lognormal")  #main effects
    stg1int <- survreg(Surv(y1new, delta1) ~ x1 + b1 + z1 + x1*b1 +x1*z1 + b1*z1 + a1 + a1*x1 + a1*b1 + a1*z1, data=data2int, dist="lognormal")  #interaction terms
#################################
##  Test set rbind            ##
#################################
    x.mat <- as.matrix(rbind(data[,c("x1","sinx1_2","x1_3","x1_4","x1b1","b1","z1")], newdata[,c("x1","sinx1_2","x1_3","x1_4","x1b1","b1","z1")]))
    x.test.stg1 <- data.frame(x.mat[rep(seq_len(nrow(x.mat)), each=a1n),], a1=rep(a1option,n+newn))  #stg1 actions for everyone
#################################
## Predict opt stg1 act & time ##
#################################
    res1t <- pred(stg1t, x.test.stg1, n, newn, a1n, a1option)
    res1lin <- pred(stg1lin, x.test.stg1, n, newn, a1n, a1option)
    res1int <- pred(stg1int, x.test.stg1, n, newn, a1n, a1option)
#################################
##    Save the results         ##
#################################
    qoracle[[i]] <- list(a2.opt=resstg2t$a.opt, yhat2optmean=resstg2t$yhatoptmean, newa2.opt=resstg2t$newa.opt, newyhat2optmean=resstg2t$newyhatoptmean, a1.opt=res1t$a.opt, yhat1optmean=res1t$yhatoptmean, newa1.opt=res1t$newa.opt, newyhat1optmean=res1t$newyhatoptmean)
    qmislin[[i]] <- list(a2.opt=resstg2lin$a.opt, yhat2optmean=resstg2lin$yhatoptmean, newa2.opt=resstg2lin$newa.opt, newyhat2optmean=resstg2lin$newyhatoptmean, a1.opt=res1lin$a.opt, yhat1optmean=res1lin$yhatoptmean, newa1.opt=res1lin$newa.opt, newyhat1optmean=res1lin$newyhatoptmean)
    qmisint[[i]] <- list(a2.opt=resstg2int$a.opt, yhat2optmean=resstg2int$yhatoptmean, newa2.opt=resstg2int$newa.opt, newyhat2optmean=resstg2int$newyhatoptmean, a1.opt=res1int$a.opt, yhat1optmean=res1int$yhatoptmean, newa1.opt=res1int$newa.opt, newyhat1optmean=res1int$newyhatoptmean)
}

saveRDS(qoracle, file="~/BART/a3/results/qoracle.RDS")
saveRDS(qmislin, file="~/BART/a3/results/qmislin.RDS")
saveRDS(qmisint, file="~/BART/a3/results/qmisint.RDS")
