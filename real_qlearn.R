realdata <- read.csv(file="~/BART/a3/data/all.comp.cases.csv")
library(survival)
library(ggplot2)
varnames <- names(realdata)
varlist <- varnames[c(8,9,11,13,14,15,16)]
realdata[,varlist] <- lapply(realdata[,varlist], as.factor)

x1 <- varnames[6:16]  # covariates in Stage 1
x2 <- varnames[6:18]  # covariates in Stage 2

set.seed(73022)
nboot <- 1000
seedset <- runif(nboot)*100000000
################
## Q-learning
################
n <- nrow(realdata)
real.stg2 <- realdata[realdata$Int2==1,]
real2n <- nrow(real.stg2)
## Stg 2 regression
a2option <- unique(realdata$A2.NHTL)
a2n <- length(a2option)
qcov2 <- paste(paste("A2.NHTL*",c(x2,"A1.NHTL")),collapse="+")
fml2 <- as.formula(paste("Surv(Y2,delta2)", qcov2, sep="~"))
qstg2 <- survreg(fml2, data=real.stg2, dist="lognormal")
stg2.est <- summary(qstg2)$table[,1]
x.mat <- real.stg2[,c(x2,"A1.NHTL")]
x.test.stg2 <- data.frame(x.mat[rep(seq_len(nrow(x.mat)), each=a2n),], A2.NHTL=rep(a2option,real2n))
pred.stg2 <- predict(qstg2, newdata=x.test.stg2)
a2names <- paste0("A2.NHTL_",a2option)
for (i in 1:a2n) assign(a2names[i], pred.stg2[1:real2n*a2n-(a2n-i)])
a2opt <- sapply(1:real2n*a2n-1, function(column) {which.max(pred.stg2[column:(column+a2n-1)])})  #opt col at stg2
a2.opt <- sapply(a2opt, function(x) a2option[x])  #opt action at stg2
yhat2optmean <- sapply(1:real2n*a2n-1, function(column) {max(pred.stg2[column:(column+a2n-1)])})
## Stg 1 regression
a1option <- unique(realdata$A1.NHTL)  #all possible stage1 actions
a1n <- length(a1option)  #total number of stage1 actions
data.stg1 <- realdata
data.stg1$y2opt <- 0
data.stg1[data.stg1$Int2==1,]$y2opt <- yhat2optmean
data.stg1$Y1 <- data.stg1$Y1.tilde+data.stg1$y2opt
qcov1 <- paste(paste("A1.NHTL*",x1),collapse="+")
fml1 <- as.formula(paste("Surv(Y1,delta1)", qcov1, sep="~"))
qstg1 <- survreg(fml1, data=data.stg1, dist="lognormal")
stg1.est <- summary(qstg1)$table[,1]
x.mat <- data.stg1[,x1]
x.test.stg1 <- data.frame(x.mat[rep(seq_len(nrow(x.mat)), each=a1n),], A1.NHTL=rep(a1option,n))
pred.stg1 <- predict(qstg1, newdata=x.test.stg1)
a1names <- paste0("A1.NHTL_",a1option)
for (i in 1:a1n) assign(a1names[i], pred.stg1[1:n*a1n-(a1n-i)])
a1opt <- sapply(1:n*a1n-1, function(column) {which.max(pred.stg1[column:(column+a1n-1)])})  #opt col at stg1
a1.opt <- sapply(a1opt, function(x) a1option[x])  #opt action at stg1
yhat1optmean <- sapply(1:n*a1n-1, function(column) {max(pred.stg1[column:(column+a1n-1)])})
orig <- list(a2.opt=a2.opt, yhat2optmean=yhat2optmean, a1.opt=a1.opt, yhat1optmean=yhat1optmean, A2.NHTL_0=A2.NHTL_0, A2.NHTL_1=A2.NHTL_1, A1.NHTL_0=A1.NHTL_0, A1.NHTL_1=A1.NHTL_1)
saveRDS(orig, file="~/BART/a3/results/real_q.RDS")

A1.NHTL0 <- A1.NHTL1 <- A2.NHTL0 <- A2.NHTL1 <- stg2boot <- stg1boot <- sd1 <- sd2 <- NULL
for (i in 1:nboot){
    set.seed(seedset[i])
    if (i %% 100 == 0) print(i)
    bootdata <- realdata[sample(n, n, replace=TRUE),]
    ## Stg 2 regression
    boot.stg2 <- bootdata[bootdata$Int2==1,]
    stg2n <- nrow(boot.stg2)
    qstg2 <- survreg(fml2, data=boot.stg2, dist="lognormal")
    sd2 <- c(sd2, qstg2$scale)
    stg2boot <- rbind(stg2boot, summary(qstg2)$table[,1])
    x.mat <- rbind(boot.stg2,real.stg2)[,c(x2,"A1.NHTL")]
    x.test.stg2 <- data.frame(x.mat[rep(seq_len(nrow(x.mat)), each=a2n),], A2.NHTL=rep(a2option,stg2n+real2n))
    pred.stg2 <- predict(qstg2, newdata=x.test.stg2)
    for (i in 1:a2n) assign(a2names[i], pred.stg2[(stg2n+1:real2n)*a2n-(a2n-i)])
    a2opt <- sapply(1:(stg2n+real2n)*a2n-1, function(column) {which.max(pred.stg2[column:(column+a2n-1)])})  #opt col at stg2
    a2.opt <- sapply(a2opt, function(x) a2option[x])  #opt action at stg2
    yhat2optmean <- sapply(1:(stg2n+real2n)*a2n-1, function(column) {max(pred.stg2[column:(column+a2n-1)])})
    A2.NHTL0 <- cbind(A2.NHTL0, A2.NHTL_0)
    A2.NHTL1 <- cbind(A2.NHTL1, A2.NHTL_1)
    ## Stg 1 regression
    bootdata$y2opt <- 0
    bootdata[bootdata$Int2==1,]$y2opt <- yhat2optmean[stg2n]
    bootdata$Y1 <- bootdata$Y1.tilde+bootdata$y2opt
    qstg1 <- survreg(fml1, data=bootdata, dist="lognormal")
    sd1 <- c(sd1, qstg1$scale)
    stg1boot <- rbind(stg1boot, summary(qstg1)$table[,1])
    x.mat <- rbind(bootdata[,x1], realdata[,x1])
    x.test.stg1 <- data.frame(x.mat[rep(seq_len(nrow(x.mat)), each=a1n),], A1.NHTL=rep(a1option,n+n))
    pred.stg1 <- predict(qstg1, newdata=x.test.stg1)
    for (i in 1:a1n) assign(a1names[i], pred.stg1[(n+1:n)*a1n-(a1n-i)])
    a1opt <- sapply(1:(n+n)*a1n-1, function(column) {which.max(pred.stg1[column:(column+a1n-1)])})  #opt col at stg1
    a1.opt <- sapply(a1opt, function(x) a1option[x])  #opt action at stg1
    yhat1optmean <- sapply(1:(n+n)*a1n-1, function(column) {max(pred.stg1[column:(column+a1n-1)])})
    A1.NHTL0 <- cbind(A1.NHTL0, A1.NHTL_0)
    A1.NHTL1 <- cbind(A1.NHTL1, A1.NHTL_1)
}

stg2sd <- apply(stg2boot, 2, sd)
plot(colMeans(stg2boot), summary(qstg2)$table[,1])
abline(a=0, b=1, col=2)
plot(stg2sd, summary(qstg2)$table[,2])
abline(a=0, b=1, col=2)
stg2lq <- stg2.est + qnorm(0.025)*stg2sd
stg2uq <- stg2.est + qnorm(0.975)*stg2sd
stg2ci <- apply(stg2boot, 2, function(x) quantile(x, probs=c(0.025,0.975)))
par(mfrow=c(1,2))
plot(stg2ci[1,], stg2lq)
abline(a=0, b=1, col=2)
plot(stg2ci[2,], stg2uq)
abline(a=0, b=1, col=2)

stg1sd <- apply(stg1boot, 2, sd)
plot(colMeans(stg1boot), summary(qstg1)$table[,1])
abline(a=0, b=1, col=2)
plot(stg1sd, summary(qstg1)$table[,2])
abline(a=0, b=1, col=2)
stg1lq <- stg1.est + qnorm(0.025)*stg1sd
stg1uq <- stg1.est + qnorm(0.975)*stg1sd
stg1ci <- apply(stg1boot, 2, function(x) quantile(x, probs=c(0.025,0.975)))
par(mfrow=c(1,2))
plot(stg1ci[1,], stg1lq)
abline(a=0, b=1, col=2)
plot(stg1ci[1,], stg1uq)
abline(a=0, b=1, col=2)

stg2.est <- cbind(stg2.est, colMeans(stg2boot), t(stg2ci))
write.csv(stg2.est, file="~/BART/a3/outputs/stg2est.csv")
stg1.est <- cbind(stg1.est, colMeans(stg1boot), t(stg1ci))
write.csv(stg1.est, file="~/BART/a3/outputs/stg1est.csv")

yhatdiff1 <- A1.NHTL1-A1.NHTL0
sigma1 <- apply(yhatdiff1, 1, sd)
origdiff1 <- orig$A1.NHTL_1-orig$A1.NHTL_0
lowerq <- origdiff1 + qnorm(0.025)*sigma1
upperq <- origdiff1 + qnorm(0.975)*sigma1
q1 <- origdiff1 + qnorm(0.25)*sigma1
q3 <- origdiff1 + qnorm(0.75)*sigma1

output <- data.frame(stg1_diff=origdiff1, lowerq, upperq, q1, q3)
sortout <- output[order(output$stg1_diff, decreasing=TRUE),]
sortout$x <- 1:nrow(sortout)

b <- ggplot(sortout, aes(x=x, y=stg1_diff)) + theme(axis.text=element_text(size=16,face="bold"),axis.title=element_text(size=18,face="bold"),title=element_text(size=18,face="bold")) + labs(list(x = "Patient Number \n(Ordered by risk difference)", y = "Risk difference")) 
smoother <- geom_line(color="black",size=2)
quartiles <- geom_pointrange(data=sortout, aes(x=x, y=stg1_diff, ymax=q3, ymin=q1), colour="gray30")
p1a <- b+ geom_pointrange(data=sortout, aes(x=x, y=stg1_diff, ymax=upperq,ymin=lowerq),colour="gray50")+quartiles+smoother+scale_y_continuous("Bootstrap median difference in predicted DFS time")+scale_x_continuous("Patient number (Ordered)")+ggtitle("Stage 1")
## +coord_cartesian(ylim=c(-0.25,0.25))
## ,breaks=c(-0.2,-0.1,0,0.1,0.2)
pdf(file="~/BART/a3/outputs/real_stg1dfs_q.pdf", width=8,height=8)
p1a
dev.off()


yhatdiff2 <- A2.NHTL1-A2.NHTL0
sigma2 <- apply(yhatdiff2, 1, sd)
origdiff2 <- orig$A2.NHTL_1-orig$A2.NHTL_0
lowerq <- origdiff2 + qnorm(0.025)*sigma2
upperq <- origdiff2 + qnorm(0.975)*sigma2
q1 <- origdiff2 + qnorm(0.25)*sigma2
q3 <- origdiff2 + qnorm(0.75)*sigma2

output <- data.frame(stg2_diff=origdiff2, lowerq, upperq, q1, q3)
sortout <- output[order(output$stg2_diff, decreasing=TRUE),]
sortout$x <- 1:nrow(sortout)

b <- ggplot(sortout, aes(x=x, y=stg2_diff)) + theme(axis.text=element_text(size=16,face="bold"),axis.title=element_text(size=18,face="bold"),title=element_text(size=18,face="bold")) + labs(list(x = "Patient Number \n(Ordered by risk difference)", y = "Risk difference")) 
smoother <- geom_line(color="black",size=2)
quartiles <- geom_pointrange(data=sortout, aes(x=x, y=stg2_diff, ymax=q3, ymin=q1), colour="gray30")
p2a <- b+ geom_pointrange(data=sortout, aes(x=x, y=stg2_diff, ymax=upperq,ymin=lowerq),colour="gray50")+quartiles+smoother+scale_y_continuous("Bootstrap median difference in predicted DFS time")+scale_x_continuous("Patient number (Ordered)")+ggtitle("Stage 2")
## +coord_cartesian(ylim=c(-0.25,0.25))
## ,breaks=c(-0.2,-0.1,0,0.1,0.2)
pdf(file="~/BART/a3/outputs/real_stg2dfs_q.pdf", width=8,height=8)
p2a
dev.off()

yhatdiff1.lg <- log(A1.NHTL1)-log(A1.NHTL0)
sigma1.lg <- apply(yhatdiff1.lg, 1, sd)
origdiff1.lg <- log(orig$A1.NHTL_1)-log(orig$A1.NHTL_0)
lowerq <- origdiff1.lg + qnorm(0.025)*sigma1.lg
upperq <- origdiff1.lg + qnorm(0.975)*sigma1.lg
q1 <- origdiff1.lg + qnorm(0.25)*sigma1.lg
q3 <- origdiff1.lg + qnorm(0.75)*sigma1.lg

output <- data.frame(stg1_diff=log(orig$A1.NHTL_1)-log(orig$A1.NHTL_0), lowerq, upperq, q1, q3)
sortout <- output[order(output$stg1_diff, decreasing=TRUE),]
sortout$x <- 1:nrow(sortout)

b <- ggplot(sortout, aes(x=x, y=stg1_diff)) + theme(axis.text=element_text(size=16,face="bold"),axis.title=element_text(size=18,face="bold"),title=element_text(size=18,face="bold")) + labs(list(x = "Patient Number \n(Ordered by risk difference)", y = "Risk difference")) 
smoother <- geom_line(color="black",size=2)
quartiles <- geom_pointrange(data=sortout, aes(x=x, y=stg1_diff, ymax=q3, ymin=q1), colour="gray30")
p1a <- b+ geom_pointrange(data=sortout, aes(x=x, y=stg1_diff, ymax=upperq,ymin=lowerq),colour="gray50")+quartiles+smoother+scale_y_continuous("Bootstrap mean difference in predicted log(DFS time)")+scale_x_continuous("Patient number (Ordered)")+ggtitle("Stage 1")
pdf(file="~/BART/a3/outputs/real_stg1dfs_qlg.pdf", width=8,height=8)
p1a
dev.off()

yhatdiff2.lg <- log(A2.NHTL1)-log(A2.NHTL0)
sigma2.lg <- apply(yhatdiff2.lg, 1, sd)
origdiff2.lg <- log(orig$A2.NHTL_1)-log(orig$A2.NHTL_0)
lowerq <- origdiff2.lg + qnorm(0.025)*sigma2.lg
upperq <- origdiff2.lg + qnorm(0.975)*sigma2.lg
q1 <- origdiff2.lg + qnorm(0.25)*sigma2.lg
q3 <- origdiff2.lg + qnorm(0.75)*sigma2.lg

output <- data.frame(stg2_diff=log(orig$A2.NHTL_1)-log(orig$A2.NHTL_0), lowerq, upperq, q1, q3)
sortout <- output[order(output$stg2_diff, decreasing=TRUE),]
sortout$x <- 1:nrow(sortout)

b <- ggplot(sortout, aes(x=x, y=stg2_diff)) + theme(axis.text=element_text(size=16,face="bold"),axis.title=element_text(size=18,face="bold"),title=element_text(size=18,face="bold")) + labs(list(x = "Patient Number \n(Ordered by risk difference)", y = "Risk difference")) 
smoother <- geom_line(color="black",size=2)
quartiles <- geom_pointrange(data=sortout, aes(x=x, y=stg2_diff, ymax=q3, ymin=q1), colour="gray30")
p2a <- b+ geom_pointrange(data=sortout, aes(x=x, y=stg2_diff, ymax=upperq,ymin=lowerq),colour="gray50")+quartiles+smoother+scale_y_continuous("Bootstrap mean difference in predicted log(DFS time)")+scale_x_continuous("Patient number (Ordered)")+ggtitle("Stage 2")
pdf(file="~/BART/a3/outputs/real_stg2dfs_qlg.pdf", width=8,height=8)
p2a
dev.off()

t <- log(2*12)  #original Ys are months
## survival probability for each treatment
s1 <- s0 <- NULL
for (i in 1:nboot){
    s1 <- cbind(s1, pnorm(q=t, mean=log(A1.NHTL1)[,i], sd=sd1[i], lower.tail=FALSE))
    s0 <- cbind(s0, pnorm(q=t, mean=log(A1.NHTL0)[,i], sd=sd1[i], lower.tail=FALSE))
}
sdiff1 <- s1-s0
sdiff1mean <- rowMeans(sdiff1)
sdiff1ci <- apply(sdiff1, 1, quantile, probs=c(0.025,0.975,0.25,0.75))
sdiff1sum <- data.frame(mean=sdiff1mean, lowerq=sdiff1ci[1,], upperq=sdiff1ci[2,], q1=sdiff1ci[3,], q3=sdiff1ci[4,])
sdiff1sum <- sdiff1sum[order(sdiff1sum$mean, decreasing=TRUE),]
sdiff1sum$x <- 1:nrow(sdiff1sum)

b <- ggplot(sdiff1sum, aes(x=x, y=mean)) + theme(axis.text=element_text(size=16,face="bold"),axis.title=element_text(size=18,face="bold"),title=element_text(size=18,face="bold"))
smoother <- geom_line(color="black",size=2)
quartiles <- geom_pointrange(data=sdiff1sum, aes(x=x, y=mean, ymax=q3, ymin=q1), colour="gray30")
p2a <- b+ geom_pointrange(data=sdiff1sum, aes(x=x, y=mean, ymax=upperq,ymin=lowerq),colour="gray50")+quartiles+smoother+scale_y_continuous("Difference in 2-yr survival probability")+scale_x_continuous("Patient number (Ordered)")+ggtitle("Stage 1")
pdf(file="~/BART/a3/outputs/real_stg1_2yr_q.pdf", width=8,height=8)
p2a
dev.off()

s2 <- s0 <- NULL
for (i in 1:nboot){
    s2 <- cbind(s2, pnorm(q=t, mean=log(A2.NHTL1)[,i], sd=sd2[i], lower.tail=FALSE))
    s0 <- cbind(s0, pnorm(q=t, mean=log(A2.NHTL0)[,i], sd=sd2[i], lower.tail=FALSE))
}
sdiff2 <- s2-s0
sdiff2mean <- rowMeans(sdiff2)
sdiff2ci <- apply(sdiff2, 1, quantile, probs=c(0.025,0.975,0.25,0.75))
sdiff2sum <- data.frame(mean=sdiff2mean, lowerq=sdiff2ci[1,], upperq=sdiff2ci[2,], q1=sdiff2ci[3,], q3=sdiff2ci[4,])
sdiff2sum <- sdiff2sum[order(sdiff2sum$mean, decreasing=TRUE),]
sdiff2sum$x <- 1:nrow(sdiff2sum)

b <- ggplot(sdiff2sum, aes(x=x, y=mean)) + theme(axis.text=element_text(size=16,face="bold"),axis.title=element_text(size=18,face="bold"),title=element_text(size=18,face="bold"))
smoother <- geom_line(color="black",size=2)
quartiles <- geom_pointrange(data=sdiff2sum, aes(x=x, y=mean, ymax=q3, ymin=q1), colour="gray30")
p2a <- b+ geom_pointrange(data=sdiff2sum, aes(x=x, y=mean, ymax=upperq,ymin=lowerq),colour="gray50")+quartiles+smoother+scale_y_continuous("Difference in 2-yr survival probability")+scale_x_continuous("Patient number (Ordered)")+ggtitle("Stage 2")
pdf(file="~/BART/a3/outputs/real_stg2_2yr_q.pdf", width=8,height=8)
p2a
dev.off()
