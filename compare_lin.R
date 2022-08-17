
library(ggplot2)
library(ggpubr)
q1t2t <- readRDS(file="~/BART/a3/results/q1t2t.RDS")
q1f2t <- readRDS(file="~/BART/a3/results/q1f2t.RDS")
q1t2f <- readRDS(file="~/BART/a3/results/q1t2f.RDS")
q1f2f <- readRDS(file="~/BART/a3/results/q1f2f.RDS")
## qmis2 <- readRDS(file="~/BART/a3/results/qmis2.RDS")
## qmis3 <- readRDS(file="~/BART/a3/results/qmis3.RDS")
## qmis4 <- readRDS(file="~/BART/a3/results/qmis4.RDS")
## qmis5 <- readRDS(file="~/BART/a3/results/qmis5.RDS")
newdata <- read.csv(file="~/BART/a3/test/data_lin.csv")
a2true <- as.numeric(-0.7+0.5*newdata$x2-0.9*newdata$b2>0)
t2opt <- log(newdata$t2opt_true)
a1true <- as.numeric(0.1-0.2*newdata$x1+0.6*newdata$b1>0)
t <- log(newdata$t_true)
rep <- length(q1t2t)
a2mat <- matrix(replicate(rep,a2true), nrow=rep, byrow=TRUE)
a1mat <- matrix(replicate(rep,a1true), nrow=rep, byrow=TRUE)
t2mat <- matrix(replicate(rep,t2opt), nrow=rep, byrow=TRUE)
tmat <- matrix(replicate(rep,t), nrow=rep, byrow=TRUE)

sumr <- function(rep, res){
    a2comb <- yhat2comb <- a1comb <- yhat1comb <- NULL
    for (i in 1:rep){
        a2comb <- rbind(a2comb, res[[i]]$newa2.opt)
        yhat2comb <- rbind(yhat2comb, res[[i]]$newyhat2optmean)
        a1comb <- rbind(a1comb, res[[i]]$newa1.opt)
        yhat1comb <- rbind(yhat1comb, res[[i]]$newyhat1optmean)
    }
    return(list(a2comb=a2comb, yhat2comb=yhat2comb, a1comb=a1comb, yhat1comb=yhat1comb))
}

res.1t2t <- sumr(rep, q1t2t)
res.1f2t <- sumr(rep, q1f2t)
res.1t2f <- sumr(rep, q1t2f)
res.1f2f <- sumr(rep, q1f2f)
## res.mis2 <- sumr(rep, qmis2)
## res.mis3 <- sumr(rep, qmis3)
## res.mis4 <- sumr(rep, qmis4)
## res.mis5 <- sumr(rep, qmis5)

a2.b <- y2.b <- cr2.b <- a1.b <- y1.b <- cr1.b <- NULL
for (i in 1:rep){
    bart1 <- readRDS(file=paste0("~/BART/a3/results/bart/bart1_",i,".RDS"))
    a2.b <- rbind(a2.b, as.numeric(apply(bart1$newa2.opt, 1, function(x) names(sort(table(x), decreasing=TRUE)[1]))))
    y2.ci <- apply(bart1$newyhat2optmean, 1, function(x) quantile(x, probs=c(0.025, 0.975)))
    y2.b <- rbind(y2.b, rowMeans(bart1$newyhat2optmean))
    cr2.b <- rbind(cr2.b, y2.ci[1,]<t2opt & y2.ci[2,]>t2opt)
    a1.b <- rbind(a1.b, as.numeric(apply(bart1$newa1.opt, 1, function(x) names(sort(table(x), decreasing=TRUE)[1]))))
    y1.ci <- apply(bart1$newyhat1optmean, 1, function(x) quantile(x, probs=c(0.025, 0.975)))
    y1.b <- rbind(y1.b, rowMeans(bart1$newyhat1optmean))
    cr1.b <- rbind(cr1.b, y1.ci[1,]<t & y1.ci[2,]>t)
}
res.b <- list(a2comb=a2.b, yhat2comb=y2.b, a1comb=a1.b, yhat1comb=y1.b, cr2comb=cr2.b, cr1comb=cr1.b)

crdata <- data.frame(grp=rep(c("Stage 2","Stage 1"),each=ncol(res.b$cr2comb)), cr=c(colMeans(res.b$cr2comb),colMeans(res.b$cr1comb)))
jpeg(file="~/BART/a3/outputs/lin_cr.jpg", width=480, height=480)
boxplot(cr~grp, data=crdata, boxwex=0.4, pch=20, col=c("white","gray"), frame.plot=TRUE, xlab="", ylab="95% Coverage Rate", ylim=c(0,1))
abline(h=0.95, lty=2, col=2, lwd=2)
dev.off()

pot.1t2t <- data.frame(stg2=mean(res.1t2t$a2comb==a2mat), stg1=mean(res.1t2t$a1comb==a1mat))
bias.1t2ts2 <- log(res.1t2t$yhat2comb)-t2mat
bias.1t2ts1 <- log(res.1t2t$yhat1comb)-tmat
mse.1t2t <- data.frame(stg2=colMeans(bias.1t2ts2^2), stg1=colMeans(bias.1t2ts1^2))
pot.1f2t <- data.frame(stg2=mean(res.1f2t$a2comb==a2mat), stg1=mean(res.1f2t$a1comb==a1mat))
bias.1f2ts2 <- log(res.1f2t$yhat2comb)-t2mat
bias.1f2ts1 <- log(res.1f2t$yhat1comb)-tmat
mse.1f2t <- data.frame(stg2=colMeans(bias.1f2ts2^2), stg1=colMeans(bias.1f2ts1^2))
pot.1t2f <- data.frame(stg2=mean(res.1t2f$a2comb==a2mat), stg1=mean(res.1t2f$a1comb==a1mat))
bias.1t2fs2 <- log(res.1t2f$yhat2comb)-t2mat
bias.1t2fs1 <- log(res.1t2f$yhat1comb)-tmat
mse.1t2f <- data.frame(stg2=colMeans(bias.1t2fs2^2), stg1=colMeans(bias.1t2fs1^2))
pot.1f2f <- data.frame(stg2=mean(res.1f2f$a2comb==a2mat), stg1=mean(res.1f2f$a1comb==a1mat))
bias.1f2fs2 <- log(res.1f2f$yhat2comb)-t2mat
bias.1f2fs1 <- log(res.1f2f$yhat1comb)-tmat
mse.1f2f <- data.frame(stg2=colMeans(bias.1f2fs2^2), stg1=colMeans(bias.1f2fs1^2))
## pot.mis2 <- data.frame(stg2=mean(res.mis2$a2comb==a2mat), stg1=mean(res.mis2$a1comb==a1mat))
## bias.mis2s2 <- log(res.mis2$yhat2comb)-t2mat
## bias.mis2s1 <- log(res.mis2$yhat1comb)-tmat
## mse.mis2 <- data.frame(stg2=colMeans(bias.mis2s2^2), stg1=colMeans(bias.mis2s1^2))
## pot.mis3 <- data.frame(stg2=mean(res.mis3$a2comb==a2mat), stg1=mean(res.mis3$a1comb==a1mat))
## bias.mis3s2 <- log(res.mis3$yhat2comb)-t2mat
## bias.mis3s1 <- log(res.mis3$yhat1comb)-tmat
## mse.mis3 <- data.frame(stg2=colMeans(bias.mis3s2^2), stg1=colMeans(bias.mis3s1^2))
## pot.mis4 <- data.frame(stg2=mean(res.mis4$a2comb==a2mat), stg1=mean(res.mis4$a1comb==a1mat))
## bias.mis4s2 <- log(res.mis4$yhat2comb)-t2mat
## bias.mis4s1 <- log(res.mis4$yhat1comb)-tmat
## mse.mis4 <- data.frame(stg2=colMeans(bias.mis4s2^2), stg1=colMeans(bias.mis4s1^2))
## pot.mis5 <- data.frame(stg2=mean(res.mis5$a2comb==a2mat), stg1=mean(res.mis5$a1comb==a1mat))
## bias.mis5s2 <- log(res.mis5$yhat2comb)-t2mat
## bias.mis5s1 <- log(res.mis5$yhat1comb)-tmat
## mse.mis5 <- data.frame(stg2=colMeans(bias.mis5s2^2), stg1=colMeans(bias.mis5s1^2))
pot.b <- data.frame(stg2=mean(res.b$a2comb==a2mat), stg1=mean(res.b$a1comb==a1mat))
bias.b2 <- res.b$yhat2comb-t2mat
bias.b1 <- res.b$yhat1comb-tmat
mse.b <- data.frame(stg2=colMeans(bias.b2^2), stg1=colMeans(bias.b1^2))

## subj1y2 <- res.q$yhat2comb[,1]
## subj1true2 <- exp(t2opt[1])
## var(subj1y2)*(rep-1)/rep   #Var
## mean((subj1y2-subj1true2)^2)  #MSE
## (mean(subj1y2)-subj1true2)^2  #Bias^2
## (mean(subj1y2)-subj1true2)^2 + var(subj1y2)*(rep-1)/rep


method <- rep(c("Oracle","1F2Tru","1Tru2F","1F2F","BART"), each=2)  #,"mis2","mis3","mis4","mis5"
split <- factor(0:1, label=c("Bias^2","Variance"))
value <- rbind(c(mean(colMeans(bias.1t2ts2)^2),mean(colMeans(bias.1t2ts1)^2)), colMeans(mse.1t2t)-c(mean(colMeans(bias.1t2ts2)^2),mean(colMeans(bias.1t2ts1)^2)), c(mean(colMeans(bias.1f2ts2)^2),mean(colMeans(bias.1f2ts1)^2)), colMeans(mse.1f2t)-c(mean(colMeans(bias.1f2ts2)^2),mean(colMeans(bias.1f2ts1)^2)), c(mean(colMeans(bias.1t2fs2)^2),mean(colMeans(bias.1t2fs1)^2)), colMeans(mse.1t2f)-c(mean(colMeans(bias.1t2fs2)^2),mean(colMeans(bias.1t2fs1)^2)), c(mean(colMeans(bias.1f2fs2)^2),mean(colMeans(bias.1f2fs1)^2)), colMeans(mse.1f2f)-c(mean(colMeans(bias.1f2fs2)^2),mean(colMeans(bias.1f2fs1)^2)),  c(mean(colMeans(bias.b2)^2),mean(colMeans(bias.b1)^2)), colMeans(mse.b)-c(mean(colMeans(bias.b2)^2),mean(colMeans(bias.b1)^2)))  #, c(mean(colMeans(bias.mis2s2)^2),mean(colMeans(bias.mis2s1)^2)), colMeans(mse.mis2)-c(mean(colMeans(bias.mis2s2)^2),mean(colMeans(bias.mis2s1)^2)), c(mean(colMeans(bias.mis3s2)^2),mean(colMeans(bias.mis3s1)^2)), colMeans(mse.mis3)-c(mean(colMeans(bias.mis3s2)^2),mean(colMeans(bias.mis3s1)^2)), c(mean(colMeans(bias.mis4s2)^2),mean(colMeans(bias.mis4s1)^2)), colMeans(mse.mis4)-c(mean(colMeans(bias.mis4s2)^2),mean(colMeans(bias.mis4s1)^2)), c(mean(colMeans(bias.mis5s2)^2),mean(colMeans(bias.mis5s1)^2)), colMeans(mse.mis5)-c(mean(colMeans(bias.mis5s2)^2),mean(colMeans(bias.mis5s1)^2))
pot <- rbind(pot.1t2t, pot.1f2t, pot.1t2f, pot.1f2f, pot.b)  #, pot.mis2, pot.mis3, pot.mis4, pot.mis5
colnames(pot) <- paste0("pot.", colnames(pot))
pot <- pot[rep(1:nrow(pot), each=2),]
pot.1t2tall <- mean(res.1t2t$a2comb==a2mat & res.1t2t$a1comb==a1mat)
pot.1f2tall <- mean(res.1f2t$a2comb==a2mat & res.1f2t$a1comb==a1mat)
pot.1t2fall <- mean(res.1t2f$a2comb==a2mat & res.1t2f$a1comb==a1mat)
pot.1f2fall <- mean(res.1f2f$a2comb==a2mat & res.1f2f$a1comb==a1mat)
pot.ball <- mean(res.b$a2comb==a2mat & res.b$a1comb==a1mat)
pot.all <- rep(c(pot.1t2tall, pot.1f2tall, pot.1t2fall, pot.1f2fall, pot.ball), each=2)
bardata <- data.frame(method, split, value, pot, pot.all)

stg2plot <- ggplot(bardata)  +
    theme(axis.text=element_text(size=16,face="bold"), axis.title=element_text(size=18,face="bold"), title=element_text(size=18,face="bold"), legend.position="bottom", legend.text=element_text(size=12,face="bold")) +
    geom_bar(aes(fill=split, y=stg2, x=method), position="stack", stat="identity") + xlab("") + ylab("MSE") + labs(fill="") +
    geom_line(aes(y=pot.stg2/16, x=method, group=1), size=1, color="blue") +
    geom_point(aes(y=pot.stg2/16, x=method, color="Stage 2"), size=3, shape=16) +
    geom_line(aes(y=pot.all/16, x=method, group=1), size=1, color="darkgreen") +
    geom_point(aes(y=pot.all/16, x=method, color="Overall"), size=3, shape=16) +
    scale_y_continuous(sec.axis = sec_axis(~.*16, name="POT")) +
    scale_color_manual(values=c("Stage 2"="blue", "Overall"="darkgreen")) + labs(color="POT") +
    ggtitle("Stage 2") 
stg1plot <- ggplot(bardata)  +
    theme(axis.text=element_text(size=16,face="bold"), axis.title=element_text(size=18,face="bold"), title=element_text(size=18,face="bold"), legend.position="bottom", legend.text=element_text(size=12,face="bold")) +
    geom_bar(aes(fill=split, y=stg1, x=method), position="stack", stat="identity") + xlab("") + ylab("MSE") + labs(fill="") +
    geom_line(aes(y=pot.stg1/25, x=method, group=1), size=1, color="red") +
    geom_point(aes(y=pot.stg1/25, x=method, color="Stage 1"), size=3, shape=16) + labs(fill="") +
    geom_line(aes(y=pot.all/25, x=method, group=1), size=1, color="darkgreen") +
    geom_point(aes(y=pot.all/25, x=method, color="Overall"), size=3, shape=16) +
    scale_y_continuous(sec.axis = sec_axis(~.*25, name="POT")) +
    scale_color_manual(values=c("Stage 1"="red", "Overall"="darkgreen")) + labs(color="POT") +
    ggtitle("Stage 1")

pdf(file="~/BART/a3/outputs/lin.pdf", width=18, height=9)
ggarrange(stg2plot, stg1plot, ncol=2, nrow=1)
dev.off()

## meth <- rep(factor(0:1, label=c("Q-learning","BART wrapper")), each=nrow(pot.b))
## pot.stg2 <- data.frame(meth=meth, y=c(pot.q[,1],pot.b[,1]))
## pot.stg1 <- data.frame(meth=meth, y=c(pot.q[,2],pot.b[,2]))
## bias.stg2 <- data.frame(meth=meth, y=c(bias.q[,1],bias.b[,1]))
## bias.stg1 <- data.frame(meth=meth, y=c(bias.q[,2],bias.b[,2]))
## rmse.stg2 <- data.frame(meth=meth, y=c(rmse.q[,1],rmse.b[,1]))
## rmse.stg1 <- data.frame(meth=meth, y=c(rmse.q[,2],rmse.b[,2]))

## par(mfrow=c(2,3))
## boxplot(y~meth, data=pot.stg2, pch=20, col=c("white","gray"), frame.plot=TRUE, ylim=c(0,1), xlab="", ylab="POT", main="Stg2")
## boxplot(y~meth, data=bias.stg2, pch=20, col=c("white","gray"), frame.plot=TRUE, xlab="", ylab="Bias", main="Stg2")
## abline(h=0, col=2)
## boxplot(y~meth, data=rmse.stg2, pch=20, col=c("white","gray"), frame.plot=TRUE, xlab="", ylab="RMSE", main="Stg2")
## boxplot(y~meth, data=pot.stg1, pch=20, col=c("white","gray"), frame.plot=TRUE, ylim=c(0,1), xlab="", ylab="POT", main="Stg1")
## boxplot(y~meth, data=bias.stg1, pch=20, col=c("white","gray"), frame.plot=TRUE, xlab="", ylab="Bias", main="Stg1")
## abline(h=0, col=2)
## boxplot(y~meth, data=rmse.stg1, pch=20, col=c("white","gray"), frame.plot=TRUE, xlab="", ylab="RMSE", main="Stg1")
