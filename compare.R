
library(ggplot2)
library(ggpubr)
qlearn <- readRDS("~/BART/a3/results/qlearn.RDS")
w2 <- readRDS("~/BART/a3/results/w2.RDS")
w3 <- readRDS("~/BART/a3/results/w3.RDS")
w4 <- readRDS("~/BART/a3/results/w4.RDS")
w5 <- readRDS("~/BART/a3/results/w5.RDS")
newdata <- read.csv(file="~/BART/a3/test/data_lin.csv")
a2true <- as.numeric(-0.7+0.5*newdata$x2-0.9*newdata$b2>0)
t2opt <- log(newdata$t2opt_true)
a1true <- as.numeric(0.1-0.2*newdata$x1+0.6*newdata$b1>0)
t <- log(newdata$t_true)
rep <- length(qlearn)
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

res.q <- sumr(rep, qlearn)
res.w2 <- sumr(rep, w2)
res.w3 <- sumr(rep, w3)
res.w4 <- sumr(rep, w4)

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


pot.q <- data.frame(stg2=mean(res.q$a2comb==a2mat), stg1=mean(res.q$a1comb==a1mat))
bias.q2 <- log(res.q$yhat2comb)-t2mat
bias.q1 <- log(res.q$yhat1comb)-tmat
mse.q <- data.frame(stg2=colMeans(bias.q2^2), stg1=colMeans(bias.q1^2))
pot.w2 <- data.frame(stg2=mean(res.w2$a2comb==a2mat), stg1=mean(res.w2$a1comb==a1mat))
bias.w22 <- log(res.w2$yhat2comb)-t2mat
bias.w21 <- log(res.w2$yhat1comb)-tmat
mse.w2 <- data.frame(stg2=colMeans(bias.w22^2), stg1=colMeans(bias.w21^2))
pot.w3 <- data.frame(stg2=mean(res.w3$a2comb==a2mat), stg1=mean(res.w3$a1comb==a1mat))
bias.w32 <- log(res.w3$yhat2comb)-t2mat
bias.w31 <- log(res.w3$yhat1comb)-tmat
mse.w3 <- data.frame(stg2=colMeans(bias.w32^2), stg1=colMeans(bias.w31^2))
pot.w4 <- data.frame(stg2=mean(res.w4$a2comb==a2mat), stg1=mean(res.w4$a1comb==a1mat))
bias.w42 <- log(res.w4$yhat2comb)-t2mat
bias.w41 <- log(res.w4$yhat1comb)-tmat
mse.w4 <- data.frame(stg2=colMeans(bias.w42^2), stg1=colMeans(bias.w41^2))
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


method <- rep(c("Q-learning","w2","w3","w4"), each=2) #"BART wrapper",
split <- factor(0:1, label=c("Bias^2","Variance"))
value <- rbind(c(mean(colMeans(bias.q2)^2),mean(colMeans(bias.q1)^2)), colMeans(mse.q)-c(mean(colMeans(bias.q2)^2),mean(colMeans(bias.q1)^2)), c(mean(colMeans(bias.w22)^2),mean(colMeans(bias.w21)^2)), colMeans(mse.w2)-c(mean(colMeans(bias.w22)^2),mean(colMeans(bias.w21)^2)), c(mean(colMeans(bias.w32)^2),mean(colMeans(bias.w31)^2)), colMeans(mse.w3)-c(mean(colMeans(bias.w32)^2),mean(colMeans(bias.w31)^2)), c(mean(colMeans(bias.w42)^2),mean(colMeans(bias.w41)^2)), colMeans(mse.w4)-c(mean(colMeans(bias.w42)^2),mean(colMeans(bias.w41)^2)))
#c(mean(colMeans(bias.b2)^2),mean(colMeans(bias.b1)^2)), colMeans(mse.b)-c(mean(colMeans(bias.b2)^2),mean(colMeans(bias.b1)^2)), 
pot <- rbind(pot.q, pot.w2, pot.w3, pot.w4) #pot.b, 
colnames(pot) <- paste0("pot.", colnames(pot))
pot <- pot[rep(1:nrow(pot), each=2),]
bardata <- data.frame(method, split, value, pot)

stg2plot <- ggplot(bardata)  +
    theme(axis.text=element_text(size=16,face="bold"), axis.title=element_text(size=18,face="bold"), title=element_text(size=18,face="bold"), legend.position="bottom", legend.text=element_text(size=12,face="bold")) +
    geom_bar(aes(fill=split, y=stg2, x=method), position="stack", stat="identity") + xlab("") + ylab("MSE") +
    theme(legend.position="bottom") + labs(fill="") +
    geom_line(aes(y=pot.stg2/20, x=method, group=1), size=1, color="darkblue") +
    geom_point(aes(y=pot.stg2/20, x=method), size=3, color="red", shape=16) +
    scale_y_continuous(sec.axis = sec_axis(~.*20, name="POT")) +
    ggtitle("Stage 2") 
stg1plot <- ggplot(bardata)  +
    theme(axis.text=element_text(size=16,face="bold"), axis.title=element_text(size=18,face="bold"), title=element_text(size=18,face="bold"), legend.position="bottom", legend.text=element_text(size=12,face="bold")) +
    geom_bar(aes(fill=split, y=stg1, x=method), position="stack", stat="identity") + xlab("") + ylab("MSE") +
    labs(fill="") +
    geom_line(aes(y=pot.stg1/10, x=method, group=1), size=1, color="darkblue") +
    geom_point(aes(y=pot.stg1/10, x=method), size=3, color="red", shape=16) +
    scale_y_continuous(sec.axis = sec_axis(~.*10, name="POT")) +
    ggtitle("Stage 1")
#pdf(file="~/BART/a3/outputs/lin.pdf", width=18, height=9)
ggarrange(stg2plot, stg1plot, ncol=2, nrow=1)
dev.off()

meth <- rep(factor(0:1, label=c("Q-learning","BART wrapper")), each=nrow(pot.b))
pot.stg2 <- data.frame(meth=meth, y=c(pot.q[,1],pot.b[,1]))
pot.stg1 <- data.frame(meth=meth, y=c(pot.q[,2],pot.b[,2]))
bias.stg2 <- data.frame(meth=meth, y=c(bias.q[,1],bias.b[,1]))
bias.stg1 <- data.frame(meth=meth, y=c(bias.q[,2],bias.b[,2]))
rmse.stg2 <- data.frame(meth=meth, y=c(rmse.q[,1],rmse.b[,1]))
rmse.stg1 <- data.frame(meth=meth, y=c(rmse.q[,2],rmse.b[,2]))

par(mfrow=c(2,3))
boxplot(y~meth, data=pot.stg2, pch=20, col=c("white","gray"), frame.plot=TRUE, ylim=c(0,1), xlab="", ylab="POT", main="Stg2")
boxplot(y~meth, data=bias.stg2, pch=20, col=c("white","gray"), frame.plot=TRUE, xlab="", ylab="Bias", main="Stg2")
abline(h=0, col=2)
boxplot(y~meth, data=rmse.stg2, pch=20, col=c("white","gray"), frame.plot=TRUE, xlab="", ylab="RMSE", main="Stg2")
boxplot(y~meth, data=pot.stg1, pch=20, col=c("white","gray"), frame.plot=TRUE, ylim=c(0,1), xlab="", ylab="POT", main="Stg1")
boxplot(y~meth, data=bias.stg1, pch=20, col=c("white","gray"), frame.plot=TRUE, xlab="", ylab="Bias", main="Stg1")
abline(h=0, col=2)
boxplot(y~meth, data=rmse.stg1, pch=20, col=c("white","gray"), frame.plot=TRUE, xlab="", ylab="RMSE", main="Stg1")
