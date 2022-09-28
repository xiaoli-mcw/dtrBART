
library(ggplot2)
library(ggpubr)
qoracle <- readRDS(file="~/BART/a3/results/qoracle.RDS")
qmislin <- readRDS(file="~/BART/a3/results/qmislin.RDS")
qmisint <- readRDS(file="~/BART/a3/results/qmisint.RDS")
newdata <- read.csv(file="~/BART/a3/test/data_nln.csv")
t2opt <- log(newdata$t2opt_true)
t <- log(newdata$t_true)
rep <- length(qoracle)
a2mat <- matrix(replicate(rep,newdata$a2true), nrow=rep, byrow=TRUE)
a1mat <- matrix(replicate(rep,newdata$a1true), nrow=rep, byrow=TRUE)
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

res.ora <- sumr(rep, qoracle)
res.lin <- sumr(rep, qmislin)
res.int <- sumr(rep, qmisint)

a2.b <- y2.b <- cr2.b <- a1.b <- y1.b <- cr1.b <- NULL
for (i in 1:rep){
    bart1 <- readRDS(file=paste0("~/BART/a3/results/bart2/bart2_",i,".RDS"))
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
jpeg(file="~/BART/a3/outputs/nln_cr2.jpg", width=480, height=480)
boxplot(cr~grp, data=crdata, boxwex=0.4, pch=20, col=c("white","gray"), frame.plot=TRUE, xlab="", ylab="95% Coverage Rate", ylim=c(0,1))
abline(h=0.95, lty=2, col=2, lwd=2)
dev.off()

pot.ora <- data.frame(stg2=mean(res.ora$a2comb==a2mat), stg1=mean(res.ora$a1comb==a1mat))
bias.oras2 <- log(res.ora$yhat2comb)-t2mat
bias.oras1 <- log(res.ora$yhat1comb)-tmat
mse.ora <- data.frame(stg2=colMeans(bias.oras2^2), stg1=colMeans(bias.oras1^2))
pot.lin <- data.frame(stg2=mean(res.lin$a2comb==a2mat), stg1=mean(res.lin$a1comb==a1mat))
bias.lins2 <- log(res.lin$yhat2comb)-t2mat
bias.lins1 <- log(res.lin$yhat1comb)-tmat
mse.lin <- data.frame(stg2=colMeans(bias.lins2^2), stg1=colMeans(bias.lins1^2))
pot.int <- data.frame(stg2=mean(res.int$a2comb==a2mat), stg1=mean(res.int$a1comb==a1mat))
bias.ints2 <- log(res.int$yhat2comb)-t2mat
bias.ints1 <- log(res.int$yhat1comb)-tmat
mse.int <- data.frame(stg2=colMeans(bias.ints2^2), stg1=colMeans(bias.ints1^2))
pot.b <- data.frame(stg2=mean(res.b$a2comb==a2mat), stg1=mean(res.b$a1comb==a1mat))
bias.b2 <- res.b$yhat2comb-t2mat
bias.b1 <- res.b$yhat1comb-tmat
mse.b <- data.frame(stg2=colMeans(bias.b2^2), stg1=colMeans(bias.b1^2))

method <- rep(factor(1:4, label=c("Linear","Interaction","BART","Oracle")), each=2) #
split <- factor(0:1, label=c("Bias^2","Variance"))
value <- rbind(c(mean(colMeans(bias.lins2)^2),mean(colMeans(bias.lins1)^2)), colMeans(mse.lin)-c(mean(colMeans(bias.lins2)^2),mean(colMeans(bias.lins1)^2)), c(mean(colMeans(bias.ints2)^2),mean(colMeans(bias.ints1)^2)), colMeans(mse.int)-c(mean(colMeans(bias.ints2)^2),mean(colMeans(bias.ints1)^2)), c(mean(colMeans(bias.b2)^2),mean(colMeans(bias.b1)^2)), colMeans(mse.b)-c(mean(colMeans(bias.b2)^2),mean(colMeans(bias.b1)^2)), c(mean(colMeans(bias.oras2)^2),mean(colMeans(bias.oras1)^2)), colMeans(mse.ora)-c(mean(colMeans(bias.oras2)^2),mean(colMeans(bias.oras1)^2))) #
pot <- rbind(pot.lin, pot.int, pot.b, pot.ora) #
colnames(pot) <- paste0("pot.", colnames(pot))
pot <- pot[rep(1:nrow(pot), each=2),]
pot.oraall <- mean(res.ora$a2comb==a2mat & res.ora$a1comb==a1mat)
pot.linall <- mean(res.lin$a2comb==a2mat & res.lin$a1comb==a1mat)
pot.intall <- mean(res.int$a2comb==a2mat & res.int$a1comb==a1mat)
pot.ball <- mean(res.b$a2comb==a2mat & res.b$a1comb==a1mat)
pot.all <- rep(c(pot.linall, pot.intall, pot.ball, pot.oraall), each=2)
bardata <- data.frame(method, split, value, pot, pot.all)

stg2plot <- ggplot(bardata)  +
    theme(axis.text=element_text(size=16,face="bold"), axis.title=element_text(size=18,face="bold"), title=element_text(size=18,face="bold"), legend.position="bottom", legend.text=element_text(size=12,face="bold")) +
    geom_bar(aes(fill=split, y=stg2, x=method), position="stack", stat="identity") + scale_fill_grey(start=0.4, end=0.7) + xlab("") + xlab("") + ylab("MSE") + labs(fill="") +
    geom_line(aes(y=pot.stg2/1.2, x=method, group=1), size=1) +
    geom_point(aes(y=pot.stg2/1.2, x=method, shape="Stage 2"), size=3) +
    geom_line(aes(y=pot.all/1.2, x=method, group=1), size=1, linetype="dashed") +
    geom_point(aes(y=pot.all/1.2, x=method, shape="Overall"), size=3) +
    scale_y_continuous(sec.axis = sec_axis(~.*1.2, name="POT")) +
    scale_shape_manual(values=c("Stage 2"=2, "Overall"=1)) + labs(shape="POT") +
    ggtitle("Stage 2") 
stg1plot <- ggplot(bardata)  +
    theme(axis.text=element_text(size=16,face="bold"), axis.title=element_text(size=18,face="bold"), title=element_text(size=18,face="bold"), legend.position="bottom", legend.text=element_text(size=12,face="bold")) +
    geom_bar(aes(fill=split, y=stg1, x=method), position="stack", stat="identity") + scale_fill_grey(start=0.4, end=0.7) + xlab("") + xlab("") + ylab("MSE") + labs(fill="") +
    geom_line(aes(y=pot.stg1/5.5, x=method, group=1), size=1) +
    geom_point(aes(y=pot.stg1/5.5, x=method, shape="Stage 1"), size=3) +
    geom_line(aes(y=pot.all/5.5, x=method, group=1), size=1, linetype="dashed") +
    geom_point(aes(y=pot.all/5.5, x=method, shape="Overall"), size=3) +
    scale_y_continuous(sec.axis = sec_axis(~.*5.5, name="POT")) +
    scale_shape_manual(values=c("Stage 1"=2, "Overall"=1)) + labs(shape="POT") +
    ggtitle("Stage 1")

pdf(file="~/BART/a3/outputs/nln2.pdf", width=18, height=9)
ggarrange(stg2plot, stg1plot, ncol=2, nrow=1)
dev.off()
