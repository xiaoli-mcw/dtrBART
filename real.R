library(BART3)
library(survival)
library(parallel)
library(doRNG)
library(doParallel)
source("~/BART/a3/codes/wrapper.R")
realdata <- read.csv(file="~/BART/a3/data/all.comp.cases.csv")
## > str(realdata)
## 'data.frame':	4171 obs. of  21 variables:
##  $ dummy_id             : num  1.69e+07 2.75e+07 8.72e+07 1.38e+08 1.38e+08 ...
##  $ Int2    treated in the 2nd stage             : int  0 1 0 0 0 1 1 0 0 1 ...
##  $ Y2     dfs              : num  1.51 1.15 12.96 52.6 0.69 ...
##  $ delta2               : int  1 1 1 0 1 1 0 1 1 1 ...
##  $ A2.NHTL  Action in the 2nd stage            : int  0 0 0 0 0 0 0 0 0 0 ...
##  $ newage    Age(0-9;10-39;40+)          : int  0 1 0 1 2 0 2 2 2 0 ...
##  $ kpspoor   Karnofsky/Lansky score           : int  0 1 0 0 0 1 0 0 0 0 ...
##  $ disease.risk.combined   Disease status    : chr  "Early" "Advanced" "Advanced" "Intermediate" ...
##  $ related     Donor relation         : chr  "Unrelated" "Unrelated" "Unrelated" "Unrelated" ...
##  $ female.to.male  Donor-recipient sex     : int  0 0 0 0 0 0 0 0 0 0 ...
##  $ graftype    Graft type         : chr  "Bone marrow" "Cord blood" "Bone marrow" "Bone marrow" ...
##  $ rec_match    HLA-match        : int  0 2 2 1 0 0 0 1 2 1 ...
##  $ tbi       TBI           : chr  "Yes" "Yes" "Yes" "Yes" ...
##  $ cmv.pair    CMV status         : chr  "Negative-Negative" "Donor or Recipient Positive" "Donor or Recipient Positive" "Donor or Recipient Positive" ...
##  $ condint.binary  Intensity conditioning     : chr  "MA" "MA" "MA" "MA" ...
##  $ pgvhcor              : chr  "Yes" "No" "Yes" "No" ...
##  $ four.or.more.ISP_new  Stage 2 only : int  0 1 0 0 0 0 1 0 0 0 ...
##  $ time.to.acute.GVHD   <1mon or >=1mon  Stage 2 only : int  1 1 1 1 0 1 1 1 1 0 ...
##  $ A1.NHTL  Action in the 1st stage            : int  0 0 0 0 0 1 0 0 0 1 ...
##  $ Y1.tilde   time to A2          : num  2.237 0.691 13.322 53.158 5.822 ...
##  $ delta1               : int  1 1 1 0 1 1 1 1 1 1 ...

varnames <- names(realdata)
varlist <- varnames[c(8,9,11,13,14,15,16)]
realdata[,varlist] <- lapply(realdata[,varlist], as.factor)

x1 <- varnames[6:16]  # covariates in Stage 1
x2 <- varnames[6:18]  # covariates in Stage 2
##########
## BART
##########
res <- dtr1(x1=x1, a1="A1.NHTL", time1="Y1.tilde", x2=x2, a2="A2.NHTL", time2="Y2", stg2ind="Int2", delta="delta2", data=realdata, opt=FALSE)
## saveRDS(res, file="~/BART/a3/results/real.RDS")
res <- readRDS(file="~/BART/a3/results/real.RDS")  # yhat for each action at each stage + optimal action/yhat are returned

## waterfall
library(ggplot2)
yhatdiff <- exp(res$A1.NHTL_1)-exp(res$A1.NHTL_0)
lowerq <- apply(yhatdiff, 1, quantile, probs=0.025)
upperq <- apply(yhatdiff, 1, quantile, probs=0.975)
q1 <- apply(yhatdiff, 1, quantile, probs=0.25)
q3 <- apply(yhatdiff, 1, quantile, probs=0.75)

output <- data.frame(stg1_diff=rowMeans(yhatdiff), lowerq, upperq, q1, q3)
sortout <- output[order(output$stg1_diff, decreasing=TRUE),]
sortout$x <- 1:nrow(sortout)
## sortold <- sortout
## sortout <- sortout[-c(1:1705,3502:4171),]
## 1715,1727,1731,1750,1753,1769,1770,1779

sortout$width <- sortout$upperq-sortout$lowerq
sortout$ratio <- sortout$width/sortout$stg1_diff
sortout$outlier <- sortout$x %in% 1700:3500 & sortout$ratio>10

sortout <- cbind(realdata[rownames(sortout),], sortout)

b <- ggplot(sortout, aes(x=x, y=stg1_diff)) + theme(axis.text=element_text(size=16,face="bold"),axis.title=element_text(size=18,face="bold"),title=element_text(size=18,face="bold")) + labs(list(x = "Patient Number \n(Ordered by risk difference)", y = "Risk difference")) 
smoother <- geom_line(color="black",size=2)
quartiles <- geom_pointrange(data=sortout, aes(x=x, y=stg1_diff, ymax=q3, ymin=q1, colour=related))
p1a <- b+ geom_pointrange(data=sortout, aes(x=x, y=stg1_diff, ymax=upperq,ymin=lowerq),colour="gray50")+quartiles+smoother+scale_y_continuous("Posterior median difference in predicted DFS time")+scale_x_continuous("Patient number (Ordered)")+ggtitle("Stage 1")
pdf(file="~/BART/a3/outputs/real_stg1related.pdf", width=8,height=8)
p1a
dev.off()

outliers <- sortout[sortout$outlier,]
a1opt.est <- as.numeric(apply(res$a1.opt[as.numeric(rownames(outliers)),], 1, function(x) names(sort(table(x), decreasing=TRUE)[1])))
out <- cbind(realdata[as.numeric(rownames(outliers)),c("A1.NHTL","A2.NHTL","delta1","delta2")], a1opt.est, outliers)
head(out)

out <- cbind(realdata[as.numeric(rownames(outliers)),], a1opt.est, outliers)
head(out)

b <- ggplot(sortout, aes(x=x, y=stg1_diff)) + theme(axis.text=element_text(size=16,face="bold"),axis.title=element_text(size=18,face="bold"),title=element_text(size=18,face="bold")) + labs(list(x = "Patient Number \n(Ordered by risk difference)", y = "Risk difference")) 
smoother <- geom_line(color="black",size=2)
quartiles <- geom_pointrange(data=sortout, aes(x=x, y=stg1_diff, ymax=q3, ymin=q1), colour="gray30")
p1a <- b+ geom_pointrange(data=sortout, aes(x=x, y=stg1_diff, ymax=upperq,ymin=lowerq),colour="gray50")+quartiles+smoother+scale_y_continuous("Posterior median difference in predicted DFS time")+scale_x_continuous("Patient number (Ordered)")+ggtitle("Stage 1")
## +coord_cartesian(ylim=c(-0.25,0.25))
## ,breaks=c(-0.2,-0.1,0,0.1,0.2)
pdf(file="~/BART/a3/outputs/real_stg1dfs.pdf", width=8,height=8)
p1a
dev.off()

yhatdiff <- exp(res$A2.NHTL_1)-exp(res$A2.NHTL_0)
lowerq <- apply(yhatdiff, 1, quantile, probs=0.025)
upperq <- apply(yhatdiff, 1, quantile, probs=0.975)
q1 <- apply(yhatdiff, 1, quantile, probs=0.25)
q3 <- apply(yhatdiff, 1, quantile, probs=0.75)

output <- data.frame(stg2_diff=rowMeans(yhatdiff), lowerq, upperq, q1, q3)
sortout <- output[order(output$stg2_diff, decreasing=TRUE),]
sortout$x <- 1:nrow(sortout)

b <- ggplot(sortout, aes(x=x, y=stg2_diff)) + theme(axis.text=element_text(size=16,face="bold"),axis.title=element_text(size=18,face="bold"),title=element_text(size=18,face="bold")) + labs(list(x = "Patient Number \n(Ordered by risk difference)", y = "Risk difference")) 
smoother <- geom_line(color="black",size=2)
quartiles <- geom_pointrange(data=sortout, aes(x=x, y=stg2_diff, ymax=q3, ymin=q1), colour="gray30")
p2a <- b+ geom_pointrange(data=sortout, aes(x=x, y=stg2_diff, ymax=upperq,ymin=lowerq),colour="gray50")+quartiles+smoother+scale_y_continuous("Posterior median difference in predicted DFS time")+scale_x_continuous("Patient number (Ordered)")+ggtitle("Stage 2")
## +coord_cartesian(ylim=c(-0.25,0.25))
## ,breaks=c(-0.2,-0.1,0,0.1,0.2)
pdf(file="~/BART/a3/outputs/real_stg2dfs.pdf", width=8,height=8)
p2a
dev.off()

## randomly 300 posterior samples
## set.seed(71422)
## samples <- sample(1:ncol(res$A2.NHTL_0), size=300)
## diff1 <- exp(res$A1.NHTL_1[,samples])-exp(res$A1.NHTL_0[,samples])
## diff2 <- exp(res$A2.NHTL_1[,samples])-exp(res$A2.NHTL_0[,samples])

## png(file="~/BART/a3/outputs/real_bart.png", width=600, height=900)
## par(mfrow=c(2,2))
## plot(1, type="n", xlim=c(0, nrow(diff1)), ylim=range(stg1_diff), xlab="Patient Index", ylab="Difference in predicted DFS time", main="NHTL vs. Standard at Stage 1 with original sample")
## for (i in 1:length(samples)) lines(diff1[subjorder1,i], col=5)
## lines(sort(stg1_diff), type='l', lwd=2)
## abline(h=0, col=2, lty=2)
## plot(1, type="n", xlim=c(0, nrow(diff1)), ylim=range(stg1_diff), xlab="Patient Index", ylab="Difference in predicted DFS time", main="NHTL vs. Standard at Stage 1 with sorted sample")
## for (i in 1:length(samples)) lines(sort(diff1[,i]), col=5)
## lines(sort(stg1_diff), type='l', lwd=2)
## abline(h=0, col=2, lty=2)

## plot(1, type="n", xlim=c(0, nrow(diff2)), ylim=range(stg2_diff), xlab="Patient Index", ylab="Difference in predicted DFS time", main="NHTL vs. Standard at Stage 2 with original sample")
## for (i in 1:length(samples)) lines(diff2[subjorder2,i], col=5)
## lines(sort(stg2_diff), type='l', lwd=2)
## abline(h=0, col=2, lty=2)
## plot(1, type="n", xlim=c(0, nrow(diff2)), ylim=range(stg2_diff), xlab="Patient Index", ylab="Difference in predicted DFS time", main="NHTL vs. Standard at Stage 2 with sorted sample")
## for (i in 1:length(samples)) lines(sort(diff2[,i]), col=5)
## lines(sort(stg2_diff), type='l', lwd=2)
## abline(h=0, col=2, lty=2)
## dev.off()

## fit the fit
library(rpart)
fitteddiff <- data.frame(y=stg1_diff, realdata)
(fmla <- as.formula(paste("stg1_diff ~ ", paste(x1, collapse= "+"))))
fit <- rpart(fmla,data=fitteddiff,method="anova",cp=0.01)
printcp(fit)
plotcp(fit)
summary(fit)
pdf(file="~/BART/a3/outputs/Rsq.pdf", width=16, height=8)
par(mfrow=c(1,2)) # two plots on one page
rsq.rpart(fit)
dev.off()

library(rpart.plot)

#### Here is nice looking tree figure

pdf(file="~/BART/a3/outputs/fitthefit.pdf")
rpart.plot(fit,type=4,extra=1,
shadow.col="gray", # add shadows just for kicks
main="Difference in DFS time at Stage 1 \n(NHTL-Standard)\n")
dev.off()


################
## Q-learning
################
data.stg2 <- realdata[realdata$Int2==1,]
stg2n <- nrow(data.stg2)
a2option <- unique(realdata$A2.NHTL)
a2n <- length(a2option)
qcov2 <- paste(paste("A2.NHTL*",c(x2,"A1.NHTL")),collapse="+")
fml2 <- as.formula(paste("Surv(Y2,delta2)", qcov2, sep="~"))
qstg2 <- survreg(fml2, data=data.stg2, dist="lognormal")
x.mat <- data.stg2[,c(x2,"A1.NHTL")]
x.test.stg2 <- data.frame(x.mat[rep(seq_len(nrow(x.mat)), each=a2n),], A2.NHTL=rep(a2option,stg2n))
pred.stg2 <- predict(qstg2, newdata=x.test.stg2)
a2names <- paste0("A2.NHTL_",a2option)
for (i in 1:a2n) assign(a2names[i], pred.stg2[1:stg2n*a2n-(a2n-i)])
a2opt <- sapply(1:stg2n*a2n-1, function(column) {which.max(pred.stg2[column:(column+a2n-1)])})  #opt col at stg2
a2.opt <- sapply(a2opt, function(x) a2option[x])  #opt action at stg2
yhat2optmean <- sapply(1:stg2n*a2n-1, function(column) {max(pred.stg2[column:(column+a2n-1)])})

a1option <- unique(realdata$A1.NHTL)  #all possible stage1 actions
a1n <- length(a1option)  #total number of stage1 actions
data.stg1 <- realdata
n <- nrow(data.stg1)
data.stg1$y2opt <- 0
data.stg1[data.stg1$Int2==1,]$y2opt <- yhat2optmean
data.stg1$Y1 <- data.stg1$Y1.tilde+data.stg1$y2opt
qcov1 <- paste(paste("A1.NHTL*",x1),collapse="+")
fml1 <- as.formula(paste("Surv(Y1,delta1)", qcov1, sep="~"))
qstg1 <- survreg(fml1, data=data.stg1, dist="lognormal")
x.mat <- data.stg1[,x1]
x.test.stg1 <- data.frame(x.mat[rep(seq_len(nrow(x.mat)), each=a1n),], A1.NHTL=rep(a1option,n))
pred.stg1 <- predict(qstg1, newdata=x.test.stg1)
a1names <- paste0("A1.NHTL_",a1option)
for (i in 1:a1n) assign(a1names[i], pred.stg1[1:n*a1n-(a1n-i)])
a1opt <- sapply(1:n*a1n-1, function(column) {which.max(pred.stg1[column:(column+a1n-1)])})  #opt col at stg1
a1.opt <- sapply(a1opt, function(x) a1option[x])  #opt action at stg1
yhat1optmean <- sapply(1:n*a1n-1, function(column) {max(pred.stg1[column:(column+a1n-1)])})

