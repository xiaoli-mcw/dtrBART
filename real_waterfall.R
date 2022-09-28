library(ggplot2)

diff <- function(a1, a0, type="BART", log=TRUE, orig=NULL){
    ## Function to calculate the difference in DFS time
    if (type == "BART") {
        if (log) {yhatdiff <- a1-a0}
        else {yhatdiff <- exp(a1)-exp(a0)}
        lowerq <- apply(yhatdiff, 1, quantile, probs=0.025)
        upperq <- apply(yhatdiff, 1, quantile, probs=0.975)
        q1 <- apply(yhatdiff, 1, quantile, probs=0.25)
        q3 <- apply(yhatdiff, 1, quantile, probs=0.75)
        output <- data.frame(diff=rowMeans(yhatdiff), lowerq, upperq, q1, q3)
    }
    else if (type == "Qlearn") {
        if (log) {
            bootdiff <- log(get(a1))-log(get(a0))
            yhatdiff <- log(orig[[a1]])-log(orig[[a0]])
        }
        else {
            bootdiff <- get(a1)-get(a0)
            yhatdiff <- orig[[a1]]-orig[[a0]]
        }
        sigma <- apply(bootdiff, 1, sd)
        lowerq <- yhatdiff + qnorm(0.025)*sigma
        upperq <- yhatdiff + qnorm(0.975)*sigma
        q1 <- yhatdiff + qnorm(0.25)*sigma
        q3 <- yhatdiff + qnorm(0.75)*sigma
        output <- data.frame(diff=yhatdiff, lowerq, upperq, q1, q3)
    }
    sortout <- output[order(output$diff, decreasing=TRUE),]
    sortout$x <- 1:nrow(sortout)
    return(sortout)
}

sdiff <- function(t, a1, a0, sd){
    ## Function to calculate the difference in survival probability at time t
    s1 <- s0 <- NULL
    rep <- length(sd)
    for (i in 1:rep){
        s1 <- cbind(s1, pnorm(q=t, mean=a1[,i], sd=sd[i], lower.tail=FALSE))
        s0 <- cbind(s0, pnorm(q=t, mean=a0[,i], sd=sd[i], lower.tail=FALSE))
    }
    diff <- s1-s0
    diffmean <- rowMeans(diff)
    diffci <- apply(diff, 1, quantile, probs=c(0.025,0.975,0.25,0.75))
    diffsum <- data.frame(mean=diffmean, lowerq=diffci[1,], upperq=diffci[2,], q1=diffci[3,], q3=diffci[4,])
    diffsum <- diffsum[order(diffsum$mean, decreasing=TRUE),]
    diffsum$x <- 1:nrow(diffsum)
    return(diffsum)
}

##########
## BART
##########
res <- readRDS(file="~/BART/a3/results/real.RDS")

## log(dfs time) diff at stg1
bart_lg_stg1 <- diff(res$A1.NHTL_1,res$A1.NHTL_0)
## log(dfs time) diff at stg2
bart_lg_stg2 <- diff(res$A2.NHTL_1,res$A2.NHTL_0)

## dfs time diff at stg1
bart_dif_stg1 <- diff(res$A1.NHTL_1,res$A1.NHTL_0,log=FALSE)
## dfs time diff at stg2
bart_dif_stg2 <- diff(res$A2.NHTL_1,res$A2.NHTL_0,log=FALSE)

t <- log(2*12)  #original Ys are months
## survival prob diff at stg1
bart_sdif1 <- sdiff(t=t, a1=res$A1.NHTL_1, a0=res$A1.NHTL_0, sd=res$sigma1)
## survival prob diff at stg2
bart_sdif2 <- sdiff(t=t, a1=res$A2.NHTL_1, a0=res$A2.NHTL_0, sd=res$sigma2)

#############
## Q-learn
#############
orig <- readRDS(file="~/BART/a3/results/real_q.RDS")
qlearn <- readRDS(file="~/BART/a3/results/real_q_bootstrap.RDS")
A1.NHTL_0 <- qlearn$A1.NHTL0
A1.NHTL_1 <- qlearn$A1.NHTL1
A2.NHTL_0 <- qlearn$A2.NHTL0
A2.NHTL_1 <- qlearn$A2.NHTL1
sd1 <- qlearn$sd1
sd2 <- qlearn$sd2

## log(dfs time) diff at stg1
q_lg_stg1 <- diff("A1.NHTL_1","A1.NHTL_0",type="Qlearn",orig=orig)
## log(dfs time) diff at stg2
q_lg_stg2 <- diff("A2.NHTL_1","A2.NHTL_0",type="Qlearn",orig=orig)

## dfs time diff at stg1
q_dif_stg1 <- diff("A1.NHTL_1","A1.NHTL_0",type="Qlearn",log=FALSE,orig=orig)
## dfs time diff at stg2
q_dif_stg2 <- diff("A2.NHTL_1","A2.NHTL_0",type="Qlearn",log=FALSE,orig=orig)

## survival prob diff at stg1
q_sdif1 <- sdiff(t=t, a1=log(A1.NHTL_1), a0=log(A1.NHTL_0), sd=sd1)
## survival prob diff at stg1
q_sdif2 <- sdiff(t=t, a1=log(A2.NHTL_1), a0=log(A2.NHTL_0), sd=sd2)

###################
## waterfall plot
###################
ylim_lg_stg1 <- range(bart_lg_stg1$lowerq, bart_lg_stg1$upperq, q_lg_stg1$lowerq, q_lg_stg1$upperq)
b <- ggplot(bart_lg_stg1, aes(x=x, y=diff)) + theme(axis.text=element_text(size=16,face="bold"),axis.title=element_text(size=18,face="bold"),title=element_text(size=18,face="bold")) + labs(list(x = "Patient Number \n(Ordered by risk difference)", y = "Risk difference")) 
smoother <- geom_line(color="black",size=2)
quartiles <- geom_pointrange(data=bart_lg_stg1, aes(x=x, y=diff, ymax=q3, ymin=q1), colour="gray30")
p1a <- b+ geom_pointrange(data=bart_lg_stg1, aes(x=x, y=diff, ymax=upperq,ymin=lowerq),colour="gray50")+quartiles+smoother+scale_y_continuous("Posterior mean difference in predicted log(DFS time)")+scale_x_continuous("Patient number (Ordered)")+ggtitle("Stage 1")+coord_cartesian(ylim=ylim_lg_stg1)
pdf(file="~/BART/a3/outputs/real_stg1dfs_lg_align.pdf", width=8,height=8)
p1a
dev.off()

b <- ggplot(q_lg_stg1, aes(x=x, y=diff)) + theme(axis.text=element_text(size=16,face="bold"),axis.title=element_text(size=18,face="bold"),title=element_text(size=18,face="bold")) + labs(list(x = "Patient Number \n(Ordered by risk difference)", y = "Risk difference")) 
smoother <- geom_line(color="black",size=2)
quartiles <- geom_pointrange(data=q_lg_stg1, aes(x=x, y=diff, ymax=q3, ymin=q1), colour="gray30")
p1a <- b+ geom_pointrange(data=q_lg_stg1, aes(x=x, y=diff, ymax=upperq,ymin=lowerq),colour="gray50")+quartiles+smoother+scale_y_continuous("Bootstrap mean difference in predicted log(DFS time)")+scale_x_continuous("Patient number (Ordered)")+ggtitle("Stage 1")+coord_cartesian(ylim=ylim_lg_stg1)
pdf(file="~/BART/a3/outputs/real_stg1dfs_qlg_align.pdf", width=8,height=8)
p1a
dev.off()

ylim_lg_stg2 <- range(bart_lg_stg2$lowerq, bart_lg_stg2$upperq, q_lg_stg2$lowerq, q_lg_stg2$upperq)
b <- ggplot(bart_lg_stg2, aes(x=x, y=diff)) + theme(axis.text=element_text(size=16,face="bold"),axis.title=element_text(size=18,face="bold"),title=element_text(size=18,face="bold")) + labs(list(x = "Patient Number \n(Ordered by risk difference)", y = "Risk difference")) 
smoother <- geom_line(color="black",size=2)
quartiles <- geom_pointrange(data=bart_lg_stg2, aes(x=x, y=diff, ymax=q3, ymin=q1), colour="gray30")
p2a <- b+ geom_pointrange(data=bart_lg_stg2, aes(x=x, y=diff, ymax=upperq,ymin=lowerq),colour="gray50")+quartiles+smoother+scale_y_continuous("Posterior mean difference in predicted log(DFS time)")+scale_x_continuous("Patient number (Ordered)")+ggtitle("Stage 2")+coord_cartesian(ylim=ylim_lg_stg2)
pdf(file="~/BART/a3/outputs/real_stg2dfs_lg_align.pdf", width=8,height=8)
p2a
dev.off()

b <- ggplot(q_lg_stg2, aes(x=x, y=diff)) + theme(axis.text=element_text(size=16,face="bold"),axis.title=element_text(size=18,face="bold"),title=element_text(size=18,face="bold")) + labs(list(x = "Patient Number \n(Ordered by risk difference)", y = "Risk difference")) 
smoother <- geom_line(color="black",size=2)
quartiles <- geom_pointrange(data=q_lg_stg2, aes(x=x, y=diff, ymax=q3, ymin=q1), colour="gray30")
p2a <- b+ geom_pointrange(data=q_lg_stg2, aes(x=x, y=diff, ymax=upperq,ymin=lowerq),colour="gray50")+quartiles+smoother+scale_y_continuous("Bootstrap mean difference in predicted log(DFS time)")+scale_x_continuous("Patient number (Ordered)")+ggtitle("Stage 2")+coord_cartesian(ylim=ylim_lg_stg2)
pdf(file="~/BART/a3/outputs/real_stg2dfs_qlg_align.pdf", width=8,height=8)
p2a
dev.off()

ylim_stg1 <- range(bart_dif_stg1$lowerq, bart_dif_stg1$upperq, q_dif_stg1$lowerq, q_dif_stg1$upperq)
b <- ggplot(bart_dif_stg1, aes(x=x, y=diff)) + theme(axis.text=element_text(size=16,face="bold"),axis.title=element_text(size=18,face="bold"),title=element_text(size=18,face="bold")) + labs(list(x = "Patient Number \n(Ordered by risk difference)", y = "Risk difference")) 
smoother <- geom_line(color="black",size=2)
quartiles <- geom_pointrange(data=bart_dif_stg1, aes(x=x, y=diff, ymax=q3, ymin=q1), colour="gray30")
p1a <- b+ geom_pointrange(data=bart_dif_stg1, aes(x=x, y=diff, ymax=upperq,ymin=lowerq),colour="gray50")+quartiles+smoother+scale_y_continuous("Posterior median difference in predicted DFS time")+scale_x_continuous("Patient number (Ordered)")+ggtitle("Stage 1")+coord_cartesian(ylim=ylim_stg1)
## ,breaks=c(-0.2,-0.1,0,0.1,0.2)
pdf(file="~/BART/a3/outputs/real_stg1dfs_align.pdf", width=8,height=8)
p1a
dev.off()

b <- ggplot(q_dif_stg1, aes(x=x, y=diff)) + theme(axis.text=element_text(size=16,face="bold"),axis.title=element_text(size=18,face="bold"),title=element_text(size=18,face="bold")) + labs(list(x = "Patient Number \n(Ordered by risk difference)", y = "Risk difference")) 
smoother <- geom_line(color="black",size=2)
quartiles <- geom_pointrange(data=q_dif_stg1, aes(x=x, y=diff, ymax=q3, ymin=q1), colour="gray30")
p1a <- b+ geom_pointrange(data=q_dif_stg1, aes(x=x, y=diff, ymax=upperq,ymin=lowerq),colour="gray50")+quartiles+smoother+scale_y_continuous("Bootstrap median difference in predicted DFS time")+scale_x_continuous("Patient number (Ordered)")+ggtitle("Stage 1")+coord_cartesian(ylim=ylim_stg1)
## ,breaks=c(-0.2,-0.1,0,0.1,0.2)
pdf(file="~/BART/a3/outputs/real_stg1dfs_q_align.pdf", width=8,height=8)
p1a
dev.off()

ylim_stg2 <- range(bart_dif_stg2$lowerq, bart_dif_stg2$upperq, q_dif_stg2$lowerq, q_dif_stg2$upperq)
b <- ggplot(bart_dif_stg2, aes(x=x, y=diff)) + theme(axis.text=element_text(size=16,face="bold"),axis.title=element_text(size=18,face="bold"),title=element_text(size=18,face="bold")) + labs(list(x = "Patient Number \n(Ordered by risk difference)", y = "Risk difference")) 
smoother <- geom_line(color="black",size=2)
quartiles <- geom_pointrange(data=bart_dif_stg2, aes(x=x, y=diff, ymax=q3, ymin=q1), colour="gray30")
p2a <- b+ geom_pointrange(data=bart_dif_stg2, aes(x=x, y=diff, ymax=upperq,ymin=lowerq),colour="gray50")+quartiles+smoother+scale_y_continuous("Posterior median difference in predicted DFS time")+scale_x_continuous("Patient number (Ordered)")+ggtitle("Stage 2")+coord_cartesian(ylim=ylim_stg2)
## ,breaks=c(-0.2,-0.1,0,0.1,0.2)
pdf(file="~/BART/a3/outputs/real_stg2dfs_align.pdf", width=8,height=8)
p2a
dev.off()

b <- ggplot(q_dif_stg2, aes(x=x, y=diff)) + theme(axis.text=element_text(size=16,face="bold"),axis.title=element_text(size=18,face="bold"),title=element_text(size=18,face="bold")) + labs(list(x = "Patient Number \n(Ordered by risk difference)", y = "Risk difference")) 
smoother <- geom_line(color="black",size=2)
quartiles <- geom_pointrange(data=q_dif_stg2, aes(x=x, y=diff, ymax=q3, ymin=q1), colour="gray30")
p2a <- b+ geom_pointrange(data=q_dif_stg2, aes(x=x, y=diff, ymax=upperq,ymin=lowerq),colour="gray50")+quartiles+smoother+scale_y_continuous("Bootstrap median difference in predicted DFS time")+scale_x_continuous("Patient number (Ordered)")+ggtitle("Stage 2")+coord_cartesian(ylim=ylim_stg2)
## ,breaks=c(-0.2,-0.1,0,0.1,0.2)
pdf(file="~/BART/a3/outputs/real_stg2dfs_q_align.pdf", width=8,height=8)
p2a
dev.off()

ylim_s_stg1 <- range(bart_sdif1$lowerq, bart_sdif1$upperq, q_sdif1$lowerq, q_sdif1$upperq)
b <- ggplot(bart_sdif1, aes(x=x, y=mean)) + theme(axis.text=element_text(size=16,face="bold"),axis.title=element_text(size=18,face="bold"),title=element_text(size=18,face="bold"))
smoother <- geom_line(color="black",size=2)
quartiles <- geom_pointrange(data=bart_sdif1, aes(x=x, y=mean, ymax=q3, ymin=q1), colour="gray30")
p2a <- b+ geom_pointrange(data=bart_sdif1, aes(x=x, y=mean, ymax=upperq,ymin=lowerq),colour="gray50")+quartiles+smoother+scale_y_continuous("Difference in 2-yr survival probability")+scale_x_continuous("Patient number (Ordered)")+ggtitle("Stage 1")+coord_cartesian(ylim=ylim_s_stg1)
pdf(file="~/BART/a3/outputs/real_stg1_2yr_align.pdf", width=8,height=8)
p2a
dev.off()

b <- ggplot(q_sdif1, aes(x=x, y=mean)) + theme(axis.text=element_text(size=16,face="bold"),axis.title=element_text(size=18,face="bold"),title=element_text(size=18,face="bold"))
smoother <- geom_line(color="black",size=2)
quartiles <- geom_pointrange(data=q_sdif1, aes(x=x, y=mean, ymax=q3, ymin=q1), colour="gray30")
p2a <- b+ geom_pointrange(data=q_sdif1, aes(x=x, y=mean, ymax=upperq,ymin=lowerq),colour="gray50")+quartiles+smoother+scale_y_continuous("Difference in 2-yr survival probability")+scale_x_continuous("Patient number (Ordered)")+ggtitle("Stage 1")+coord_cartesian(ylim=ylim_s_stg1)
pdf(file="~/BART/a3/outputs/real_stg1_2yr_q_align.pdf", width=8,height=8)
p2a
dev.off()

ylim_s_stg2 <- range(bart_sdif2$lowerq, bart_sdif2$upperq, q_sdif2$lowerq, q_sdif2$upperq)
b <- ggplot(bart_sdif2, aes(x=x, y=mean)) + theme(axis.text=element_text(size=16,face="bold"),axis.title=element_text(size=18,face="bold"),title=element_text(size=18,face="bold"))
smoother <- geom_line(color="black",size=2)
quartiles <- geom_pointrange(data=bart_sdif2, aes(x=x, y=mean, ymax=q3, ymin=q1), colour="gray30")
p2a <- b+ geom_pointrange(data=bart_sdif2, aes(x=x, y=mean, ymax=upperq,ymin=lowerq),colour="gray50")+quartiles+smoother+scale_y_continuous("Difference in 2-yr survival probability")+scale_x_continuous("Patient number (Ordered)")+ggtitle("Stage 2")+coord_cartesian(ylim=ylim_s_stg2)
pdf(file="~/BART/a3/outputs/real_stg2_2yr_align.pdf", width=8,height=8)
p2a
dev.off()

b <- ggplot(q_sdif2, aes(x=x, y=mean)) + theme(axis.text=element_text(size=16,face="bold"),axis.title=element_text(size=18,face="bold"),title=element_text(size=18,face="bold"))
smoother <- geom_line(color="black",size=2)
quartiles <- geom_pointrange(data=q_sdif2, aes(x=x, y=mean, ymax=q3, ymin=q1), colour="gray30")
p2a <- b+ geom_pointrange(data=q_sdif2, aes(x=x, y=mean, ymax=upperq,ymin=lowerq),colour="gray50")+quartiles+smoother+scale_y_continuous("Difference in 2-yr survival probability")+scale_x_continuous("Patient number (Ordered)")+ggtitle("Stage 2")+coord_cartesian(ylim=ylim_s_stg2)
pdf(file="~/BART/a3/outputs/real_stg2_2yr_q_align.pdf", width=8,height=8)
p2a
dev.off()
