library(timeROC)
library(survival)
realdata <- read.csv(file="~/BART/a3/data/all.comp.cases.csv")
res_b <- readRDS(file="~/BART/a3/results/real.RDS")  # yhat for each action at each stage + optimal action/yhat are returned
res_q <- readRDS(file="~/BART/a3/results/real_q.RDS")  # yhat for each action at each stage + optimal action/yhat are returned

stg2 <- realdata$Int2==1
T_stg2 <- realdata[stg2,]$Y2
delta2 <- realdata[stg2,]$delta2
timepoint <- quantile(T_stg2, probs=c(0.5,0.75))
## pick the predicted times that corresponding to the treatment observed
stg2_0 <- realdata[stg2,]$A2.NHTL==0
pred_b_stg2 <- res_b$A2.NHTL_1
pred_b_stg2[stg2_0,] <- res_b$A2.NHTL_0[stg2_0,]
pred_q_stg2 <- res_q$A2.NHTL_1
pred_q_stg2[stg2_0] <- res_q$A2.NHTL_0[stg2_0]
## calculate the time dependent AUC
ROC_b_stg2 <- timeROC(T=T_stg2, delta=delta2, marker=-rowMeans(exp(pred_b_stg2)), cause=1, time=timepoint)
ROC_q_stg2 <- timeROC(T=T_stg2, delta=delta2, marker=-pred_q_stg2, cause=1, time=timepoint)
comp_stg2 <- cbind(bart=ROC_b_stg2$AUC, q=ROC_q_stg2$AUC)
(comp_stg2 <- cbind(comp_stg2, diff=comp_stg2[,1]-comp_stg2[,2]))

a2opt_b <- as.numeric(apply(res_b$a2.opt, 1, function(x) names(sort(table(x), decreasing=TRUE)[1])))
a2opt_q <- res_q$a2.opt
agree <- a2opt_b==a2opt_q
a2opt_agree <- a2opt_b[agree]
match_b <- a2opt_b==realdata[stg2,]$A2.NHTL
match_q <- a2opt_q==realdata[stg2,]$A2.NHTL
match_agree <- a2opt_agree==realdata[stg2,][agree,]$A2.NHTL

T_stg1_orig <- realdata$Y1.tilde
delta1_orig <- realdata$delta1
## set opt a2 as predicted by bart
T_stg1_b <- T_stg1_orig
T_stg1_b[stg2][match_b] <- T_stg1_orig[stg2][match_b] + T_stg2[match_b]  #entered stg2 & matched opt a2, keep Y2+Y1
delta1_b <- delta1_orig
delta1_b[stg2][match_b] <- delta2[match_b]  #entered stg2 & matched opt a2, keep delta2
delta1_b[stg2][!match_b] <- 0  #entered stg2 & unmatched opt a2, censor at Y1

## pick the predicted times that corresponding to the treatment observed
stg1_0 <- realdata$A1.NHTL==0
pred_b_stg1 <- res_b$A1.NHTL_1
pred_b_stg1[stg1_0,] <- res_b$A1.NHTL_0[stg1_0,]
pred_q_stg1 <- res_q$A1.NHTL_1
pred_q_stg1[stg1_0] <- res_q$A1.NHTL_0[stg1_0]

timepoint <- c(12,24,36)
#timepoint <- quantile(T_stg1_b, probs=seq(0,1,0.1))
ROC_b_stg1_b <- timeROC(T=T_stg1_b, delta=delta1_b, marker=-rowMeans(exp(pred_b_stg1)), cause=1, time=timepoint)
ROC_q_stg1_b <- timeROC(T=T_stg1_b, delta=delta1_b, marker=-pred_q_stg1, cause=1, time=timepoint)
comp_stg1_b <- cbind(bart=ROC_b_stg1_b$AUC, q=ROC_q_stg1_b$AUC)
(comp_stg1_b <- cbind(comp_stg1_b, diff=comp_stg1_b[,1]-comp_stg1_b[,2]))

## set opt a2 as predicted by q-learning
T_stg1_q <- T_stg1_orig
T_stg1_q[stg2][match_q] <- T_stg1_orig[stg2][match_q] + T_stg2[match_q]  #entered stg2 & matched opt a2, keep Y2+Y1
delta1_q <- delta1_orig
delta1_q[stg2][match_q] <- delta2[match_q]  #entered stg2 & matched opt a2, keep delta2
delta1_q[stg2][!match_q] <- 0  #entered stg2 & unmatched opt a2, censor at Y1

#timepoint <- quantile(T_stg1_q, probs=seq(0,1,0.1))
ROC_b_stg1_q <- timeROC(T=T_stg1_q, delta=delta1_q, marker=-rowMeans(exp(pred_b_stg1)), cause=1, time=timepoint)
ROC_q_stg1_q <- timeROC(T=T_stg1_q, delta=delta1_q, marker=-pred_q_stg1, cause=1, time=timepoint)
comp_stg1_q <- cbind(bart=ROC_b_stg1_q$AUC, q=ROC_q_stg1_q$AUC)
(comp_stg1_q <- cbind(comp_stg1_q, diff=comp_stg1_q[,1]-comp_stg1_q[,2]))

## set opt a2 as predicted by q-learning
T_stg1_agree <- T_stg1_orig
T_stg1_agree[stg2][match_agree] <- T_stg1_orig[stg2][match_agree] + T_stg2[match_agree]  #entered stg2 & matched opt a2, keep Y2+Y1
delta1_agree <- delta1_orig
delta1_agree[stg2][match_agree] <- delta2[match_agree]  #entered stg2 & matched opt a2, keep delta2
delta1_agree[stg2][!match_agree] <- 0  #entered stg2 & unmatched opt a2, censor at Y1

#timepoint <- quantile(T_stg1_agree, probs=seq(0,1,0.1))
ROC_b_stg1_agree <- timeROC(T=T_stg1_agree, delta=delta1_agree, marker=-rowMeans(exp(pred_b_stg1)), cause=1, time=timepoint)
ROC_q_stg1_agree <- timeROC(T=T_stg1_agree, delta=delta1_agree, marker=-pred_q_stg1, cause=1, time=timepoint)
comp_stg1_agree <- cbind(bart=ROC_b_stg1_agree$AUC, q=ROC_q_stg1_agree$AUC)
(comp_stg1_agree <- cbind(comp_stg1_agree, diff=comp_stg1_agree[,1]-comp_stg1_agree[,2]))
