realdata <- read.csv(file="~/BART/a3/data/all.comp.cases.csv")
library(penAFT)
library(survival)

varnames <- names(realdata)
varlist <- varnames[c(8,9,11,13,14,15,16)]
realdata[,varlist] <- lapply(realdata[,varlist], as.factor)

## x1 <- varnames[6:16]  # covariates in Stage 1
## x2 <- varnames[6:18]  # covariates in Stage 2

## convert factor var to dummy vars bc penaft only works with numeric vars
library(fastDummies)
realdata$newage <- as.factor(realdata$newage)
realdata$rec_match <- as.factor(realdata$rec_match)
varlist <- names(which(sapply(realdata, class)=="factor"))
datadum <- dummy_cols(realdata, select_columns=varlist, remove_selected_columns=T)
## newage: 1(<10yo) 0(10-39yo) 2(>=40yo)
## rec_match: 0(well-matched) 1(partially-matched) 2(mismatched)
## > levels(realdata$disease.risk.combined)
## [1] "Advanced"     "Early"        "Intermediate"
## > levels(realdata$graftype)
## [1] "Bone marrow"      "Cord blood"       "Peripheral blood"
datadum <- cbind(1, datadum)
xl.train <- datadum[, qr(datadum)$pivot[seq_len(qr(datadum)$rank)]]
dumdata <- xl.train[,-1]
names(dumdata) <- gsub(" ", ".", names(dumdata))

n <- nrow(dumdata)
dum.stg2 <- dumdata[dumdata$Int2==1,]
dum2n <- nrow(dum.stg2)
## Stg 2 regression
a2option <- unique(dumdata$A2.NHTL)
a2n <- length(a2option)
X2 <- dum.stg2[,!(colnames(dum.stg2) %in% c("dummy_id","Int2","Y2","delta2","Y1.tilde","delta1"))]  # covariates in Stage 2
crossX2 <- dum.stg2$A2.NHTL * X2[,-1]  #covariates cross stg2 trt
names(crossX2) <- paste0("A2.NHTL*", names(crossX2))
cov2 <- subset(X2, select=-c(A2.NHTL,A1.NHTL))  #main covariates
## two-way interactions among main covariates
inter2 <- inter2name <- NULL
for (i in 1:(ncol(cov2)-1)){
    if (i==1) inter2 <- cov2[,1] * cov2[,-1]
    else inter2 <- cbind(inter2, cov2[,1] * cov2[,-1])
    inter2name <- c(inter2name, paste0(names(cov2)[1], "*", names(cov2)[-1]))
    cov2 <- cov2[,-1]
}
names(inter2) <- inter2name
combX2 <- cbind(X2, crossX2, inter2)
## remove dependent columns
combX2 <- cbind(1, combX2)
xl.train <- combX2[, qr(combX2)$pivot[seq_len(qr(combX2)$rank)]]
combX2 <- xl.train[,-1]

x2names <- names(combX2)
#fit2 <- penAFT.cv(X = as.matrix(combX2), logY = log(dum.stg2$Y2), delta = dum.stg2$delta2, nlambda = 50, lambda.ratio.min = 0.01, penalty = "EN", alpha = 1)

## saveRDS(fit2, file="~/BART/a3/results/penaft2.RDS")
fit2 <- readRDS(file="~/BART/a3/results/penaft2.RDS")
beta2 <- penAFT.coef(fit2)$beta
qcov2 <- paste(x2names[beta2!=0], collapse="+")
fml2 <- as.formula(paste("Surv(Y2,delta2)", qcov2, sep="~"))
comb.stg2 <- cbind(Y2=dum.stg2$Y2, delta2=dum.stg2$delta2, combX2)
qstg2 <- survreg(fml2, data=comb.stg2, dist="lognormal")
A2.NHTL <- rep(a2option,dum2n)
X2.test <- X2[rep(seq_len(dum2n), each=a2n), -1]
crossX2.test <- A2.NHTL * X2.test
names(crossX2.test) <- paste0("A2.NHTL*", names(crossX2.test))
inter2.test <- inter2[rep(seq_len(dum2n), each=a2n), ]
combX2.test <- cbind(X2.test, crossX2.test, inter2.test)
## remove dependent columns
combX2.test <- cbind(1, combX2.test)
xl.train <- combX2.test[, qr(combX2.test)$pivot[seq_len(qr(combX2.test)$rank)]]
combX2.test <- xl.train[,-1]
pred.stg2 <- predict(qstg2, newdata=combX2.test)  #make predictions for each subj under all stg2 trt options

## compare and pick the optimal stg2 trt
a2names <- paste0("A2.NHTL_",a2option)
for (i in 1:a2n) assign(a2names[i], pred.stg2[1:dum2n*a2n-(a2n-i)])
a2opt <- sapply(1:dum2n*a2n-1, function(column) {which.max(pred.stg2[column:(column+a2n-1)])})  #opt col at stg2
a2.opt <- sapply(a2opt, function(x) a2option[x])  #opt action at stg2
yhat2optmean <- sapply(1:dum2n*a2n-1, function(column) {max(pred.stg2[column:(column+a2n-1)])})

## Stg 1 regression
a1option <- unique(dumdata$A1.NHTL)  #all possible stage1 actions
a1n <- length(a1option)  #total number of stage1 actions
data.stg1 <- dumdata
data.stg1$y2opt <- 0
data.stg1[data.stg1$Int2==1,]$y2opt <- yhat2optmean
data.stg1$Y1 <- data.stg1$Y1.tilde+data.stg1$y2opt

X1 <- data.stg1[,!(colnames(data.stg1) %in% c("dummy_id","Int2","Y2","delta2","Y1.tilde","delta1","Y1","y2opt","A2.NHTL","four.or.more.ISP_new","time.to.acute.GVHD"))]  # covariates in Stage 1
crossX1 <- data.stg1$A1.NHTL * X1[,colnames(X1)!="A1.NHTL"]  #covariates cross stg1 trt
names(crossX1) <- paste0("A1.NHTL*", names(crossX1))
cov1 <- subset(X1, select=-A1.NHTL)  #main covariates
## two-way interactions among main covariates
inter1 <- inter1name <- NULL
for (i in 1:(ncol(cov1)-1)){
    if (i==1) inter1 <- cov1[,1] * cov1[,-1]
    else inter1 <- cbind(inter1, cov1[,1] * cov1[,-1])
    inter1name <- c(inter1name, paste0(names(cov1)[1], "*", names(cov1)[-1]))
    cov1 <- cov1[,-1]
}
names(inter1) <- inter1name
combX1 <- cbind(X1, crossX1, inter1)
## remove dependent columns
combX1 <- cbind(1, combX1)
xl.train <- combX1[, qr(combX1)$pivot[seq_len(qr(combX1)$rank)]]
combX1 <- xl.train[,-1]

x1names <- names(combX1)
fit1 <- penAFT.cv(X = as.matrix(combX1), logY = log(data.stg1$Y1), delta = data.stg1$delta1, nlambda = 50, lambda.ratio.min = 0.01, penalty = "EN", alpha = 1)

saveRDS(fit1, file="~/BART/a3/results/penaft1.RDS")
## fit1 <- readRDS(file="~/BART/a3/results/penaft1.RDS")
## beta1 <- penAFT.coef(fit1)$beta
## qcov1 <- paste(x1names[beta1!=0], collapse="+")
## fml1 <- as.formula(paste("Surv(Y1,delta1)", qcov1, sep="~"))
## comb.stg1 <- cbind(Y1=data.stg1$Y1, delta1=data.stg1$delta1, combX1)
## qstg1 <- survreg(fml1, data=comb.stg1, dist="lognormal")
## A1.NHTL <- rep(a1option,n)
## X1.test <- X1[rep(seq_len(n), each=a1n),colnames(X1)!="A1.NHTL"]
## crossX1.test <- A1.NHTL * X1.test
## names(crossX1.test) <- paste0("A1.NHTL*", names(crossX1.test))
## inter1.test <- inter1[rep(seq_len(n), each=a1n), ]
## combX1.test <- cbind(X1.test, crossX1.test, inter1.test)
## ## remove dependent columns
## combX1.test <- cbind(1, combX1.test)
## xl.train <- combX1.test[, qr(combX1.test)$pivot[seq_len(qr(combX1.test)$rank)]]
## combX1.test <- xl.train[,-1]
## pred.stg1 <- predict(qstg1, newdata=combX1.test)  #make predictions for each subj under all stg1 trt options

## ## compare and pick the optimal stg1 trt
## a1names <- paste0("A1.NHTL_",a1option)
## for (i in 1:a1n) assign(a1names[i], pred.stg1[1:n*a1n-(a1n-i)])
## a1opt <- sapply(1:n*a1n-1, function(column) {which.max(pred.stg1[column:(column+a1n-1)])})  #opt col at stg1
## a1.opt <- sapply(a1opt, function(x) a1option[x])  #opt action at stg1
## yhat1optmean <- sapply(1:n*a1n-1, function(column) {max(pred.stg1[column:(column+a1n-1)])})


## qcov1 <- paste(paste("A1.NHTL*",x1),collapse="+")
## fml1 <- as.formula(paste("Surv(Y1,delta1)", qcov1, sep="~"))
## qstg1 <- survreg(fml1, data=data.stg1, dist="lognormal")
## stg1.est <- summary(qstg1)$table[,1]
## x.mat <- data.stg1[,x1]
## x.test.stg1 <- data.frame(x.mat[rep(seq_len(nrow(x.mat)), each=a1n),], A1.NHTL=rep(a1option,n))
## pred.stg1 <- predict(qstg1, newdata=x.test.stg1)
