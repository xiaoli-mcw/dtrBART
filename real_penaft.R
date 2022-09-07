realdata <- read.csv(file="~/BART/a3/data/all.comp.cases.csv")
library(penAFT)

varnames <- names(realdata)
varlist <- varnames[c(8,9,11,13,14,15,16)]
realdata[,varlist] <- lapply(realdata[,varlist], as.factor)

x1 <- varnames[6:16]  # covariates in Stage 1
x2 <- varnames[6:18]  # covariates in Stage 2

n <- nrow(realdata)
real.stg2 <- realdata[realdata$Int2==1,]
real2n <- nrow(real.stg2)
## Stg 2 regression
a2option <- unique(realdata$A2.NHTL)
a2n <- length(a2option)
X2 <- real.stg2[,colnames(real.stg2) %in% c(x2,"A1.NHTL","A2.NHTL")]
crossX1 <- real.stg2$A2.NHTL * X1  #error

fit.main <- penAFT.cv(X = X2, logY = log(real.stg2$Y2), delta = real.stg2$delta2, nlambda = 50, lambda.ratio.min = 0.01, penalty = "EN", alpha = 1)

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



set.seed(1)
genData <- genSurvData(n = 50, p = 50, s = 10, mag = 2, cens.quant = 0.6)
X <- genData$X
logY <- genData$logY
delta <- genData$status

fit.en <- penAFT.cv(X = X, logY = logY, delta = delta, nlambda = 50, lambda.ratio.min = 0.01, penalty = "EN", alpha = 1)
coef.en.10 <- penAFT.coef(fit.old, lambda = fit.old$lambda[10])
coef.en <- penAFT.coef(fit.en)
which(coef.en$beta>0)
