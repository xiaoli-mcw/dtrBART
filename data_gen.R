
n <- 800
m <- 400
rep <- 200
set.seed(4522)
seedset <- floor(rnorm(2*rep)*100000000)

data_lin <- function(n){
    expit <- function(x) {1/(1+exp(-x))}
    x1 <- runif(n, 0.1, 1.29)  #baseline covariate
    b1 <- rbinom(n, 1, 0.5)  #baseline binary covariate
    z1 <- rnorm(n, 10, 3)
    a1 <- rbinom(n, 1, expit(2*x1-1))  #stg1 action
    x2 <- runif(n, 0.9, 2)  #secondary covariate
    b2 <- rbinom(n, 1, 0.5)  #binary covariate
    z2 <- rnorm(n, 20, 4)
    a2 <- rbinom(n, 1, expit(-2*x2+2.8))  #stg2 action
    e2 <- rnorm(n, 0, 0.3)
    eta <- rbinom(n, 1, 0.6)  #stg2 indicator
    t2 <- exp(4+0.3*x2-0.6*x2*b2+b2+0.3*x1+0.4*b1-0.5*x1*b1+a2*(-0.7+0.5*x2-0.9*b2)+e2)  #stg2 time
    a2true <- 1*(-0.7+0.5*x2-0.9*b2>0)
    t2opt <- t2*exp(((-0.7+0.5*x2-0.9*b2>0)*1-a2)*(-0.7+0.5*x2-0.9*b2))  # > range(log(t2opt)) 3.729356 5.890724
    t2opt_true <- t2opt/exp(e2)
    e1 <- rnorm(n, 0, 0.3)
    t <- exp(6.3+0.7*x1+0.6*b1-0.8*x1*b1+a1*(0.1-0.2*x1+0.6*b1) + e1)  #opt time
    a1true <- 1*(0.1-0.2*x1+0.6*b1>0)
    t_true <- exp(6.3+0.7*x1+0.6*b1-0.8*x1*b1+a1true*(0.1-0.2*x1+0.6*b1))  # > range(log(t_true)) 6.450979 7.568638
    t1 <- t-t2opt
    total <- t1+t2
    cen.init <- runif(n, 100, 2000)
    cen <- cen.init+eta*t1
    delta <- ifelse(eta, ifelse(cen<total, 0, 1), ifelse(cen<t, 0, 1))  # > mean(delta) 0.786
    y2 <- ifelse(eta, ifelse(delta, t2, cen-t1), 0)
    y1 <- ifelse(eta, t1, pmin(t, cen))
    data <- data.frame(delta, eta, t1, t2, t, t2opt, t2opt_true, t_true, a2true, a1true, y1, y2, cen, x1, b1, z1, x2, b2, z2, a1, a2)
    return(data)
}

##############################################################
##  Data generation from Simoneau et al. (2019) ##
##  linear                                                  ##
##############################################################
## data_scn1 <- function(n){
##     expit <- function(x) {1/(1+exp(-x))}
##     x1 <- runif(n, 0.1, 1.29)  #baseline covariate
##     a1 <- rbinom(n, 1, expit(2*x1-1))  #stg1 action
##     x2 <- runif(n, 0.9, 2)  #secondary covariate
##     a2 <- rbinom(n, 1, expit(-2*x2+2.8))  #stg2 action
##     e2 <- rnorm(n, 0, 0.3)
##     eta <- rbinom(n, 1, 0.6)  #stg2 indicator
##     t2 <- exp(4+1.1*x2-0.2*x2^3-0.1*x1+a2*(-0.9+0.6*x2) + e2)  #stg2 time   > range(log(t2))  3.593261 5.795494
##     a2true <- 1*(-0.9+0.6*x2>0)
##     t2opt <- t2*exp(((-0.9+0.6*x2>0)*1-a2)*(-0.9+0.6*x2))  #stg2 opt time
##     t2opt_true <- exp(4+1.1*x2-0.2*x2^3-0.1*x1+a2*(-0.9+0.6*x2)+((-0.9+0.6*x2>0)*1-a2)*(-0.9+0.6*x2))
##                                         # > range(log(t2opt_true)) 4.737539 4.987389
##     e1 <- rnorm(n, 0, 0.3)
##     t <- exp(6.3+1.5*x1-0.8*x1^4+a1*(0.1-0.1*x1) + e1)  #opt time
##     a1true <- 1*(0.1-0.1*x1>0)
##     t_true <- exp(6.3+1.5*x1-0.8*x1^4+a1true*(0.1-0.1*x1))  # > range(log(t_true)) 5.998693 7.197007
##     t1 <- t-t2opt
##     total <- t1+t2
##     cen.init <- runif(n, 100, 2500)
##     cen <- cen.init+eta*t1
##     delta <- ifelse(eta, ifelse(cen<total, 0, 1), ifelse(cen<t, 0, 1))  # mean(delta)  0.845
##     y2 <- ifelse(eta, ifelse(delta, t2, cen-t1), 0)
##     y1 <- ifelse(eta, t1, pmin(t, cen))
##     data <- data.frame(delta, eta, t1, t2, t, t2opt, t2opt_true, t_true, y1, y2, cen, x1, x2, a1, a2)
##     return(data)
## }
############ training sets ######################
i <- 1
j <- 0
while (i <= rep){
    set.seed(seedset[i])
    data <- data_lin(n)
    if(min(data$t1)<=0) {print(i);seedset <- seedset[-i]; j <- j+1; next}
    write.csv(data,file=paste0("~/BART/a3/train/data_lin_",i,".csv"),row.names=FALSE)
    i <- i+1
}
set.seed(51222)
newdata <- data_lin(m)
while(min(newdata$t1)<=0) newdata <- data_lin(m)
write.csv(newdata,file=paste0("~/BART/a3/test/data_lin.csv"),row.names=FALSE)
############ OLD training sets ######################
## i <- 1
## j <- 0
## while (i <= rep){
##     set.seed(seedset[i])
##     data <- data_scn1(n)
##     if(min(data$t1)<=0) {print(i);seedset <- seedset[-i]; j <- j+1; next}
##     write.csv(data,file=paste0("~/BART/a3/train/data_scn1_",i,".csv"),row.names=FALSE)
##     i <- i+1
## }
## set.seed(51222)
## newdata <- data_scn1(m)
## while(min(newdata$t1)<=0) newdata <- data_scn1(m)
## write.csv(newdata,file=paste0("~/BART/a3/test/data_scn1.csv"),row.names=FALSE)
##############################################################
data_nln <- function(n){
    expit <- function(x) {1/(1+exp(-x))}
    x1 <- runif(n, 0.1, 1.29)  #baseline covariate
    b1 <- rbinom(n, 1, 0.5)  #baseline binary covariate
    z1 <- rnorm(n, 10, 3)
    a1 <- rbinom(n, 1, expit(2*x1-1))  #stg1 action
    x2 <- runif(n, 0.9, 2)  #secondary covariate
    b2 <- rbinom(n, 1, 0.5)  #binary covariate
    z2 <- rnorm(n, 20, 4)
    a2 <- rbinom(n, 1, expit(-2*x2+2.8))  #stg2 action
    e2 <- rnorm(n, 0, 0.1)
    eta <- rbinom(n, 1, 0.6)  #stg2 indicator
    t2 <- exp(4+cos(x2^3)-0.4*(x2*b2+0.5)^2-0.1*x1-sin(pi*x1*b1)+a2*(0.7*x2^2-1) + e2)  #stg2 time
    a2true <- 1*(0.7*x2^2-1>0)
    t2opt <- exp(4+cos(x2^3)-0.4*(x2*b2+0.5)^2-0.1*x1-sin(pi*x1*b1)+a2true*(0.7*x2^2-1) + e2)  #stg2 opt time
    t2opt_true <- t2opt/exp(e2)
    e1 <- rnorm(n, 0, 0.1)
    t <- exp(7.4+sin(x1^2)+x1^4+x1*b1+a1*(0.1-0.2*x1^3) + e1)  #opt time
    a1true <- 1*(0.1-0.2*x1^3>0)
    t_true <- exp(7.4+sin(x1^2)+x1^4+x1*b1+a1true*(0.1-0.2*x1^3))
    t1 <- t-t2opt
    total <- t1+t2
    cen.init <- runif(n, 400, 5000)
    cen <- cen.init+eta*t1
    delta <- ifelse(eta, ifelse(cen<total, 0, 1), ifelse(cen<t, 0, 1))
    y2 <- ifelse(eta, ifelse(delta, t2, cen-t1), 0)
    y1 <- ifelse(eta, t1, pmin(t, cen))
    data <- data.frame(delta, eta, t1, t2, t, t2opt, t2opt_true, t_true, a2true, a1true, y1, y2, cen, x1, b1, z1, x2, b2, z2, a1, a2)
    return(data)
}
############ training sets ######################
i <- 1
j <- 0
while (i <= rep){
    set.seed(seedset[i])
    data <- data_nln(n)
    if(min(data$t1)<=0) {print(i);seedset <- seedset[-i]; j <- j+1; next}
    write.csv(data,file=paste0("~/BART/a3/train/data_nln_",i,".csv"),row.names=FALSE)
    i <- i+1
}
set.seed(62222)
newdata <- data_nln(m)
while(min(newdata$t1)<=0) newdata <- data_nln(m)
write.csv(newdata,file=paste0("~/BART/a3/test/data_nln.csv"),row.names=FALSE)
