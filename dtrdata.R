set.seed(9722)
n <- 1000

expit <- function(x) {1/(1+exp(-x))}
x1 <- runif(n, 0.1, 1.29)  #baseline covariate
a1 <- rbinom(n, 1, expit(2*x1-1))  #stg1 action
x2 <- runif(n, 0.9, 2)  #secondary covariate
a2 <- rbinom(n, 1, expit(-2*x2+2.8))  #stg2 action
e2 <- rnorm(n, 0, 0.3)
eta <- rbinom(n, 1, 0.6)  #stg2 indicator
t2 <- exp(4+0.3*x2+0.3*x1+a2*(-0.7+0.5*x2)+e2)  #stg2 time
a2true <- 1*(-0.7+0.5*x2>0)
t2opt <- t2*exp(((-0.7+0.5*x2>0)*1-a2)*(-0.7+0.5*x2))  # > range(log(t2opt)) 3.729356 5.890724
t2opt_true <- t2opt/exp(e2)
e1 <- rnorm(n, 0, 0.3)
t <- exp(6.3+0.7*x1-+a1*(0.1-0.2*x1) + e1)  #opt time
a1true <- 1*(0.1-0.2*x1>0)
t_true <- exp(6.3+0.7*x1+a1true*(0.1-0.2*x1))  # > range(log(t_true)) 6.450979 7.568638
t1 <- t-t2opt
total <- t1+t2
cen.init <- runif(n, 100, 2000)
cen <- cen.init+eta*t1
delta <- ifelse(eta, ifelse(cen<total, 0, 1), ifelse(cen<t, 0, 1))  # > mean(delta) 0.786
y2 <- ifelse(eta, ifelse(delta, t2, cen-t1), 0)
y1 <- ifelse(eta, t1, pmin(t, cen))
data <- data.frame(delta, eta, t1, t2, t, t2opt, t2opt_true, t_true, a2true, a1true, y1, y2, cen, x1, x2, a1, a2)

write.csv(data, file="~/BART/a3/codes/dtrdata.csv", row.names=FALSE)
