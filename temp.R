n <- 400
rep <- 20
set.seed(4522)
seedset <- floor(rnorm(rep)*100000000)
i <- 20

set.seed(seedset[i])
x1 <- runif(n, 0.1, 1.29)  #baseline covariate
a1 <- rbinom(n, 1, 1/(1+exp(1-2*x1)))  #stg1 action
x2 <- runif(n, 0.9, 2)  #secondary covariate
a2 <- rbinom(n, 1, 1/(1+exp(2*x2-2.8)))  #stg2 action
e2 <- rnorm(n, 0, 0.3)
eta <- rbinom(n, 1, 0.6)  #stg2 indicator
t2 <- exp(4+1.1*x2-0.2*x2^3-0.1*x1+a2*(-0.9+0.6*x2) + e2)  #stg2 time
t2opt <- t2*exp(((-0.9+0.6*x2>0)*1-a2)*(-0.9+0.6*x2))  #stg2 opt time
t2opt_true <- t2/exp(e2)*exp(((-0.9+0.6*x2>0)*1-a2)*(-0.9+0.6*x2))  #stg2 opt time
t2a0 <- exp(4+1.1*x2-0.2*x2^3-0.1*x1)  #stg2 time with action 0
t2a1 <- exp(4+1.1*x2-0.2*x2^3-0.1*x1-0.9+0.6*x2)  #stg2 time with action 1
e1 <- rnorm(n, 0, 0.3)
t <- exp(6.3+1.5*x1-0.8*x1^4+a1*(0.1+0.1*x1) + e1)  #opt time
t_true <- exp(6.3+1.5*x1-0.8*x1^4+a1*(0.1+0.1*x1))  #opt time
t1 <- t-t2opt
total <- t1+t2
cen.init <- runif(n, 100, 2500)
cen <- cen.init+eta*t1
delta <- ifelse(eta, ifelse(cen<total, 0, 1), ifelse(cen<t, 0, 1))
y2 <- ifelse(eta, ifelse(delta, t2, cen-t1), 0)
y1 <- ifelse(eta, t1, pmin(t, cen))
data <- data.frame(delta, eta, t1, t2, t, t2opt, t2opt_true, t_true, t2a0, t2a1, y1, y2, cen, x1, x2, a1, a2)
