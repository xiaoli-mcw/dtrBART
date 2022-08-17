bmabart <- function(formula, family=gaussian, data, prior.df=Inf, prior.p=0.5,
                    x.test=matrix(0,0,0), sparse=FALSE, ndpost=1000L,
                    nskip=100L, mc.cores=2L, seed=999){
    n <- nrow(data)
    x.train <- data[,labels(terms(formula))]
    y.train <- data[,all.vars(formula)[1]]
    fitarm <- bayesglm(formula,family=family,data=data,prior.df=prior.df)
    fitpb <- mc.gbart(x.train=x.train,y.train=y.train,x.test=x.test,type="pbart",sparse=sparse,ndpost=ndpost,nskip=nskip,mc.cores=mc.cores)

    Y <- t(matrix(y.train,nrow=n,ncol=ndpost))

    postb <- coef(sim(fitarm, n.sims=ndpost))
    arm.prob.train <- apply(postb, 1, FUN=function(x) pnorm(matrix(x,nrow=1)%*%rbind(1,t(x.train))))
    arm.prob.train <- t(arm.prob.train)
    logpdf.arm <- dbinom(Y,1,arm.prob.train,TRUE)
    minlogpdf.arm <- t(matrix(apply(logpdf.arm,2,min),nrow=n,ncol=ndpost))
    logCPO.arm <- log(ndpost)+minlogpdf.arm[1,]-log(apply(exp(minlogpdf.arm-logpdf.arm),2,sum))
    LPML.arm <- sum(logCPO.arm)

    LPML.pb <- fitpb$LPML
    
    diffLPML <- LPML.arm-LPML.pb
    p.arm <- prior.p/(exp(-diffLPML)*(1-prior.p)+prior.p)

    set.seed(seed)
    index <- rbinom(ndpost,1,p.arm)
    prob.train <- rbind(arm.prob.train[index==1,],fitpb$prob.train[index==0,])
    res <- list(prior.GLM=prior.p, prior.BART=1-prior.p, post.GLM=p.arm, post.BART=1-p.arm, prob.train=prob.train)
    res$prob.train.mean <- apply(prob.train, 2, mean)

    if (length(x.test)>0) {
        arm.prob.test <- apply(postb, 1, FUN=function(x) pnorm(matrix(x,nrow=1)%*%rbind(1,t(x.test))))
        arm.prob.test <- t(arm.prob.test)
        res$prob.test <- rbind(arm.prob.test[index==1,],fitpb$prob.test[index==0,])
        res$prob.test.mean <- apply(res$prob.test, 2, mean)
    }
    return(res)
}
