## This is for testing the wrapper function of two stage DTR using BART
## The function created here will be transferred into a new R package

dtr1 <- function(x1=c("x1"), a1="a1", time1="y1",
                 x2=c("x1","x2"), a2="a2", stg2ind="eta", time2="y2", delta="delta",
                 data, newdata=NULL, opt=TRUE, mc.cores=8)
{
    eta <- data[, stg2ind]
    n <- nrow(data)
    stg2n <- sum(eta==1)
    stg2 <- which(eta==1)  #index of subj entered stg 2
    data.stg2 <- data[stg2, ]  #stg2 data
    newstg2n <- ifelse(is.null(newdata), 0, sum(newdata[,stg2ind]==1))
    a2option <- unique(data[, a2])  #all possible stg2 actions
    a2n <- length(a2option)  #total number of stg2 actions
    x.mat <- rbind(data.stg2,newdata[which(newdata[,stg2ind]==1),])[,x2,drop=F]  #x matrix for stg2
    x.test.stg2 <- data.frame(x.mat[rep(seq_len(nrow(x.mat)), each=a2n),], a2=rep(a2option,stg2n+newstg2n))  #stg2 actions for everyone
    fit2 <- mc.abart(x.train=data.stg2[,c(x2, a2)], times=data.stg2[, time2], delta=data.stg2[, delta], x.test=x.test.stg2, ndpost=1000, keepevery=5)
    burn <- (length(fit2$sigma) - fit2$ndpost)/ncol(fit2$sigma)
    sigma2 <- as.vector(fit2$sigma[-(1:burn),])
    if (!opt) {
        a2names <- paste0(a2,"_",a2option)
        for (i in 1:a2n) assign(a2names[i], fit2$yhat.test[,1:(stg2n+newstg2n)*a2n-(a2n-i)])
    }
    a2opt <- sapply(1:(stg2n+newstg2n)*a2n-1, function(column) {max.col(fit2$yhat.test[,column:(column+a2n-1)],'first')})  #opt col at stg2 for each iter
    a2.opt <- apply(a2opt, 2, function(x) a2option[x])  #opt action at stg2 in each iter
    yhat2optmean <- sapply(1:(stg2n+newstg2n)*a2n-1, function(column) {apply(fit2$yhat.test[,column:(column+a2n-1)],1,max)})  #opt y mean at stg2 for each iter
    newa2.opt <- a2.opt[,-(1:stg2n)]
    newyhat2optmean <- yhat2optmean[,-(1:stg2n)]
#######################################
##    validate est action & ymean  ##
#######################################
## a2true <- as.numeric(-0.9+0.6*x2>0)[stg2]
## a2optmean <- apply(a2.opt,2,
##                    function(x) names(sort(table(x), decreasing=TRUE)[1]))
## table(est=as.numeric(a2optmean),true=a2true)
## pdf(file="~/BART/a3/outputs/stg2.pdf",width=12,height=6)
## par(mfrow=c(1,2))
## plot(log(t2opt)[stg2],colMeans(yhat2optmean),pch=20,xlab="true opt stg2 time",ylab="est opt stg2 time mean")
## abline(a=0,b=1,col=2)
## t2bias <- colMeans(yhat2optmean)-log(t2opt)[stg2]
## boxplot(t2bias)
## abline(h=0,col=2)
## dev.off()
########################################################
    
    delta1 <- 1-(1-data[,delta])*(1-eta)
    a1option <- unique(data[, a1])  #all possible stage1 actions
    a1n <- length(a1option)  #total number of stage1 actions
    newn <- ifelse(is.null(newdata), 0, nrow(newdata))
    x.mat <- rbind(data,newdata)[,x1,drop=F]
    x.test.stg1 <- data.frame(x.mat[rep(seq_len(nrow(x.mat)), each=a1n),], a1=rep(a1option,n+newn))  #stg1 actions for everyone

    y2lst <- y1lst <- NULL
    ## rnormt <- function(n, a, b, mu, s = 1) {
  
    ##     F.a <- pnorm(a, mean = mu, sd = s)
    ##     F.b <- pnorm(b, mean = mu, sd = s)
  
    ##     u <- runif(n, min = F.a, max = F.b)
  
    ##     qnorm(u, mean = mu, sd = s)
        
    ## }

    ## library(truncnorm)
    ## library(MASS)
    stg1fit <- function(post){
        y2opt <- ifelse(post[-(1:(stg2n+1))]==data.stg2[, a2],
                 ifelse(data.stg2[, delta]==1, log(data.stg2[, time2]),
                        rtruncnorm(1,a=log(data.stg2[, time2]),mean=post[1:stg2n],sd=rep(post[stg2n+1],stg2n))),
                 mvrnorm(mu=post[1:stg2n],Sigma=diag(post[stg2n+1],nrow=stg2n)))
        ## aopt <- post[-(1:(stg2n+1))]==data.stg2[, a2]
        ## obs <- data.stg2[, delta]==1
        ## y2cen <- rnormt(stg2n, a=log(data.stg2[, time2]), b=Inf, mu=post[1:stg2n], s=post[stg2n+1]) 
        ## y2imp <- MASS::mvrnorm(mu=post[1:stg2n],Sigma=diag(post[stg2n+1],nrow=stg2n))
        ## y2opt <- log(data.stg2[,time2])
        ## y2opt[aopt & !obs] <- y2cen[aopt & !obs]
        ## y2opt[!aopt] <- y2imp[!aopt]
        time2opt <- rep(0, n)
        time2opt[stg2] <- exp(y2opt)
        y1new <- data[, time1]+time2opt
        fit1 <- abart(x.train=data[,c(x1,a1)], times=y1new, delta=delta1, x.test=x.test.stg1, ndpost=1)
        sigma1 <- tail(fit1$sigma, n=1)
        y1lst <- NULL
        if (!opt) {
            a1names <- paste0(a1,"_",a1option)
            y1lst <- vector(mode="list", length=a1n)
            names(y1lst) <- a1names
            for (i in 1:a1n) {
                assign(a1names[i], fit1$yhat.test[,1:(n+newn)*a1n-(a1n-i)])
                y1lst[[i]] <- eval(as.name(a1names[i]))
            }
        }
        a1opt <- sapply(1:(n+newn)*a1n-1, function(column) {max.col(fit1$yhat.test[,column:(column+a1n-1),drop=F],'first')})  #opt col at stg1
        a1.opt <- a1option[a1opt]  #opt action at stg1 in each iter
        yhat1optmean <- sapply(1:(n+newn)*a1n-1, function(column) {apply(fit1$yhat.test[,column:(column+a1n-1),drop=F],1,max)})  #opt y mean at stg1
        return(append(list(a1.opt=a1.opt,yhat1optmean=yhat1optmean,sigma1=sigma1),y1lst))  #post samples with opt action
    }
    
## a <- proc.time()
## cl <- makeCluster(getOption("cl.cores",cores))
## clusterEvalQ(cl, library(BART3))
## clusterExport(cl, c("n","data","stg2","a1option","a1n"))
## stg1res. <- parApply(cl, yhat2opt, 2, stg1fit)
## stopCluster(cl)
## partime <- proc.time()-a

    postsamp <- cbind(yhat2optmean[,1:stg2n],sigma2,a2.opt[,1:stg2n])
    cores <- mc.cores
    cl <- makeCluster(getOption("cl.cores",cores))
    registerDoParallel(cl)
    registerDoRNG(seed=1)
    res <- foreach(iter=1:1000, .combine=data.frame, .packages=c("BART3","truncnorm","MASS")) %dopar% stg1fit(postsamp[iter,])
    stopCluster(cl)

    ncolstg1res <- ncol(res)/1000
    a1names <- colnames(res)[4:ncolstg1res]
    res <- unname(as.matrix(res))
    
    a1.opt <- res[1:n,1:1000*ncolstg1res-(ncolstg1res-1)]
    yhat1optmean <- res[1:n,1:1000*ncolstg1res-(ncolstg1res-2)]
    newa1.opt <- res[-(1:n),1:1000*ncolstg1res-(ncolstg1res-1)]
    newyhat1optmean <- res[-(1:n),1:1000*ncolstg1res-(ncolstg1res-2)]
    sigma1 <- res[1,1:1000*ncolstg1res-(ncolstg1res-3)]

## a1optmean <- apply(a1.opt,2,
##                    function(x) names(sort(table(x), decreasing=TRUE)[1]))
## table(a1optmean)
    if (!opt) {
        y2lst <- vector(mode="list", length=a2n)
        names(y2lst) <- a2names
        for (i in 1:a2n) y2lst[[i]] <- t(eval(as.name(a2names[i])))
        y1lst <- vector(mode="list", length=a1n)
        names(y1lst) <- a1names
        for (i in 1:a1n) {
            assign(a1names[i], res[1:n,1:1000*ncolstg1res-(ncolstg1res-3-i)])
            y1lst[[i]] <- eval(as.name(a1names[i]))
        }
    }
    if (ncol(newa2.opt)==0) {final <- list(a2.opt=t(a2.opt),yhat2optmean=t(yhat2optmean),sigma2=sigma2,a1.opt=a1.opt,yhat1optmean=yhat1optmean,sigma1=sigma1)}
    else {final <- list(a2.opt=t(a2.opt[,1:stg2n]),yhat2optmean=t(yhat2optmean[,1:stg2n]),newa2.opt=t(newa2.opt),newyhat2optmean=t(newyhat2optmean),sigma2=sigma2,a1.opt=a1.opt,yhat1optmean=yhat1optmean,newa1.opt=newa1.opt,newyhat1optmean=newyhat1optmean,sigma1=sigma1)}

    final <- append(final,append(y2lst,y1lst))
    return(final)
}

## pdf(file="~/BART/a3/outputs/stg1.pdf",width=12,height=6)
## par(mfrow=c(1,2))
## plot(lgt.opt,colMeans(yhat1optmean),pch=20,xlab="true opt stg1 time",ylab="est opt stg1 time mean")
## abline(a=0,b=1,col=2)
## t1bias <- colMeans(yhat1optmean)-lgt.opt
## boxplot(t1bias)
## abline(h=0,col=2)
## dev.off()

dtr2 <- function(x1=c("x1"), a1="a1", time1="y1",
                 x2=c("x1","x2"), a2="a2", stg2ind="eta", time2="y2", delta="delta",
                 data, newdata=NULL, opt=TRUE, mc.cores=8)
{
    eta <- data[, stg2ind]
    n <- nrow(data)
    stg2n <- sum(eta==1)
    stg2 <- which(eta==1)  #index of subj entered stg 2
    data.stg2 <- data[stg2, ]  #stg2 data
    newstg2n <- ifelse(is.null(newdata), 0, sum(newdata[,stg2ind]==1))
    a2option <- unique(data[, a2])  #all possible stg2 actions
    a2n <- length(a2option)  #total number of stg2 actions
    x.mat <- rbind(data.stg2,newdata[which(newdata[,stg2ind]==1),])[,x2,drop=F]  #x matrix for stg2
    x.test.stg2 <- data.frame(x.mat[rep(seq_len(nrow(x.mat)), each=a2n),], a2=rep(a2option,stg2n+newstg2n))  #stg2 actions for everyone
    fit2 <- mc.abart(x.train=data.stg2[,c(x2, a2)], times=data.stg2[, time2], delta=data.stg2[, delta], x.test=x.test.stg2, ndpost=1000, keepevery=5)
    print("stg2reg finished")
    burn <- (length(fit2$sigma) - fit2$ndpost)/ncol(fit2$sigma)
    sigma2 <- as.vector(fit2$sigma[-(1:burn),])
    if (!opt) {
        a2names <- paste0(a2,"_",a2option)
        for (i in 1:a2n) assign(a2names[i], fit2$yhat.test[,1:(stg2n+newstg2n)*a2n-(a2n-i)])
    }
    a2opt <- sapply(1:(stg2n+newstg2n)*a2n-1, function(column) {max.col(fit2$yhat.test[,column:(column+a2n-1)],'first')})  #opt col at stg2 for each iter
    a2.opt <- apply(a2opt, 2, function(x) a2option[x])  #opt action at stg2 in each iter
    yhat2optmean <- sapply(1:(stg2n+newstg2n)*a2n-1, function(column) {apply(fit2$yhat.test[,column:(column+a2n-1)],1,max)})  #opt y mean at stg2 for each iter
    newa2.opt <- a2.opt[,-(1:stg2n)]
    newyhat2optmean <- yhat2optmean[,-(1:stg2n)]
    
    a1option <- unique(data[, a1])  #all possible stage1 actions
    a1n <- length(a1option)  #total number of stage1 actions
    newn <- ifelse(is.null(newdata), 0, nrow(newdata))
    x.mat <- rbind(data,newdata)[,x1,drop=F]
    x.test.stg1 <- data.frame(x.mat[rep(seq_len(nrow(x.mat)), each=a1n),], a1=rep(a1option,n+newn))  #stg1 actions for everyone

    y2lst <- y1lst <- NULL
    stg1fit <- function(post){
        a2opt.ind <- post[-(1:(stg2n+1))]==data.stg2[, a2]
        delta.stg2 <- ifelse(a2opt.ind, data.stg2[, delta], rep(1,stg2n))
        y2opt <- ifelse(a2opt.ind, log(data.stg2[, time2]),
                 mvrnorm(mu=post[1:stg2n],Sigma=diag(post[stg2n+1],nrow=stg2n)))
        time2opt <- rep(0, n)
        time2opt[stg2] <- exp(y2opt)
        y1new <- data[, time1]+time2opt
        delta1 <- 1-(1-data[,delta])*(1-eta)
        delta1[stg2] <- delta.stg2
        fit1 <- abart(x.train=data[,c(x1,a1)], times=y1new, delta=delta1, x.test=x.test.stg1, ndpost=1)
        sigma1 <- tail(fit1$sigma, n=1)
        y1lst <- NULL
        if (!opt) {
            a1names <- paste0(a1,"_",a1option)
            y1lst <- vector(mode="list", length=a1n)
            names(y1lst) <- a1names
            for (i in 1:a1n) {
                assign(a1names[i], fit1$yhat.test[,1:(n+newn)*a1n-(a1n-i)])
                y1lst[[i]] <- eval(as.name(a1names[i]))
            }
        }
        a1opt <- sapply(1:(n+newn)*a1n-1, function(column) {max.col(fit1$yhat.test[,column:(column+a1n-1),drop=F],'first')})  #opt col at stg1
        a1.opt <- a1option[a1opt]  #opt action at stg1 in each iter
        yhat1optmean <- sapply(1:(n+newn)*a1n-1, function(column) {apply(fit1$yhat.test[,column:(column+a1n-1),drop=F],1,max)})  #opt y mean at stg1
        return(append(list(a1.opt=a1.opt,yhat1optmean=yhat1optmean,sigma1=sigma1),y1lst))  #post samples with opt action
    }
print("parallel start")
    postsamp <- cbind(yhat2optmean[,1:stg2n],sigma2,a2.opt[,1:stg2n])
    cores <- mc.cores
    cl <- makeCluster(getOption("cl.cores",cores))
    registerDoParallel(cl)
    registerDoRNG(seed=1)
    res <- foreach(iter=1:1000, .combine=data.frame, .packages=c("BART3","MASS")) %dopar% stg1fit(postsamp[iter,])
    stopCluster(cl)
print("parallel finished")
    ncolstg1res <- ncol(res)/1000
    a1names <- colnames(res)[4:ncolstg1res]
    res <- unname(as.matrix(res))
    
    a1.opt <- res[1:n,1:1000*ncolstg1res-(ncolstg1res-1)]
    yhat1optmean <- res[1:n,1:1000*ncolstg1res-(ncolstg1res-2)]
    newa1.opt <- res[-(1:n),1:1000*ncolstg1res-(ncolstg1res-1)]
    newyhat1optmean <- res[-(1:n),1:1000*ncolstg1res-(ncolstg1res-2)]
    sigma1 <- res[1,1:1000*ncolstg1res-(ncolstg1res-3)]
    if (!opt) {
        y2lst <- vector(mode="list", length=a2n)
        names(y2lst) <- a2names
        for (i in 1:a2n) y2lst[[i]] <- t(eval(as.name(a2names[i])))
        y1lst <- vector(mode="list", length=a1n)
        names(y1lst) <- a1names
        for (i in 1:a1n) {
            assign(a1names[i], res[1:n,1:1000*ncolstg1res-(ncolstg1res-3-i)])
            y1lst[[i]] <- eval(as.name(a1names[i]))
        }
    }
    if (ncol(newa2.opt)==0) {final <- list(a2.opt=t(a2.opt),yhat2optmean=t(yhat2optmean),sigma2=sigma2,a1.opt=a1.opt,yhat1optmean=yhat1optmean,sigma1=sigma1)}
    else {final <- list(a2.opt=t(a2.opt[,1:stg2n]),yhat2optmean=t(yhat2optmean[,1:stg2n]),newa2.opt=t(newa2.opt),newyhat2optmean=t(newyhat2optmean),sigma2=sigma2,a1.opt=a1.opt,yhat1optmean=yhat1optmean,newa1.opt=newa1.opt,newyhat1optmean=newyhat1optmean,sigma1=sigma1)}

    final <- append(final,append(y2lst,y1lst))
    return(final)
}

dtr3 <- function(tdata,mc.cores=8){
    n <- nrow(tdata)
    stg2 <- which(tdata$eta==1)  #index of subj entered stg 2
    data.stg2 <- tdata[stg2, ]  #stg2 data
    stg2n <- nrow(data.stg2)
    a2option <- unique(tdata$a2)  #all possible stg2 actions
    a2n <- length(a2option)  #total number of stg2 actions
    x.mat <- data.stg2[,c("x1","x2")]  #x matrix for stg2
    x.test.stg2 <- data.frame(x.mat[rep(seq_len(nrow(x.mat)), each=a2n),], a2=rep(a2option,stg2n))  #stg2 actions for everyone
    fit2 <- mc.abart(x.train=data.stg2[,c("x1", "x2", "a2")], times=data.stg2$y2, delta=data.stg2$delta, x.test=x.test.stg2, ndpost=1000)
    burn <- (length(fit2$sigma) - fit2$ndpost)/ncol(fit2$sigma)
    sigma2 <- as.vector(fit2$sigma[-(1:burn),])  
    a2opt <- sapply(1:stg2n*a2n-1, function(column) {max.col(fit2$yhat.test[,column:(column+a2n-1)],'first')})  #opt col at stg2 for each iter
    a2.opt <- apply(a2opt, 2, function(x) sapply(x,function(x) a2option[x]))  #opt action at stg2 in each iter
    yhat2optmean <- sapply(1:stg2n*a2n-1, function(column) {apply(fit2$yhat.test[,column:(column+a2n-1)],1,max)})  #opt y mean at stg2 for each iter

    tdata$delta1 <- 1-(1-tdata$delta)*(1-tdata$eta)
    a1option <- unique(tdata$a1)  #all possible stage1 actions
    a1n <- length(a1option)  #total number of stage1 actions
    x.mat <- as.matrix(tdata[,"x1"])
    x.test.stg1 <- data.frame(x.mat[rep(seq_len(nrow(x.mat)), each=a1n),], a1=rep(a1option,n))  #stg1 actions for everyone

    stg1fit <- function(post){
        y2opt <- ifelse(post[-(1:(stg2n+1))]==data.stg2$a2,
                 ifelse(data.stg2$delta==1, log(data.stg2$y2),
                        truncnorm::rtruncnorm(1,a=log(data.stg2$y2),mean=post[1:stg2n],sd=rep(post[stg2n+1],stg2n))),
                 MASS::mvrnorm(mu=post[1:stg2n],Sigma=diag(post[stg2n+1],nrow=stg2n)))
        tdata$y2opt <- 0
        tdata[stg2,]$y2opt <- exp(y2opt)
        tdata$y1new <- tdata$y1+tdata$y2opt
        capture.output(fit1 <- abart(x.train=tdata[,c("x1","a1")], times=tdata$y1new, delta=tdata$delta1, x.test=x.test.stg1, ndpost=1))
        burn <- length(fit1$sigma) - nrow(fit1$yhat.train)
        sigma1 <- fit1$sigma[-(1:burn)]   
        a1opt <- sapply(1:n*a1n-1, function(column) {max.col(fit1$yhat.test[,column:(column+a1n-1),drop=F],'first')})  #opt col at stg1 for each iter
        a1.opt <- apply(as.matrix(a1opt), 2, function(x) sapply(x,function(x) a1option[x]))  #opt action at stg1 in each iter
        yhat1optmean <- sapply(1:n*a1n-1, function(column) {apply(fit1$yhat.test[,column:(column+a1n-1),drop=F],1,max)})  #opt y mean at stg2 for each iter
        return(list(a1.opt=a1.opt,yhat1optmean=yhat1optmean,sigma1=sigma1))  #post samples with opt action
    }

postsamp <- cbind(yhat2optmean,sigma2,a2.opt)
cores <- mc.cores
    
cl <- makeCluster(getOption("cl.cores",cores))
    registerDoParallel(cl)
    registerDoRNG(seed=1)
    res <- foreach(iter=1:1000, .export=c("tdata"), .combine=data.frame, .packages="BART3") %dopar% stg1fit(postsamp[iter,])
    stopCluster(cl)

    res <- as.matrix(res)
    res <- unname(res)

    a1.opt <- res[,1:1000*3-2]
    yhat1optmean <- res[,1:1000*3-1]
    sigma1 <- res[1,1:1000*3]

    return(list(a2.opt=t(a2.opt),yhat2optmean=t(yhat2optmean),sigma2=sigma2,a1.opt=a1.opt,yhat1optmean=yhat1optmean,sigma1=sigma1))
}
