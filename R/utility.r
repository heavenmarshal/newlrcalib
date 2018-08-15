evallrt <- function(delta,design,reps,cond,lrnuginit,nthread)
{
    ndes <- ncol(delta)
    N <- nrow(delta)
    if(ndes<=5) nthread <- 1
    if(ndes<nthread) nthread <- ndes
    logl0 <- -N*log(colMeans(delta^2))
    if(nthread<=1)
        uloglik <- apply(delta,2,maximizelikelihood.merr.lite,design,1,
                         reps=reps,conditioning=cond,lrnuginit=lrnuginit)
    else
    {
        cl <- makeCluster(nthread)
        uloglik <- tryCatch(parApply(cl,delta,2,maximizelikelihood.merr.lite,design,1,
                                     reps=reps,conditioning=cond,lrnuginit=lrnuginit),
                            finally=stopCluster(cl))
    }
    LRT <- -2*(logl0-uloglik)
}
calibfromLRT <- function(LRT,design,cands,alpha,reps,cond,bound1d)
{
    model <- maximizelikelihood(LRT,design,alpha,reps=reps,
                                conditioning=cond,bound1d=bound1d)
    new <- rbind(design,cands)
    LRT.hat <- rhodacepred.merr2(model,new)
    ests <- new[which.min(LRT.hat),]
    ret <- list(xopt=ests,model=model)
    return(ret)
}
calibtrace <- function(LRT,design,newLRT,knew,cands,alpha,reps,cond,bound1d)
{
    if(is.null(newLRT))
    {
        ret <- calibfromLRT(LRT,design,cands,alpha,reps,cond,bound1d)
        return(ret)
    }
    len <- length(newLRT)
    dd <- ncol(design)
    xopt <- matrix(nrow=len,ncol=dd)
    for(i in 1:len)
    {
        cdes <- rbind(design,knew[1:i,,drop=FALSE])
        cLRT <- c(LRT,newLRT[1:i])
        res <- calibfromLRT(cLRT,cdes,cands,alpha,reps,cond,bound1d)
        xopt[i,] <- res$xopt
    }
    ret <- list(xopt=xopt,model=res$m)
}

euclidean <- function(x,y)
{
    sqrt(sum((x-y)^2))
}
nearestPointIdx <- function(x,cand,exclude)
{
    edist <- apply(cand,1,euclidean,x)
    edist[exclude] <- Inf
    return(which.min(edist))
}
## the difference with calibfromLRT is this function do NOT merge
## design and cands
calibfromLRTset <- function(LRT,design,cands,alpha,reps,cond,bound1d)
{
    model <- maximizelikelihood(LRT,design,alpha,reps=reps,
                                conditioning=cond,bound1d=bound1d)
    LRT.hat <- rhodacepred.merr2(model,cands)
    optidx <- which.min(LRT.hat)
    ests <- cands[optidx,]
    ret <- list(xopt=ests,optidx=optidx,model=model)
    return(ret)
}

calibtraceset <- function(LRT,design,newLRT,knew,cands,alpha,reps,cond,bound1d)
{
    if(is.null(newLRT))
    {
        ret <- calibfromLRTset(LRT,design,cands,alpha,reps,cond,bound1d)
        return(ret)
    }
    len <- length(newLRT)
    dd <- ncol(design)
    xopt <- matrix(nrow=len,ncol=dd)
    optidx <- rep(NA,len)
    for(i in 1:len)
    {
        cdes <- rbind(design,knew[1:i,,drop=FALSE])
        cLRT <- c(LRT,newLRT[1:i])
        res <- calibfromLRTset(cLRT,cdes,cands,alpha,reps,cond,bound1d)
        xopt[i,] <- res$xopt
        optidx[i] <- res$optidx
    }
    ret <- list(xopt=xopt,optidx=optidx,model=res$model)
}

initQueue <- function(capacity)
{
    if(capacity<=0) error("capacity of a queue must be positive!")
    qv <- rep(NA,capacity)
    queue <- list(qv=qv,capacity=capacity,length=0)
    return(queue)
}
isFull <- function(queue)
{
    return(queue$capacity==queue$length)
}
## if is full automatically discard the end item
enQueue <- function(queue,x)
{
    if(queue$length==0 || queue$capacity==1)
    {
        queue$qv[1]=x
        queue$length=1
        return(queue)
    }
    qv <- queue$qv
    if(isFull(queue))
    {
        shiftidx <- 1:(queue$capacity-1)
        qv[shiftidx+1] <- qv[shiftidx]
        qv[1] <- x
        queue$qv <- qv
        return(queue)
    }
    len <- queue$length
    qv[2:(len+1)] <- qv[1:len]
    qv[1] <- x
    queue$qv <- qv
    queue$length <- len+1
    return(queue)
}
evalRatio <- function(queue,newx)
{
    if(queue$length==0) return(Inf)
    qv <- queue$qv
    denominator <- sum(abs(qv),na.rm=TRUE)
    extqv <- c(newx,qv)
    difeqv <- abs(diff(extqv))
    numerator <- sum(difeqv,na.rm=TRUE)
    return(numerator/denominator)
}

batchStop <- function(difqueue,info,difthres)
{
    len <- length(info)
    for(i in 1:len)
    {
        rr <- evalRatio(difqueue,info[i])
        if(isFull(difqueue) && rr < difthres)
        {
            ret <- list(difqueue=difqueue,stopp=TRUE)
            return(ret)
        }
        difqueue <- enQueue(difqueue,info[i])
    }
    ret <- list(difqueue=difqueue,stopp=FALSE)
}

batchCount <- function(difqueue,info,difthres)
{
    len <- length(info)
    count <- 0
    for(i in 1:len)
    {
        rr <- evalRatio(difqueue,info[i])
        if(isFull(difqueue) && rr < difthres)
            count <- count+1
        difqueue <- enQueue(difqueue,info[i])
    }
    ret <- list(difqueue=difqueue,count=count)
}
