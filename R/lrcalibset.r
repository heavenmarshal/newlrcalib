lrcalibseqset <- function(xi,yi,yobs,timepoints,nadd,nbatch,feasible,cands,
                          fearesp,alpha=2,lwr=0,upr=1,repsext=20,repsei=100,
                          cond=1e-10,lrnuginit=1e-14,relthres=0,nthread=4,bound1d=FALSE)
{
    xi <- as.matrix(xi)
    cands <- as.matrix(cands)
    feasible <- as.matrix(feasible)
    timepoints <- as.matrix(timepoints)
    delta <- yobs-yi
    LRT <- evallrt(delta,timepoints,repsext,cond,lrnuginit,nthread)
    nrem <- nadd
    maxinfo <- NULL
    iter <- 1
    thres <- 0
    while(nrem>0)
    {
        model <- maximizelikelihood(LRT,xi,alpha,reps=repsext,
                                    conditioning=cond,bound1d=bound1d)
        if(nrem<nbatch) nbatch <- nrem
        newret <- ei.batch.set(model,nbatch,feasible,lwr=lwr,upr=upr,reps=repsei)
        knew <- newret$xnew
        idxnew <- newret$idxnew
        mm <- newret$eivals
        maxinfo <- c(maxinfo,mm)
        if(iter == 1) thres <- max(mm)*relthres
        if(min(mm)<thres)
        {
            idx <- min(which(mm<thres))-1
            if(idx<=0) break
            knew <- knew[1:idx,,drop=FALSE]
            idxnew <- idxnew[1:idx]
            ynew <- fearesp[,idxnew,drop=FALSE]
            deltanew <- as.matrix(yobs-ynew)
            LRTnew <- evallrt(deltanew,timepoints,repsext,cond,lrnuginit,nthread)
            xi <- rbind(xi,knew)
            yi <- cbind(yi,ynew)
            LRT <- c(LRT,LRTnew)
            break
        }
        ynew <- fearesp[,idxnew,drop=FALSE]
        deltanew <- as.matrix(yobs-ynew)
        LRTnew <- evallrt(deltanew,timepoints,repsext,cond,lrnuginit,nthread)
        xi <- rbind(xi,knew)
        yi <- cbind(yi,ynew)
        LRT <- c(LRT,LRTnew)
        feasible <- feasible[-idxnew,,drop=FALSE]
        fearesp <- fearesp[,-idxnew,drop=FALSE]
        nrem <- nrem-nbatch
        iter <- iter+1
    }
    ret <- calibfromLRTset(LRT,xi,cands,alpha,repsext,cond,bound1d)
    out <- list(xx=xi,yy=yi,xopt=ret$xopt,maxinfo=maxinfo)
    return(out)
}
lrcalibseqsettr <- function(xi,yi,yobs,timepoints,nadd,nbatch,feasible,cands,
                            fearesp,candsresp,alpha=2,lwr=0,upr=1,repsext=20,repsei=100,
                            cond=1e-10,lrnuginit=1e-14,relthres=0,difthres=0, difstep=1,
                            dxithres=0, dxistep=1, dxihittime=1, nthread=4,bound1d=FALSE)
{
    xi <- as.matrix(xi)
    feasible <- as.matrix(feasible)
    cands <- as.matrix(cands)
    feasible <- as.matrix(feasible)
    timepoints <- as.matrix(timepoints)
    delta <- yobs-yi
    tlen <- length(yobs)
    LRT <- evallrt(delta,timepoints,repsext,cond,lrnuginit,nthread)
    nrem <- nadd
    knew <- LRTnew<- NULL
    ret <- calibtraceset(LRT,xi,LRTnew,knew,cands,alpha,repsext,cond,bound1d)
    xoptr <- ret$xopt
    maxinfo <- NULL
    thres <- 0
    iter <- 1
    difqueue <- initQueue(difstep)
    dxiqueue <- initQueue(dxistep)
    val <- var(yobs)*(tlen-1)
    stopflag <- 0
    hitcount <- 0
    while(nrem>0)
    {
        model <- ret$model
        if(nrem<nbatch) nbatch <- nrem
        newret <- ei.batch.set(model,nbatch,feasible,lwr=lwr,upr=upr,reps=repsei)
        mm <- newret$eivals
        maxinfo <- c(maxinfo,mm)
        knew <- newret$xnew
        idxnew <- newret$idxnew
        if(iter == 1) thres <- relthres*max(mm)
        if(min(mm)<thres)
        {
            stopflag <- 1
            idx <- min(which(mm<thres))-1
            if(idx<=0) break
            knew <- knew[1:idx,,drop=FALSE]
            idxnew <- idxnew[1:idx]
            ynew <- fearesp[,idxnew,drop=FALSE]
            deltanew <- as.matrix(yobs-ynew)
            LRTnew <- evallrt(deltanew,timepoints,repsext,cond,lrnuginit,nthread)
            ret <- calibtraceset(LRT,xi,LRTnew,knew,cands,alpha,repsext,cond,bound1d)
            xoptr <- rbind(xoptr,ret$xopt)
            xi <- rbind(xi,knew)
            yi <- cbind(yi,ynew)
            LRT <- c(LRT,LRTnew)
            break
        }
        bres <- batchStop(difqueue,mm,difthres)
        if(bres$stopp)
        {
            stopflag <- 2
            break
        }
        difqueue <- bres$difqueue

        ynew <- fearesp[,idxnew,drop=FALSE]
        deltanew <- as.matrix(yobs-ynew)
        LRTnew <- evallrt(deltanew,timepoints,repsext,cond,lrnuginit,nthread)
        ret <- calibtraceset(LRT,xi,LRTnew,knew,cands,alpha,repsext,cond,bound1d)
        xoptr <- rbind(xoptr,ret$xopt)
        xi <- rbind(xi,knew)
        yi <- cbind(yi,ynew)
        LRT <- c(LRT,LRTnew)

        yopts <- candsresp[,ret$optidx,drop=FALSE]
        cdevs <- log(apply((yobs-yopts)^2,2,sum)/val)
        dxires <- batchCount(dxiqueue,cdevs,dxithres)
        hitcount <- hitcount + dxires$count
        if(hitcount >= dxihittime)
        {
            stopflag <- 3
            break
        }
        dxiqueue <- dxires$difqueue
        feasible <- feasible[-idxnew,,drop=FALSE]
        fearesp <- fearesp[,-idxnew,drop=FALSE]
        nrem <- nrem-nbatch
        iter <- iter+1
    }
    trlen <- nrow(xoptr)
    xopt <- xoptr[trlen,]
    out <- list(xx=xi,yy=yi,xopt=xopt,xoptr=xoptr,maxinfo=maxinfo,stopflag=stopflag)
    return(out)
}
