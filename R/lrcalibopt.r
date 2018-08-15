lrdev <- function(x,model)
{
    x <- matrix(x,nrow=1)
    yhat <- drop(rhodacepred.merr2(model,x))
}
calibfromLRTopt <- function(LRT,design,alpha,reps,cond,bound1d,gactl=list())
{
    ndim <- ncol(design)
    gactl.default <- list(pop.size=1000, max.generations = 100,
                          wait.generations = 10)
    remnames <- setdiff(names(gactl.default),names(gactl))
    gactl <- c(gactl,gactl.default[remnames])
    model <- maximizelikelihood(LRT,design,alpha,reps=reps,
                                conditioning=cond,bound1d=bound1d)
    opt <- genoud(lrdev,ndim,pop.size=gactl$pop.size,
                  max.generations=gactl$max.generations,
                  wait.generations=gactl$wait.generations,
                  Domains=cbind(rep(0,ndim),rep(1,ndim)),
                  boundary.enforcement=2,print.level=0,
                  model=model)
    xopt <- opt$par
    ret <- list(xopt=xopt,model=model)
    return(ret)
}
calibfromLRTliteopt <- function(model,ndim,gactl=list())
{
    gactl.default <- list(pop.size=1000, max.generations = 100,
                          wait.generations = 10)
    remnames <- setdiff(names(gactl.default),names(gactl))
    gactl <- c(gactl,gactl.default[remnames])
    opt <- genoud(lrdev,ndim,pop.size=gactl$pop.size,
                  max.generations=gactl$max.generations,
                  wait.generations=gactl$wait.generations,
                  Domains=cbind(rep(0,ndim),rep(1,ndim)),
                  boundary.enforcement=2,print.level=0,
                  model=model)
    xopt <- opt$par
    ret <- list(xopt=xopt,model=model)
}
lrcalibseqopt <- function(xi,yi,yobs,timepoints,nadd,nbatch,
                          func,...,alpha=2,lwr=0,upr=1,repsext=20,repsei=100,
                          cond=1e-10,lrnuginit=1e-14,
                          nthread=4,bound1d=FALSE,gactl=list())
{
    xi <- as.matrix(xi)
    timepoints <- as.matrix(timepoints)
    delta <- yobs-yi
    LRT <- evallrt(delta,timepoints,repsext,cond,lrnuginit,nthread)
    nrem <- nadd
    maxinfo <- NULL
    thres <- 0
    iter <- 1
    while(nrem>0)
    {
        model <- maximizelikelihood(LRT,xi,alpha,reps=repsext,
                                    conditioning=cond,bound1d=bound1d)
        if(nrem<nbatch) nbatch <- nrem
        eiobj <- ei.batch(model,nbatch,lwr=lwr,upr=upr,reps=repsei)
        knew <- eiobj$xnew
        mm <- eiobj$eivals
        maxinfo <- c(maxinfo,mm)
        ynew <- apply(knew,1,func,...)
        deltanew <- as.matrix(yobs-ynew)
        LRTnew <- evallrt(deltanew,timepoints,repsext,cond,lrnuginit,nthread)
        xi <- rbind(xi,knew)
        yi <- cbind(yi,ynew)
        LRT <- c(LRT,LRTnew)
        nrem <- nrem-nbatch
        iter <- iter+1
    }
    ret <- calibfromLRTopt(LRT,xi,alpha,repsext,cond,bound1d)
    out <- list(xx=xi,yy=yi,xopt=ret$xopt,maxinfo=maxinfo)
    return(out)
}

lrcalibseqopttr <- function(xi,yi,yobs,timepoints,nadd,nbatch,
                            func,...,alpha=2,lwr=0,upr=1,repsext=20,repsei=100,
                            cond=1e-10,lrnuginit=1e-14,
                            nthread=4,bound1d=FALSE,gactl=list())
{
    xi <- as.matrix(xi)
    ndim <- ncol(xi)
    timepoints <- as.matrix(timepoints)
    delta <- yobs-yi
    LRT <- evallrt(delta,timepoints,repsext,cond,lrnuginit,nthread)
    nrem <- nadd
    maxinfo <- NULL
    xoptr <- NULL
    thres <- 0
    iter <- 1
    while(nrem>0)
    {
        model <- maximizelikelihood(LRT,xi,alpha,reps=repsext,
                                    conditioning=cond,bound1d=bound1d)
        cxopt <- calibfromLRTliteopt(model,ndim,gactl)
        xoptr <- rbind(xoptr,cxopt$xopt)
        if(nrem<nbatch) nbatch <- nrem
        eiobj <- ei.batch(model,nbatch,lwr=lwr,upr=upr,reps=repsei)
        knew <- eiobj$xnew
        mm <- eiobj$eivals
        maxinfo <- c(maxinfo,mm)
        ynew <- apply(knew,1,func,...)
        deltanew <- as.matrix(yobs-ynew)
        LRTnew <- evallrt(deltanew,timepoints,repsext,cond,lrnuginit,nthread)
        xi <- rbind(xi,knew)
        yi <- cbind(yi,ynew)
        LRT <- c(LRT,LRTnew)
        nrem <- nrem-nbatch
        iter <- iter+1
    }
    ret <- calibfromLRTopt(LRT,xi,alpha,repsext,cond,bound1d)
    xoptr <- rbind(xoptr,ret$xopt)
    out <- list(xx=xi,yy=yi,xopt=ret$xopt,xoptr=xoptr,maxinfo=maxinfo)
    return(out)
}
