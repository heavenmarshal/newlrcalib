ei.batch.set <- function(model,batchsize,feasible,lwr=0,upr=1,reps=100)
{
    di <- ncol(model$design)
    xnew <- matrix(NA,ncol=di,nrow=batchsize)
    ynew <- rep(NA,batchsize)
    eivals <- rep(NA,batchsize)
    excidx <- NULL
    for(i in 1:batchsize)
    {
        temp <- ei(model,lwr=lwr,upr=upr,reps)
        txnew <- temp$bestx
        newidx <- nearestPointIdx(txnew,feasible,excidx)
        xnew[i,] <- feasible[newidx,]
        excidx <- c(excidx,newidx)
        eivals[i] <- temp$bestEI
        fmin <- min(model$response)
        ynew[i] <- calc.ei.detailed(xnew[i,],fmin,model)$yhat
        model <- updateGP(model,xnew[i,],ynew[i])
    }
    return(list(xnew=xnew,ynew=ynew,eivals=eivals,idxnew=excidx))
}
