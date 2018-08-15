ei.batch <- function(model,batchsize,lwr=0,upr=1,reps=100)
{
    di <- ncol(model$design)
    xnew <- matrix(NA,ncol=di,nrow=batchsize)
    ynew <- rep(NA,batchsize)
    eivals <- rep(NA,batchsize)
    for(i in 1:batchsize)
    {
        temp <- ei(model,lwr=lwr,upr=upr,reps)
        xnew[i,] <- temp$bestx
        eivals[i] <- temp$bestEI
        fmin <- min(model$response)
        ynew[i] <- calc.ei.detailed(xnew[i,],fmin,model)$yhat
        model <- updateGP(model,xnew[i,],ynew[i])
    }
    return(list(xnew=xnew,ynew=ynew,eivals=eivals))
}
ei <- function(model,lwr=0,upr=1,reps=100)
{
    nparms <- ncol(model$design)+1
    bestx <- NULL
    bestEI <- -1e6
    fmin <- min(model$response)
    x.mat <- runif.sobol(reps,nparms-1)*(upr-lwr)+lwr
    for(i in 1:reps)
    {
        x.init <- x.mat[i,]
        maximized <- optim(x.init,calc.ei,method="L-BFGS-B",control=list(ndeps=rep(1e-20,nparms-1)),
                           lower=rep(lwr,nparms-1), upper=rep(upr,nparms-1), fmin=fmin,model=model)
        if(maximized$convergence!=0) warning("Optim did not converge!\n converge=",maximized$convergence, " ",maximized$message)
        value <- -maximized$value
        if(value>bestEI)
        {
            bestx <- maximized$par
            bestEI <- value
        }
    }
    return(list(bestx=bestx,bestEI=bestEI))
}

calc.ei <- function(x,fmin,model)
{
    N <- nrow(model$design)
    p <- ncol(model$design)
    out <- .C("calcEI",as.double(x),as.double(model$design),
              as.double(model$solres), as.double(model$maximized),
              as.double(model$betahat),as.double(model$sigma2hat),
              as.double(model$sigma2epsilonhat), as.double(model$alpha),
              as.double(model$Einv), as.double(fmin), as.integer(N),
              as.integer(p), yhat = double(1), mse = double(1),
              nei = double(1),PACKAGE="newlrcalib")
    return(out$nei)
}

calc.ei.detailed <- function(x,fmin,model)
{
    N <- nrow(model$design)
    p <- ncol(model$design)
    out <- .C("calcEI",as.double(x),as.double(model$design),
              as.double(model$solres), as.double(model$maximized),
              as.double(model$betahat),as.double(model$sigma2hat),
              as.double(model$sigma2epsilonhat), as.double(model$alpha),
              as.double(model$Einv), as.double(fmin), as.integer(N),
              as.integer(p), yhat = double(1), mse = double(1),
              nei = double(1),PACKAGE="newlrcalib")
    return(list(EI=out$nei,yhat=out$yhat,shat=sqrt(out$mse)))
}
