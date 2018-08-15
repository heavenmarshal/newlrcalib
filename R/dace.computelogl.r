corrmat <- function(rho, design, alpha, sigma2,nug)
{
    n <- nrow(design)
    p <- ncol(design)
    out <- .C("corrmat",as.double(design), as.double(rho),
              as.double(alpha), as.double(sigma2),
              as.double(nug), as.integer(n),
              as.integer(p), ans = double(n*n),
              PACKAGE="newlrcalib")
    corr <- matrix(out$ans,nrow=n)
}
crosscorrmat <- function(rho,design1,design2,alpha)
{
    n1 <- nrow(design1)
    n2 <- nrow(design2)
    p <- ncol(design1)
    out <- .C("crosscorrmat",as.double(design1), as.double(design2),
              as.double(rho), as.double(alpha), as.integer(n1),
              as.integer(n2), as.integer(p), ans = double(n1*n2),
              PACKAGE="newlrcalib")
    corr <- matrix(out$ans,nrow=n1)
}
computelogl.1D <- function(rho,design,response,alpha,N,conditioning)
{
    out <- .C("computelogl",as.double(rho), as.double(design),
              as.double(response), as.double(alpha),
              as.double(conditioning), as.integer(N),
              as.integer(1), ans = double(1),
              PACKAGE = "newlrcalib")
    return(out$ans)
}
computelogl.pdim <- function(theta.tilde,design,response,alpha,N,conditioning)
{
    rho <- exp(-1/4*exp(theta.tilde))
    p <- ncol(design)
    out <- .C("computelogl",as.double(rho), as.double(design),
              as.double(response), as.double(alpha),
              as.double(conditioning), as.integer(N),
              as.integer(p), ans = double(1),
              PACKAGE = "newlrcalib")
    return(out$ans)
}
maximizelikelihood <- function(response,design,alpha=2,reps=100,
                               conditioning=.Machine$double.eps,rho.ini=1,
                               bound1d=FALSE)
{
    design <- as.matrix(design)
    N <- length(response)
    pdims <- ncol(design)
    if(pdims==1) #1D CASE
    {
        upbound <- if(bound1d) exp(-0.25/darg(NULL,design)$max) else (1-.Machine$double.eps)
        maximized <- optimize(computelogl.1D,interval=c(.Machine$double.eps,upbound),
                              design=design,response=response,
                              alpha=alpha,N=N,conditioning=conditioning)
        bestlogl <- -maximized$objective
        maximized <- maximized$minimum
        maximized <- exp(-1/4*exp(maximized))
    }
    else # P-DIM CASE
    {
        bestlogl <- -1e16
        bestrho <- NULL
        rho.mat <- runif.sobol(reps,pdims)*rho.ini
        for(i in 1:reps)
        {
            rho.init <- rho.mat[i,] # 0<rho<1, so init to 0.5 is right in the middle.
            theta.tilde.init <- log(-4*log(rho.init)) # optim will optimize over theta.tilde instead of rho directly.
            maximized <- optim(theta.tilde.init,computelogl.pdim,method="Nelder-Mead",
                               control=list(alpha=1.0,beta=0.5,gamma=2.0),design=design,
                               response=response,alpha=alpha,N=N,conditioning=conditioning)
            if(maximized$convergence!=0) warning("Optim did not converge!\n")
            value <- -maximized$value
            if(value>bestlogl)
            {
                bestrho <- maximized$par
                bestlogl <- value
            }
        }
        maximized <- bestrho
        maximized <- exp(-1/4*exp(maximized)) #convert back to rho-scale, since that's what we're interested in.
    }
    R <- corrmat(maximized,design,alpha,1.0,conditioning)
    cholR <- chol(R)
    Rinv <- chol2inv(cholR)
    betahat <- sum(Rinv%*%response)/sum(Rinv)
    dev <- response-betahat
    sigma2hat <- drop(1/N*t(dev)%*%Rinv%*%dev)
    Einv <- (1/sigma2hat)*Rinv
    solres <- Einv%*%dev
    ret <- list(maximized=maximized,alpha=alpha,
                design=design,response=response,solres=solres,
                betahat=betahat, sigma2hat=sigma2hat,
                sigma2epsilonhat=0,
                conditioning=conditioning, Einv=Einv,
                bestlogl=bestlogl)
    class(ret) <- "gpmodel"
    return(ret)
}

computelogl.merr.pdim <- function(theta.tilde,design,response,alpha,N,conditioning)
{
    np <- length(theta.tilde)
    out <- .C("computeloglMerr",as.double(theta.tilde),as.double(design),
              as.double(response), as.double(alpha), as.double(conditioning),
              as.integer(N), as.integer(np), ans = double(1),
              PACKAGE="newlrcalib")
    return(out$ans)
}
maximizelikelihood.merr <- function(response,design,alpha=2,reps=100,
                                    conditioning=.Machine$double.eps,
                                    lrnuginit=1e-14,rho.ini=1)
{
    N <- length(response)
    design <- as.matrix(design)
    pdims <- ncol(design)+2
    rho.mat <- runif.sobol(reps,pdims-2)*rho.ini
    bestlogl <- -1e16
    bestrho <- NULL
    for(i in 1:reps)
    {
        rho.init <- rho.mat[i,]
        theta.tilde.init <- log(-4*log(rho.init)) # optim will optimize over theta.tilde instead of rho directly.
        theta.tilde.init <- c(theta.tilde.init,log(var(response)),log(lrnuginit))
        maximized <- optim(theta.tilde.init,computelogl.merr.pdim,method="Nelder-Mead",
                           control=list(alpha=1,beta=0.5,gamma=2,maxit=5000),design=design,
                           response=response,alpha=alpha,N=N,conditioning=conditioning)
        if(maximized$convergence!=0) warning("Optim did not converge!\n converge=",maximized$convergence)
        value <- -maximized$value
        if(value>bestlogl)
        {
            bestrho <- maximized$par
            bestlogl <- value
        }
    }
    if(is.null(bestrho)) stop("rho=NULL, intermediate step converged to rho=0!")
    maximized <- bestrho
    sigma2hat <- exp(maximized[pdims-1])
    sigma2epsilonhat <- exp(maximized[pdims])
    maximized <- maximized[1:(pdims-2)]
    maximized <- exp(-1/4*exp(maximized)) #convert back to rho-scale, since that's what we're interested in.
    nug <- sigma2epsilonhat/sigma2hat + conditioning
    E <- corrmat(maximized,design,alpha,1.0,nug)
    cholE <- chol(E)
    Einv <- (1/sigma2hat)*chol2inv(cholE)
    betahat <- sum(Einv%*%response)/sum(Einv)
    solres <- Einv%*%(response-betahat)
    ret <- list(maximized=maximized,alpha=alpha,
                design=design,response=response,solres=solres,
                betahat=betahat, sigma2hat=sigma2hat,
                sigma2epsilonhat=sigma2epsilonhat,
                conditioning=conditioning, Einv=Einv,
                bestlogl=bestlogl)
    class(ret) <- "gpmodel"
    return(ret)
}
rhodacepred.merr2<-function(model,newdesign)
{
    betahat <- model$betahat
    Rpred <- crosscorrmat(model$maximized,newdesign,model$design,model$alpha)
    resid <- model$response-betahat
    resid <- model$Einv%*%resid
    yhat <- betahat+model$sigma2hat*Rpred%*%resid
    return(yhat)
}
predict.gpmodel <- function(model,newdesign)
{
    betahat <- model$betahat
    s2hat <- model$sigma2hat
    se2hat <- model$sigma2epsilonhat
    Rpred <- crosscorrmat(model$maximized,newdesign,model$design,model$alpha)
    resid <- model$response-betahat
    resid <- model$Einv%*%resid
    yhat <- betahat+s2hat*Rpred%*%resid
    solER <- Rpred%*%model$Einv
    quadER <- apply(solER*Rpred,1,sum)
    bilER1 <- apply(solER,1,sum)
    quadE1 <- sum(model$Einv)
    part1 <- s2hat^2*quadER
    part2 <- (1-s2hat*bilER1)^2/quadE1
    mse <- s2hat+se2hat-part1+part2
    mse <- ifelse(mse>0,mse,0)
    ret <- list(yhat=yhat,mse=mse)
    return(ret)
}
## there is space to be optimized
updateGP <- function(model,newdes,newresp)
{
    design <- rbind(newdes,model$design)
    N <- nrow(design)
    response <- c(newresp,model$response)
    nug <- model$sigma2epsilonhat/model$sigma2hat + model$conditioning
    E <- corrmat(model$maximized,design,model$alpha,1.0,nug)
    cholE <- chol(E)
    Einv <- (1/model$sigma2hat)*chol2inv(cholE)
    betahat <- sum(Einv%*%response)/sum(Einv)
    model$design <- design
    model$response <- response
    model$betahat <- betahat
    model$Einv <- Einv
    model$solres <- Einv%*%(response-betahat)
    return(model)
}
maximizelikelihood.merr.lite <- function(response,design,pdims,alpha=2,reps=100,
                                         conditioning=.aMchine$double.eps,lrnuginit=1e-14,
                                         rho.ini=1)
{
    N <- length(response)
    rho.mat <- runif.sobol(reps,pdims)*rho.ini
    bestlogl <- -1e16
    bestrho <- NULL
    for(i in 1:reps)
    {
        rho.init <- rho.mat[i,]
        theta.tilde.init <- log(-4*log(rho.init)) # optim will optimize over theta.tilde instead of rho directly.
        theta.tilde.init <- c(theta.tilde.init,log(var(response)),log(lrnuginit))
        maximized <- optim(theta.tilde.init,computelogl.merr.pdim,method="Nelder-Mead",
                           control=list(alpha=1,beta=0.5,gamma=2,maxit=5000),design=design,
                           response=response,alpha=alpha,N=N,conditioning=conditioning)
        if(maximized$convergence!=0) warning("Optim did not converge!\n converge=",maximized$convergence)
        value <- -maximized$value
        if(value>bestlogl)
        {
            bestrho <- maximized$par
            bestlogl <- value
        }
    }
    if(is.null(bestrho)) stop("rho=NULL, intermediate step converged to rho=0!")
    return(bestlogl)
}
