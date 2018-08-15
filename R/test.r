computelogl.1D.R <- function(rho,design,response,alpha,N,conditioning)
{
    R <- corrmat(rho,design,alpha,1.0,conditioning)
    cholR <- chol(R)
    Rinv <- chol2inv(cholR)
    logdetR=2*sum(log(diag(cholR)))
    betahat <- sum(Rinv%*%response)/sum(Rinv)
    dev <- response - betahat
    sigma2hat <- 1/N*t(dev)%*%Rinv%*%dev
    logl <- -N*log(sigma2hat)-logdetR
    -logl
}
computelogl.pdim.R <- function(theta.tilde,design,response,alpha,N,conditioning)
{
    rho <- exp(-1/4*exp(theta.tilde))
    R <- corrmat(rho,design,alpha,1.0,conditioning)
    cholR <- chol(R)
    Rinv <- chol2inv(cholR)
    logdetR <- 2*sum(log(diag(cholR)))
    betahat <- sum(Rinv%*%response)/sum(Rinv)
    dev <- response - betahat
    sigma2hat <- 1/N*t(dev)%*%Rinv%*%dev
    logl <- -N*log(sigma2hat)-logdetR
    -logl
}

maximizelikelihood.R <- function(response,design,alpha=2,reps=100,
                               conditioning=.Machine$double.eps,rho.ini=1,
                               bound1d=FALSE)
{
    design <- as.matrix(design)
    N <- length(response)
    pdims <- ncol(design)
    if(pdims==1) #1D CASE
    {
        upbound <- if(bound1d) exp(-0.25/darg(NULL,design)$max) else (1-.Machine$double.eps)
        maximized <- optimize(computelogl.1D.R,interval=c(.Machine$double.eps,upbound),
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
            maximized <- optim(theta.tilde.init,computelogl.pdim.R,method="Nelder-Mead",
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
computelogl.merr.pdim.R <- function(theta.tilde,design,response,alpha,N,conditioning)
{
    np=length(theta.tilde)
    sigma2hat=exp(theta.tilde[np-1])       #optim optimizes sigma2, sigma2epsilon on log scale which is -inf..inf, so
    sigma2epsilonhat=exp(theta.tilde[np])  #here we conver it back to 0..inf scale by exponentiating
    theta.tilde=theta.tilde[1:(np-2)]
    rho=exp(-1/4*exp(theta.tilde)) #optim is optimizes the log(theta) on -inf, inf, so here convert back to true rho
    R <- corrmat(rho,design,alpha,sigma2hat,sigma2epsilonhat+conditioning)
    cholR=chol(R)
    Rinv=chol2inv(cholR)
    logdetR=2*sum(log(diag(cholR)))
    betahat=sum(Rinv%*%response)/sum(Rinv)
    dev=response-betahat
    logl=-1/2*t(dev)%*%Rinv%*%dev-1/2*logdetR
    -logl #return negative logl since optim minimizes function rather than maximizes
}
