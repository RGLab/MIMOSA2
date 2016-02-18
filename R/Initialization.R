#' Initialize model parameters
#' @details The \code{thetahat} are the initial hyperparameter estimates. The \code{pi_est} are the initial mixing proportions. The \code{inds} are the initial assignments of data to model components.
#' @param P \code{numeric} The number of observations / subjects.
#' @param Nott
#' @param ns1
#' @param nu1
#' @param ns0
#' @param nu0
#' @return \code{list} of initialized parameter estimates with components \code{inds} \code{theta_hat} \code{pi_est}.
initialize = function(P,Ntot,ns1,nu1,ns0,nu0,K) {
  reps=ifelse(ncol(Ntot)>2,4,2)
  K = ifelse(ncol(Ntot)>2,K,2)
  thetahat = rep(c(logit(0.5), log(10000)), reps)
  #empirical estimates:
  deltap  = (ns1/Ntot[,"ns1"]-nu1/Ntot[,"nu1"]-ns0/Ntot[,"ns0"]+nu0/Ntot[,"nu0"])>0
  r1 = ns1/Ntot[,"ns1"]>nu1/Ntot[,"nu1"]
  r0 = ns0/Ntot[,"ns0"]>nu0/Ntot[,"nu0"]
  indicator = deltap>0&(r1&!r0)
    inds = matrix(0,ncol=K,nrow=nrow(Ntot))
    if(sum(indicator)>0){
    		for(i in which(indicator)){
    			inds[i,1]=1
    		}
    }
    for(i in which(!indicator)){
    	inds[i,5]=1
    }
   inds2 = rowSums(inds[,1:4] )>0
   thetahat[1] = logit(mean((ns1/Ntot[,"ns1"])[inds2]))
   thetahat[5] = logit(mean((nu1/Ntot[,"nu1"])[inds2]))
   thetahat[3] = logit(mean((ns0/Ntot[,"ns0"])[inds2]))
   thetahat[7] = logit(mean((nu0/Ntot[,"nu0"])[inds2]))
   thetahat[2] = 2
   thetahat[6] = 2

  pi_est = colMeans(inds)

  #initialize thetahat from inds and data.
  # make several initialization runs and take the best sometimes it fails miserably.
  bestest=NULL
  for(k in 1:10){
    est = optim(
      par = thetahat,
      fn = sumcll,
      pi_est = pi_est,
      inds = inds,
      Ntot = Ntot,
      ns1 = ns1,
      nu1 = nu1,
      ns0 = ns0,
      nu0 = nu0,
      control = list(fnscale = -1, maxit = 100000)
    )
    if(k==1&est$convergence==0)
      bestest=est
    if(est$value>bestest$value)
      bestest=est
  }

  thetahat = bestest$par

  #'Estimate pi, mixing proportions
  return(list(thetahat=thetahat,pi_est=pi_est,inds=inds))
}
