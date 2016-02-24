#' Initialize model parameters
#' @details The \code{thetahat} are the initial hyperparameter estimates. The \code{pi_est} are the initial mixing proportions. The \code{inds} are the initial assignments of data to model components.
#' @param P \code{numeric} The number of observations / subjects.
#' @param Ntot matrix of total counts
#' @param ns1 vector of counts at time 1 stimulated
#' @param nu1 vector of counts at time 1 unstimulated
#' @param ns0 vector of counts at time 0 stimulated
#' @param nu0 vector of counts at time 0 unstimulated
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
    			inds[i,]=1/K
    		}
    }
    for(i in which(!indicator)){
    	inds[i,]=1/K
    }
   inds2 = rowSums(inds[,1:4] )>0
   thetahat[1] = logit(mean((ns1/Ntot[,"ns1"])[indicator]))
   thetahat[5] = logit(mean((nu1/Ntot[,"nu1"])))
   thetahat[3] = logit(mean((ns0/Ntot[,"ns0"])[indicator]))
   thetahat[7] = logit(mean((nu0/Ntot[,"nu0"])))
   thetahat[2] = 2
   thetahat[4] = 2
   thetahat[8] = 2
   thetahat[6] = 2
   thetahat[c(1,3,5,7)][is.infinite(thetahat[c(1,3,5,7)])]=log(1e-6)

  pi_est = colMeans(inds)

  #initialize thetahat from inds and data.
  # make several initialization runs and take the best sometimes it fails miserably.
  bestest=NULL
  bestest$value=-Inf
  for(k in 1:5){
    est = try(optim(
      par = thetahat,
      fn = sumcll,
      pi_est = pi_est,
      z = inds,
      Ntot = Ntot,
      ns1 = ns1,
      nu1 = nu1,
      ns0 = ns0,
      nu0 = nu0
    ))
    if(inherits(est,"try-error")){browser()}
    if(k==1&est$convergence==0&est$value>bestest$value)
      bestest=est
    else if (est$convergence!=0){
      est = optimx(
        par = thetahat,
        fn = sumcll,
        method="bobyqa",
        pi_est = pi_est,
        z = inds,
        Ntot = Ntot,
        ns1 = ns1,
        nu1 = nu1,
        ns0 = ns0,
        nu0 = nu0)
      if(est$convcode!=0){
        stop("failed to initialize")
      }else if(est$value>bestest$value){
        bestest$par=unlist(est[1:8])
        bestest$value = est$value
      }
    }
  }

  thetahat = bestest$par

  #'Estimate pi, mixing proportions
  return(list(thetahat=thetahat,pi_est=pi_est,inds=inds))
}
