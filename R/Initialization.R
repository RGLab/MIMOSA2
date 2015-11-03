#' Initialize model parameters
#' @details The \code{thetahat} are the initial hyperparameter estimates. The \code{pi_est} are the initial mixing proportions. The \code{inds} are the initial assignments of data to model components.
#' @param P \code{numeric} The number of observations / subjects.
#' @param Nott
#' @param ns1
#' @param nu1
#' @param ns0
#' @param nu0
#' @return \code{list} of initialized parameter estimates with components \code{inds} \code{theta_hat} \code{pi_est}.
initialize = function(P,Ntot,ns1,nu1,ns0,nu0,random=FALSE,K) {
  reps=ifelse(ncol(Ntot)>2,4,2)
  K = ifelse(ncol(Ntot)>2,K,2)
  thetahat = rep(c(logit(0.5), log(100)), reps)
    #'Initialize random hard assignments
  if(random){
    inds = t(sapply(sample(1:K, P, replace = TRUE), function(x) {
     y = rep(0, K)
      y[x] = 1
      y
    }))
  }else{
    ortest = ORTest(Ntot,ns1,nu1,ns0,nu0)
    o=order(ortest,decreasing=FALSE)
    inds = matrix(0,ncol=K,nrow=nrow(Ntot))
    for(i in seq_along(1:round(0.2*length(ortest)))){
      inds[o[i],sample(1:4,1)]=1
    }
    for(i in (round(0.2*length(ortest))+1):length(ortest)){
      inds[o[i],sample(c(5:8,4))]=1
    }
  }

  pi_est = colMeans(inds)

  #initialize thetahat from inds and data.
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

  #'Estimate pi, mixing proportions
  return(list(thetahat=thetahat,pi_est=pi_est,inds=inds))
}
