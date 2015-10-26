#' Fit a MIMOSA model with baseline.
#'
#' @name MIMOSA2
#' @param Ntot \code{matrix} integer vector of total trials. One row per subject. Should have four columns named "ns1" "ns0" "nu1" and "nu0"
#' @param ns1 \code{numeric} integer vector of successes in condition 1 treatment s.
#' @param nu1 \code{numeric} integer vector of successes in condition 1 treatment u.
#' @param ns0 \code{numeric} integer vector of successes in condition 0 treatment s.
#' @param nu0 \code{numeric} integer vector of successes in condition 0 treatment u.
#' @param tol \code{numeric} tolerance for stopping criteria, change in relative log-likelihood.
#' @param inds \code{matrix} numeric indicator of assignments to initialize with.
#' @param maxit \code{numeric} maximum number of iterations
#' @usage MIMOSA2(Ntot,ns1,nu1,ns0,nu0)
#' @return \code{list} of fitted model parameters with components \code{z} \code{inds} \code{thetahat} \code{pi_est}
#' @export
#' @importFrom matrixStats logSumExp
#' @import data.table ggplot2
#' @seealso \link{ORTest} \link{ROC} \link{ROCPlot} \link{Boxplot} \link{simulate_MIMOSA2} \link{logit} \link{invlogit}
#' @examples
#' s = simulate_MIMOSA2()
#' R = MIMOSA2(Ntot=s$Ntot, ns1 = s$ns1, nu1 = s$nu1, nu0 = s$nu0, ns0 = s$ns0)
#'
MIMOSA2 = function(Ntot,ns1,nu1,ns0,nu0,tol=1e-10,inds=NULL,maxit=100){
  K=5
  #' Get the number of observations from the data.
  P = nrow(Ntot)

  #' Estimate proportions from data
  pu1_hat = prop.table(cbind(Ntot[, 1], nu1), 1)[, 2]
  ps1_hat = prop.table(cbind(Ntot[, 2], ns1), 1)[, 2]
  pu0_hat = prop.table(cbind(Ntot[, 3], nu0), 1)[, 2]
  ps0_hat = prop.table(cbind(Ntot[, 4], ns0), 1)[, 2]

  #' Indices of potential responders
  flag_ind = (ps1_hat-pu1_hat)>(ps0_hat-pu0_hat)
  
  #'Initialize parameter estimates
  inits = initialize(P,inds)
  thetahat = inits$thetahat
  pi_est = inits$pi_est
  inds = inits$ind

  #' Per observation complete data log likelihood matrix
  mat = cll(thetahat,
            Ntot,
            ns1,
            nu1,
            ns0,
            nu0)

  #'Current complete data log-likelihood
  llold = sum(t(t(mat) + log(pi_est)) * inds)

  #' Difference
  ldiff = Inf

  #' Maximum itertions 100 (will be configurable).Current iteration 0
  maxiter=maxit
  iter=0

  #'Fitting loop, alternate E and M steps,
  #' stop when relative change in ll is 1e-5
  brk=FALSE
  while (ldiff > tol) {
    iter=iter+1
    if(iter>maxiter){
      break;
    }
    #' optimize fnscale -1 for maximization.
    est = try(optim(
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
    ),silent = TRUE)
    if(!inherits(est,"try-error")){
      if(est$value<0){
        thetahat = est$par
        ldiff = abs(est$value-llold)/abs(est$value)
        llold = est$value
        #' New parameter estimates from optim
        #' Relative change in log-likelihood
        #' Update current log-likelihood and print it.
        if(maxiter>10)
          message("-Log-Likelihood: ", -llold)
      }else{
        brk=TRUE
      }
    }else{
      brk=TRUE
    }

    #' Calculate matrix of complete-data log likelihood for each observation.
    #' ns1, ks1, ns0, ks0, nu, ku, ns, ks, n1, k1,n0,k0,n,k
    mat = cll(thetahat,
              Ntot,
              ns1,
              nu1,
              ns0,
              nu0)

    #' Don't forget the mixing proportions.
    #' Should zero out the likelihood for !ind_flag to be responders
    mat = t(t(mat) + log(pi_est))
    mat[!flag_ind,c(1,5)]=-.Machine$integer.max

    #' Update the z's
    z = exp(mat - apply(mat, 1, function(x)
      matrixStats::logSumExp(x)))

    #' Components 2:K are non-responders, component 1 and 5 are  responders.
    #' Assign hierarchically to either the
    #' responder or non-responder components
    inds_resp = max.col(cbind(rowSums(z[, c(1,5)]), 1-rowSums(z[, c(1,5)])))

    inds = matrix(0, nrow = P, ncol = K)
    for (i in 1:nrow(z)) {
      v = rep(0, K)
      if ((inds_resp[i] == 1) #&
          #((ps1_hat > pu1_hat) & (ps1_hat > ps0_hat))[i]
          ) {
        v[c(1,5)][which.max(z[i, c(1,5)])] = 1
      } else{
        v[-c(1,5)][which.max(z[i, -c(1,5)])] = 1
      }
      inds[i, ] = v
    }

    #' update mixing proportions
    pi_est = colMeans(z)
    if(brk){
      break
    }
  }
  if(maxiter>10)
    cat("done\n")
  colnames(inds) = c("R","NR","NV1","NV2","R")
  l = list(z, inds, pi_est, thetahat,ps1_hat,ps0_hat,pu1_hat,pu0_hat,Ntot,ns1,nu1,ns0,nu0,-llold)
  names(l) = c("z","inds","pi_est","thetahat","ps1_hat","ps0_hat","pu1_hat","pu0_hat","Ntot","ns1","nu1","ns0","nu0","ll")
  return(l)
}

#' Fit a MIMOSA model with baseline with random restarts.
#'
#' @name MIMOSA2fit
#' @param Ntot \code{matrix} integer vector of total trials. One row per subject. Should have four columns named "ns1" "ns0" "nu1" and "nu0"
#' @param ns1 \code{numeric} integer vector of successes in condition 1 treatment s.
#' @param nu1 \code{numeric} integer vector of successes in condition 1 treatment u.
#' @param ns0 \code{numeric} integer vector of successes in condition 0 treatment s.
#' @param nu0 \code{numeric} integer vector of successes in condition 0 treatment u.
#' @param tol \code{numeric} tolerance for stopping criteria, change in relative log-likelihood.
#' @usage MIMOSA2(Ntot,ns1,nu1,ns0,nu0)
#' @return \code{list} of fitted model parameters with components \code{z} \code{inds} \code{thetahat} \code{pi_est}
#' @export
#' @importFrom matrixStats logSumExp
#' @import data.table ggplot2
#' @seealso \link{ORTest} \link{ROC} \link{ROCPlot} \link{Boxplot} \link{simulate_MIMOSA2} \link{logit} \link{invlogit}
#' @examples
#' s = simulate_MIMOSA2()
#' R = MIMOSA2(Ntot=s$Ntot, ns1 = s$ns1, nu1 = s$nu1, nu0 = s$nu0, ns0 = s$ns0)
#'
MIMOSA2fit = function(Ntot,ns1,nu1,ns0,nu0,tol=1e-10,restarts=10){
  best = Inf
  keep = NULL
  message("Performing ",restarts," random restarts.")
  for(i in 1:restarts){
   fit = MIMOSA2(Ntot,ns1,nu1,ns0,nu0,tol=1e-10,maxit=5)
   if(fit$ll<best){
     keep=fit
   }
  }
  keep = MIMOSA2(Ntot,ns1,nu1,ns0,nu0,tol=1e-10,inds=keep$inds)
  return(keep)
}
