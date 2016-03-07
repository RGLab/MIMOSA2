#' Fit a MIMOSA model with baseline.
#'
#' @name MIMOSA2
#' @param Ntot \code{matrix} integer vector of total trials. One row per subject. Should have four columns named "ns1" "ns0" "nu1" and "nu0"
#' @param ns1 \code{numeric} integer vector of successes in condition 1 treatment s.
#' @param nu1 \code{numeric} integer vector of successes in condition 1 treatment u.
#' @param ns0 \code{numeric} integer vector of successes in condition 0 treatment s.
#' @param nu0 \code{numeric} integer vector of successes in condition 0 treatment u.
#' @param tol \code{numeric} tolerance for stopping criteria, change in relative log-likelihood.
#' @param maxit \code{numeric} maximum number of iterations
#' @param verbose \code{logical} print the absolute difference of the sum of successive estimates of parameters. Defaults to FALSE
#' @usage MIMOSA2(Ntot,ns1,nu1,ns0,nu0,tol=1e-8,maxit=100,verbose=FALSE)
#' @return \code{list} of fitted model parameters with components \code{z} \code{inds} \code{thetahat} \code{pi_est}
#' @export
#' @importFrom matrixStats logSumExp
#' @import data.table ggplot2 optimx
#' @seealso \link{ORTest} \link{ROC} \link{ROCPlot} \link{Boxplot} \link{simulate_MIMOSA2} \link{logit} \link{invlogit}
#' @examples
#' s = simulate_MIMOSA2();
#' R = MIMOSA2(Ntot=s$Ntot, ns1 = s$ns1, nu1 = s$nu1, nu0 = s$nu0, ns0 = s$ns0,maxit=10)
MIMOSA2 = function(Ntot,ns1,nu1,ns0,nu0,tol=1e-8,maxit=100,verbose=FALSE){
  K=11
  rcomps = c(1:4)
  # Get the number of observations from the data.
  P = nrow(Ntot)

  # Estimate proportions from data
  pu1_hat = prop.table(cbind(Ntot[, "nu1"], nu1), 1)[, 2]
  ps1_hat = prop.table(cbind(Ntot[, "ns1"], ns1), 1)[, 2]
  pu0_hat = prop.table(cbind(Ntot[, "nu0"], nu0), 1)[, 2]
  ps0_hat = prop.table(cbind(Ntot[, "ns0"], ns0), 1)[, 2]

  # Initialize parameter estimates
  inits = initialize(P,Ntot=Ntot,ns1=ns1,nu1=nu1,ns0=ns0,nu0=nu0,K=K)
  thetahat = inits$thetahat
  thetahat[c(1,3,5,7)]=sort(thetahat[c(1,3,5,7)],decreasing=TRUE)
  pi_est = inits$pi_est
  z=inits$inds

  ps1=ns1/Ntot[,"ns1"]
  pu1=nu1/Ntot[,"nu1"]
  ps0=ns0/Ntot[,"ns0"]
  pu0=nu0/Ntot[,"nu0"]
  dp = ps1-pu1-ps0+pu0
  dp1 = ps1-pu1
  dp0 = ps0-pu0
  dpu = pu0-pu1

  # Per observation complete data log likelihood matrix
  mat = cll(par=thetahat,
            Ntot=Ntot,
            ns1=ns1,
            nu1=nu1,
            ns0=ns0,
            nu0=nu0)
  if(length(which((dp<0&dp1<0)|dp<0))>0)
    mat[(dp<0&dp1<0)|dp<0,]=t(apply(mat[(dp<0&dp1<0)|dp<0,,drop=FALSE],1,function(x)c(rep(-10000,4),x[5:11])))      # for(i in which (dpu<0&dp<0)){
  # for(i in which (dpu<0&dp<0)){
  #   mat[i,3]=min(mat[i,])
  # }
  mat = t(t(mat)+log1p(pi_est))
  mx = apply(mat,1,max)
  tmp = exp(mat-mx)
  z = (tmp/rowSums(tmp))
  pi_est = colMeans(z)
  #Current complete data log-likelihood

  ldiff=Inf
  # Maximum itertions 100 (will be configurable).Current iteration 0
  maxiter=maxit
  iter=0

  #Fitting loop, alternate E and M steps,
  # stop when relative change in ll is 1e-5
  brk=FALSE
  while (all(ldiff > tol)) {
    iter=iter+1
    if(iter>maxiter){
      break;
    }
    est = try(optimx(
      par = thetahat,
      fn = sumcll,
      pi_est = pi_est,
      z = z,
      Ntot = Ntot,
      ns1 = ns1,
      nu1 = nu1,
      ns0 = ns0,
      nu0 = nu0,
      method=c("newuoa","bobyqa")),silent = TRUE)
    if(!inherits(est,"try-error")){
      est=est[order(est[,"convcode"],est[,"value"],decreasing=FALSE)[1],,drop=FALSE]
      est = unlist(est[1:8])
      est[c(1,3,5,7)]=sort(est[c(1,3,5,7)],decreasing=TRUE)
      mat_new = cll(par=unlist(est[1:8]),
                Ntot=Ntot,
                ns1=ns1,
                nu1=nu1,
                ns0=ns0,
                nu0=nu0)
      if(length(which((dp<0&dp1<0)|dp<0))>0)
        mat_new[(dp<0&dp1<0)|dp<0,]=t(apply(mat_new[(dp<0&dp1<0)|dp<0,,drop=FALSE],1,function(x)c(rep(-10000,4),x[5:11])))      # for(i in which (dpu<0&dp<0)){
      #   mat[i,3]=min(mat[i,])
      # }
      mat_new=t(t(mat_new)+log1p(pi_est))
      mx = apply(mat_new,1,max)
      tmp = exp(mat_new-mx)
      z_new = (tmp/rowSums(tmp))
      pi_new = colMeans(z_new)
      if(verbose)
        cat(sum(abs(unlist(est[c(1,2,3,4,5,6,7,8)])-thetahat[c(1,2,3,4,5,6,7,8)])),"\n")
      ldiff = abs(c(unlist(est[c(1,2,3,4,5,6,7,8)]),pi_new)-c(thetahat[c(1,2,3,4,5,6,7,8)],pi_est))
      thetahat = unlist(est[1:8])
      z=z_new
      pi_est = colMeans(z)
    }else{
      brk=TRUE
    }



    # Assign hierarchically to either the
    # responder or non-responder components
    inds_resp = max.col(cbind(rowSums(z[, rcomps,drop=FALSE]), 1-rowSums(z[, rcomps,drop=FALSE])))

    inds = matrix(0, nrow = P, ncol = K)
    for (i in 1:nrow(z)) {
      v = rep(0, K)
      if ((inds_resp[i] == 1)) {
        v[rcomps][which.max(z[i, rcomps])] = 1
      } else{
        v[-rcomps][which.max(z[i, -rcomps])] = 1
      }
      inds[i, ] = v
    }

    # update mixing proportions
    if(brk){
      break
    }

  }
  #if(maxiter>10)
   # cat("done\n")
  colnames(inds) = 1:K
  l = list(z, inds, pi_est, thetahat,ps1_hat,ps0_hat,pu1_hat,pu0_hat,Ntot,ns1,nu1,ns0,nu0)
  names(l) = c("z","inds","pi_est","thetahat","ps1_hat","ps0_hat","pu1_hat","pu0_hat","Ntot","ns1","nu1","ns0","nu0")
  return(l)
}
