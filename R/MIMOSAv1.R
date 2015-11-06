MIMOSAv1 = function(Ntot,ns,nu,tol=1e-10,inds=NULL,maxit=100){
  K=2
  #' Get the number of observations from the data.
  P = nrow(Ntot)

  #' Estimate proportions from data
  pu_hat = prop.table(cbind(Ntot[, "nu"], nu), 1)[, 2]
  ps_hat = prop.table(cbind(Ntot[, "ns"], ns), 1)[, 2]

  #' Indices of potential responders
  flag_ind = ps_hat>pu_hat

  #'Initialize parameter estimates
  inits = initialize(P,inds,K=2)
  thetahat = inits$thetahat
  pi_est = inits$pi_est
  inds = inits$ind

  #' Per observation complete data log likelihood matrix
  mat = cllm0m1(par=thetahat,
            Ntot=Ntot,
            ns=ns,
            nu=nu)
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
      fn = sumcllm0m1,
      pi_est = pi_est,
      inds = inds,
      Ntot = Ntot,
      ns = ns,
      nu = nu,
      control = list(fnscale = -1, maxit = 100000)
    ),silent = TRUE)
    if(!inherits(est,"try-error")){
      if(est$value<0&est$value>llold){
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
    mat = cllm0m1(par=thetahat,
                  Ntot=Ntot,
                  ns=ns,
                  nu=nu)

    #' Don't forget the mixing proportions.
    #' Should zero out the likelihood for !ind_flag to be responders
    mat[!flag_ind,c(2)]=-.Machine$integer.max
    mat = t(t(mat) + log(pi_est))


    #' Update the z's
    z = exp(mat - apply(mat, 1, function(x)
      matrixStats::logSumExp(x)))

    #' Components 2:K are non-responders, component 1 and 5 are  responders.
    #' Assign hierarchically to either the
    #' responder or non-responder components
    inds_resp = max.col(cbind(rowSums(z[, c(2),drop=FALSE]), 1-rowSums(z[, c(2),drop=FALSE])))

    inds = matrix(0, nrow = P, ncol = K)
    for (i in 1:nrow(z)) {
      v = rep(0, K)
      if ((inds_resp[i] == 1) #&
          #((ps1_hat > pu1_hat) & (ps1_hat > ps0_hat))[i]
      ) {
        v[c(2)][which.max(z[i, c(2)])] = 1
      } else{
        v[-c(2)][which.max(z[i, -c(2)])] = 1
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
  colnames(inds) = c("NR","NR")
  l = list(z, inds, pi_est, thetahat,ps_hat,pu_hat,Ntot,ns,nu,-llold)
  names(l) = c("z","inds","pi_est","thetahat","ps_hat","pu_hat","Ntot","ns","nu","ll")
  return(l)
}
