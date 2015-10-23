#' Initialize model parameters
#' @details The \code{thetahat} are the initial hyperparameter estimates. The \code{pi_est} are the initial mixing proportions. The \code{inds} are the initial assignments of data to model components.
#' @param P \code{numeric} The number of observations / subjects.
#'
#' @return \code{list} of initialized parameter estimates with components \code{inds} \code{theta_hat} \code{pi_est}.
initialize = function(P) {
  K=5
  #'Initialize random hard assignments
  inds = t(sapply(sample(1:K, P, replace = TRUE), function(x) {
    y = rep(0, K)
    y[x] = 1
    y
  }))
  
  #'Initialize parameter estimates.
  thetahat = rep(c(logit(0.5), log(100)), 3)
  
  #'Estimate pi, mixing proportions
  pi_est = colMeans(inds)
  return(list(thetahat=thetahat,pi_est=pi_est,inds=inds))
}
