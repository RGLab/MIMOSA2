#' Simulate data for MIMOSA2 model, ICS with baseline.
#'
#' @param effect \code{numeric} effect size for ps1 - pu1 - ps0 - pu0
#' @param phi \code{numeric} integer, the precision.
#' @param P \code{numeric} number of subjects.
#' @param seed \code{numeric} random seed
#'
#' @return \code{list} with components "Ntot" "ns0" "nu0" "ns1" "nu1" "truth"
#' @export
#' @seealso \link{MIMOSA2}
#' @examples
#' s = simulate_MIMOSA2()
simulate_MIMOSA2 = function(effect = 5e-4,
                    phi = 10000,
                    P = 100) {

  #' Run 100 simulations, dropping any errors (i.e. simulations that error out)
  K = 5

  #'proportion of vaccine specific responses
  pis = prop.table(rgamma(K, 1))
  pis = sort(pis, decreasing = TRUE)
  R = NULL
  D = 4

  #' Total observations simulated
  #' from uniform with between 50,000 and 100,000 cells
  Ntot = matrix(round(runif(P * D, 50000, 100000)), ncol = D, nrow = P)

  #' Hyperprior mean for stimulated time 0
  MS0 = 2.5e-4

  #' Mean of hyperprior (stimulated time 1) $\alpha/(\alphaa+\beta)$
  MS1 = MS0 + effect

  #' Non Stimulated
  MU = 1e-4

  #'Precision of hyperprior
  #'$\alpha+\beta = \phi$
  PHI = phi

  #' When there is a response ps1 > ps0 & ps1 > pu
  #' When there is no response ps0 = ps1 = pu
  #' When there is no vaccine-specific
  #' response ps0 = ps1, and one or both of
  #' ps1>pu,   ps0 > pu may also be true

  #' number of responders / non-responders / non-specific responders
  nresp = round(P * pis[1])
  n_non_resp = round(pis[2] * P)
  n_novaccine1 = round(P * pis[3])
  n_novaccine2 = round(P * pis[4])
  nresp2 = round(P - (nresp+n_non_resp+n_novaccine1+n_novaccine2))

  if((nresp2+nresp+n_non_resp+n_novaccine1+n_novaccine2)!=P){
    stop("try again")
  }
  #'Simulate responders
  ps0 = rep(0, nresp)
  ps1 = rep(0, nresp)
  pu =  rbeta(nresp, MU * PHI, (1 - MU) * PHI)
  while (any(ps1 <= ps0 | ps1 < pu)) {
    bar = ps1 <= ps0 | ps1 < pu
    foo = sum(bar)
    ps0[bar] = rbeta(foo, MS0 * PHI, (1 - MS0) * PHI)
    ps1[bar] = rbeta(foo, MS1 * PHI, (1 - MS1) * PHI)
  }
  if (effect == 0) {
    ps1 = ps0
  }

  #'Simulate non-responders. ps0=ps1=pu
  foo = rbeta(n_non_resp, MU * PHI, (1 - MU) * PHI)
  ps0 = c(ps0, foo)
  ps1 = c(ps1, foo)
  pu =  c(pu, foo)

  #' Simulate from component 3. S0 is a response but s1 is not
  foo = rbeta(n_novaccine1, MU * PHI, (1 - MU) * PHI)
  bar = rbeta(n_novaccine1, MS0 * PHI, (1 - MS0) * PHI)
  while (any(bar <= foo)) {
    w = bar <= foo
    l =  sum(w)
    bar[w] = rbeta(l, MS0 * PHI, (1 - MS0) * PHI)
  }
  ps0 = c(ps0, bar)
  ps1 = c(ps1, foo)
  pu =  c(pu, foo)

  #' Simulate non-specific responses
  #' where basline = post-vaccine simulated
  #' from hyper-prior with $\mu_{s0}$.
  foo = rbeta(n_novaccine2, MS0 * PHI, (1 - MS0) * PHI)
  bar = rbeta(n_novaccine2, MU * PHI, (1 - MU) * PHI)
  while (any(foo <= bar)) {
    w = foo <= bar
    l =  sum(w)
    foo[w] = rbeta(l, MS1 * PHI, (1 - MS1) * PHI)
  }
  ps0 = c(ps0, foo)
  ps1 = c(ps1, foo)
  pu = c(pu, bar)

  #'simulate other responder group
  baz =  rbeta(nresp2, MU * PHI, (1 - MU) * PHI)
  foo = rbeta(nresp2, MS1 * PHI, (1 - MS1) * PHI)
  bar = baz
  while (any(foo <= baz)) {
    w = foo <= baz
    l = sum(w)
    foo[w] = rbeta(l, MS1 * PHI, (1 - MS1) * PHI)
  }
  if (effect == 0) {
    foo = bar
  }
  pu = c(pu, baz)
  ps0 = c(ps0, bar)
  ps1 = c(ps1, foo)

  #' Simulate count observations from binomial.
  colnames(Ntot) = c("nu1", "ns1", "nu0", "ns0")
  nu1 = rbinom(P, Ntot[, 1], pu)
  ns1 = rbinom(P, Ntot[, 2], ps1)
  nu0 = rbinom(P, Ntot[, 3], pu)
  ns0 = rbinom(P, Ntot[, 4], ps0)


  #' Empirical estimates of proportions
  pu1_hat = prop.table(cbind(Ntot[, 1], nu1), 1)[, 2]
  ps1_hat = prop.table(cbind(Ntot[, 2], ns1), 1)[, 2]
  pu0_hat = prop.table(cbind(Ntot[, 3], nu0), 1)[, 2]
  ps0_hat = prop.table(cbind(Ntot[, 4], ns0), 1)[, 2]

  #' True response categories
  truth = rep(
    c("R", "NR", "NV1", "NV2", "R"),
    c(nresp, n_non_resp, n_novaccine1, n_novaccine2, nresp2)
  )
  l=list("Ntot"=Ntot, "ns0"=ns0, "ns1"=ns1, "nu0"=nu0, "nu1"=nu1, "truth"=truth)
  names(l) = c("Ntot","ns0","ns1","nu0","nu1","truth")
  return(l)
}
