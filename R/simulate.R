#' Simulate data for MIMOSA2 model, ICS with baseline.
#'
#' @param effect \code{numeric} effect size for ps1 - pu1 - ps0 + pu0. Default 5e-4
#' @param bg_effect \code{numeric} effect size for background pu0-pu1. Default 0
#' @param baseline_stim_effect \code{numeric} baseline stimulation effect. Default 2.5e-4
#' @param baseline_background \code{numeric} baseline background effect. Default 1e-4
#' @param phi \code{numeric} integer, the precision. Default 5000
#' @param P \code{numeric} number of subjects.'
#' @param rng \code{numeric} vector of length 2. The range of total cell counts to simulate from.
#' @return \code{list} with components "Ntot" "ns0" "nu0" "ns1" "nu1" "truth"
#' @export
#' @seealso \link{MIMOSA2}
#' @examples
#' s = simulate_MIMOSA2()
simulate_MIMOSA2 = function(effect = 5e-4, bg_effect = 0,baseline_stim_effect=2.5e-4,baseline_background=1e-4,
                    phi = 5000,
                    P = 100,rng = c(100000,150000)) {

  #' Run 100 simulations, dropping any errors (i.e. simulations that error out)
  K = 8
  n=rep(0,K)
  while(sum(n)!=P|any(P==0)){
  #'proportion of vaccine specific responses
  pis = (prop.table(runif(K)))
  #pis = sort(pis, decreasing = TRUE)
  R = NULL
  D = 4

  #' Total observations simulated
  #' from uniform with between 100,000 and 150,000 cells
  Ntot = matrix(round(runif(P * D, rng[1], rng[2])), ncol = D, nrow = P)

  #' Hyperprior mean for stimulated time 0
  MU0 = baseline_background

  MS0 = baseline_stim_effect+MU0

  #' Non Stimulated post vaccine
  MU1 = MU0 + bg_effect

  #' Mean of hyperprior (stimulated time 1) $\alpha/(\alpha+\beta)$
  MS1 = MU1 + effect


  #'Precision of hyperprior
  #'$\alpha+\beta = \phi$
  PHI = rep(phi,4)

  #' There are 8 model components.
  #' 1. all different
  #' 2. s0=u0
  #' 3. s1 = s0
  #' 4. u1 = u0
  #' 5. s0 = u0, s1 = u1
  #' 6. s1 = u1
  #' 7. s1 = u1 = s0 = u0
  #' 8. s1 = s0, u1 = u0

  #' number of responders / non-responders / non-specific responders
  n = round(P * pis)
  }

  PS0=PS1=PU0=PU1=NULL
  #'Simulate from component 1
  k=1
  ps0 = rep(0, n[k])
  ps1 = rep(0, n[k])
  pu1 = rbeta(n[k], MU1 * PHI[1], (1 - MU1) * PHI[1])
  pu0 = rbeta(n[k], MU0 * PHI[2], (1-MU0) * PHI[2])
  while (any(ps1 <= pu1 | ps0 <= pu0 | ps1-pu1 <= ps0 - pu0)) {
    bar = ps1 <= pu1 | ps0 <= pu0 | ps1-pu1 <= ps0 - pu0
    foo = sum(bar)
    ps0[bar] = rbeta(foo, MS0 * PHI[3], (1 - MS0) * PHI[4])
    ps1[bar] = rbeta(foo, MS1 * PHI[4], (1 - MS1) * PHI[4])
  }
  if (effect == 0) {
    ps1 = ps0
  }
  PU1=c(PU1,pu1)
  PU0=c(PU0,pu0)
  PS1=c(PS1,ps1)
  PS0=c(PS0,ps0)

  #'Simulate from component 2
  k=2
  pu1 = rbeta(n[k], MU1 * PHI[1], (1 - MU1) * PHI[1])
  pu0 = rbeta(n[k], MU0 * PHI[2], (1 - MU0) * PHI[2])
  ps0 = pu0
  ps1 = rbeta(n[k],MS1*PHI[4], (1-MS1)*PHI[4])

  while(any(ps1-pu1 <= ps0 - pu0|ps1 <= pu1)){
    bar = ps1-pu1 <= ps0 - pu0|ps1 <= pu1
    foo = sum(bar)
    ps1[bar] = rbeta(foo, MS1 * PHI[4], (1 - MS1) * PHI[4])
    pu1[bar] = rbeta(foo, MU1 * PHI[1], (1 - MU1) * PHI[1])
  }

  PU1=c(PU1,pu1)
  PU0=c(PU0,pu0)
  PS1=c(PS1,ps1)
  PS0=c(PS0,ps0)


  #' Simulate from component 3.
  k=3
  ps0 = ps1 = rbeta(n[k],MS1*PHI[4],(1-MS1)*PHI[4])
  pu0 = rbeta(n[k],MU0*PHI[2],(1-MU0)*PHI[2])
  pu1 = rbeta(n[k],MU1*PHI[1],(1-MU1)*PHI[1])

  while(any(ps1-pu1 <= ps0 - pu0|ps1<=pu1| pu0<=pu1)){
    bar = ps1-pu1 <= ps0 - pu0|ps1<=pu1| pu0<=pu1
    foo = sum(bar)
    pu0[bar] = rbeta(foo, MU0 * PHI[2], (1 - MU0) * PHI[2])
    pu1[bar] = rbeta(foo, MU1 * PHI[1], (1 - MU1) * PHI[1])
  }
  PU1=c(PU1,pu1)
  PU0=c(PU0,pu0)
  PS1=c(PS1,ps1)
  PS0=c(PS0,ps0)

  #' Simulate from component 4.
  k=4
  pu0 = pu1 = rbeta(n[k],MU0*PHI[2],(1-MU0)*PHI[2])
  ps1 = rbeta(n[k],MS1*PHI[4],(1-MS1)*PHI[4])
  ps0 = rbeta(n[k],MS0*PHI[3],(1-MS0)*PHI[3])

  iter = 0
  while(any(ps1-pu1 <= ps0 - pu0|ps1<=pu1| ps1<=ps0)){
    bar = ps1-pu1 <= ps0 - pu0|ps1<=pu1| ps1<=ps0
    foo = sum(bar)
    ps0[bar] = rbeta(foo, MS0 * PHI[3], (1 - MS0) * PHI[3])
    ps1[bar] = rbeta(foo, MS1 * PHI[4], (1 - MS1) * PHI[4])
  }
  PU1=c(PU1,pu1)
  PU0=c(PU0,pu0)
  PS1=c(PS1,ps1)
  PS0=c(PS0,ps0)

  #'simulate from component 5
  k=5
  ps1 = pu1 = rbeta(n[k],MU0*PHI[2],(1-MU0)*PHI[2])
  ps0 = pu0 = rbeta(n[k],MU1*PHI[1],(1-MU1)*PHI[1])
  PU1=c(PU1,pu1)
  PU0=c(PU0,pu0)
  PS1=c(PS1,ps1)
  PS0=c(PS0,ps0)

  #'component 6
  k=6
  ps1 = pu1 = rbeta(n[k],MU1*PHI[2],(1-MU1)*PHI[2])
  ps0  = rbeta(n[k],MS0*PHI[3],(1-MS0)*PHI[3])
  pu0  = rbeta(n[k],MU0*PHI[2],(1-MU0)*PHI[2])
  #s0 should be greater than u0
  while(any(ps0<pu0)){
    bar = ps0<pu0
    foo = sum(bar)
    ps0[bar] = rbeta(foo, MS0 * PHI[3], (1 - MS0) * PHI[3])
    pu0[bar] = rbeta(foo, MU0 * PHI[2], (1 - MU0) * PHI[2])
  }
  PU1=c(PU1,pu1)
  PU0=c(PU0,pu0)
  PS1=c(PS1,ps1)
  PS0=c(PS0,ps0)

  #'component 7
  k=7
  ps1=ps0=pu1=pu0 = rbeta(n[k],MU0*PHI[2],(1-MU0)*PHI[2])
  PU1=c(PU1,pu1)
  PU0=c(PU0,pu0)
  PS1=c(PS1,ps1)
  PS0=c(PS0,ps0)

  #'component 8
  k=8
  if(n[k]>0){
  ps0 = ps1 = rbeta(n[k],MS0*PHI[3],(1-MS0)*PHI[3])
  pu0 = pu1 = rbeta(n[k],MU0*PHI[2],(1-MU0)*PHI[2])
  PU1=c(PU1,pu1)
  PU0=c(PU0,pu0)
  PS1=c(PS1,ps1)
  PS0=c(PS0,ps0)
}

  #' Simulate count observations from binomial.
  colnames(Ntot) = c("nu1", "ns1", "nu0", "ns0")
  nu1 = rbinom(P, Ntot[, "nu1"], PU1)
  ns1 = rbinom(P, Ntot[, "ns1"], PS1)
  nu0 = rbinom(P, Ntot[, "nu0"], PU0)
  ns0 = rbinom(P, Ntot[, "ns0"], PS0)


  #' Empirical estimates of proportions
  pu1_hat = prop.table(cbind(Ntot[, "nu1"], nu1), 1)[, 2]
  ps1_hat = prop.table(cbind(Ntot[, "ns1"], ns1), 1)[, 2]
  pu0_hat = prop.table(cbind(Ntot[, "nu0"], nu0), 1)[, 2]
  ps0_hat = prop.table(cbind(Ntot[, "ns0"], ns0), 1)[, 2]

  #' True response categories
  truth = rep(
    c("R1","R2","R3","R4","NR1","NR2","NR3","NSR"),
    n
  )
  l=list("Ntot"=Ntot, "ns0"=ns0, "ns1"=ns1, "nu0"=nu0, "nu1"=nu1, "truth"=truth)
  names(l) = c("Ntot","ns0","ns1","nu0","nu1","truth")
  return(l)
}
