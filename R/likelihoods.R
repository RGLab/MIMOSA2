.getpars = function(thetahat,Ntot,ns1,nu1,ns0,nu0){
  mu=.getpars2(thetahat)[c(1:4)]
  nu = .getpars2(thetahat)[c(5:6)]
  v = (mu*(1-mu))/(1+c(nu,nu))
  pb=rbind(mu,v)
  return(pb)
}


.getpars2 = function(thetahat){
  c(invlogit(thetahat[c(1,3,5,7)]),exp(thetahat[c(2,6)]))
}
.clamp <- function(x,lwr=1e-10,upr=9e-10) {
  pmin(upr,pmax(lwr,x))
}

#' Beta-binomial log likelihood
#'
#' @name bbll
#' @param par vector of beta-binomial parameters \code{par[1] = logit(a / (a+b)), par[2] = log(a+b)}.
#' @param n \code{numeric} integer, the total trials. Can be a vector.
#' @param k \code{numeric} integer, the number of successes. Can be a vector.
#'
#' @return vector of beta-binomial log likelihoods for data \code{n,k} with parameters \code{par}
#' @export
#'
#' @examples
#' bbll(c(logit(0.5),log(100)),100,50)
bbll = function(par, n, k) {
  alpha = invlogit(par[1]) * exp(par[2])
  beta = (1 - invlogit(par[1])) * exp(par[2])
   lgamma(k + alpha) + lgamma(n - k + beta) -
    lgamma(n + alpha + beta) +
    lgamma(alpha + beta) - lgamma(alpha) -
    lgamma(beta)
}

const = function(n,k){
  lgamma(n+1)-lgamma(k+1)-lgamma(n-k+1)
}

.ll1 = function(par, Ntot, ns1, nu1, ns0, nu0) {
#s1 > u1 and s0>u0
  (bbll(par[c(1,2)], Ntot[, "ns1"], ns1) +
    bbll(par[c(3,2)],  Ntot[, "ns0"], ns0)+
    bbll(par[c(5,6)],  Ntot[, "nu1"], nu1)+
    bbll(par[c(7,8)],   Ntot[,"nu0"],nu0))#+log(ns1/Ntot[, "ns1"]-nu1/Ntot[, "nu1"]-ns0/Ntot[, "ns0"]+nu0/Ntot[, "nu0"]+1)+log(ns1/Ntot[, "ns1"]-nu1/Ntot[, "nu1"]+1)+log(ns0/Ntot[, "ns0"]+nu0/Ntot[, "nu0"]+1)
}

.ll2 = function(par, Ntot, ns1, nu1, ns0, nu0) {
#s1 > u1
  bbll(par[c(1, 2)], Ntot[, "ns1"] , ns1)+
    bbll(par[c(5,6)], Ntot[, "nu1"],  nu1)+
  bbll(par[c(7,8)],Ntot[,"nu0"]+Ntot[,"ns0"],ns0+nu0)#+log(ns1/Ntot[, "ns1"]-nu1/Ntot[, "nu1"]-ns0/Ntot[, "ns0"]+nu0/Ntot[, "nu0"]+1)+log(ns1/Ntot[, "ns1"]-nu1/Ntot[, "nu1"]+1)#+
}

.ll3 = function(par, Ntot, ns1, nu1, ns0, nu0) {
    bbll(par[c(1,2)], Ntot[, "ns1"] + Ntot[, "ns0"], ns1+ns0)+
    bbll(par[c(7,8)], Ntot[, "nu1"],nu1)+
    bbll(par[c(5,6)],Ntot[, "nu0"], nu0)
  }

.ll4 = function(par, Ntot, ns1, nu1, ns0, nu0) {
  bbll(par[c(1,2)],(Ntot[,c("ns1")]),ns1)+
    bbll(par[c(7,8)], rowSums(Ntot[,c("nu0","nu1")]), nu0 + nu1)+
    bbll(par[c(3,4)],(Ntot[,c("ns0")]),ns0)#+log(ns1/Ntot[, "ns1"]-nu1/Ntot[, "nu1"]-ns0/Ntot[, "ns0"]+nu0/Ntot[, "nu0"]+1)+log(ns1/Ntot[, "ns1"]-nu1/Ntot[, "nu1"]+1)
  }

.ll5 = function(par, Ntot, ns1, nu1, ns0, nu0) {
  bbll(par[c(5,6)],rowSums(Ntot[,c("ns1","nu1")]),ns1+nu1)+
    bbll(par[c(7,8)], rowSums(Ntot[,c("ns0","nu0")]), ns0 + nu0)
}
.ll6 = function(par, Ntot, ns1, nu1, ns0, nu0) {
  bbll(par[c(5,6)],rowSums(Ntot[,c("ns1","nu1")]),ns1+nu1)+
    bbll(par[c(7,8)], (Ntot[,c("nu0")]), nu0) +
    bbll(par[c(3,4)], (Ntot[,c("ns0")]), ns0)
}
.ll7 = function(par, Ntot, ns1, nu1, ns0, nu0) {
    bbll(par[c(7,8)], rowSums(Ntot[,c("nu1","ns1")])+rowSums(Ntot[,c("nu0","ns0")]), ns1+nu1+nu0+ns0)
}

.ll8 = function(par, Ntot, ns1, nu1, ns0, nu0) {
  bbll(par[c(7,8)],rowSums(Ntot[,c("ns1","ns0")]),ns1+ns0)+
    bbll(par[c(1,2)], rowSums(Ntot[,c("nu1","nu0")]), nu0+nu1)#+log(nu1/Ntot[, "nu1"]+nu0/Ntot[, "nu0"]+1)+log(-ns1/Ntot[, "ns1"]+nu1/Ntot[, "nu1"]+ns0/Ntot[, "ns0"]-nu0/Ntot[, "nu0"]+1)
}
# .ll9 = function(par, Ntot, ns1, nu1, ns0, nu0) {
#   bbll(par[c(1,2)],(Ntot[,c("nu1")]),nu1)+
#     bbll(par[c(3,4)], rowSums(Ntot[,c("ns1","nu0","ns0")]), ns1+nu0+ns0)+log(nu1/Ntot[, "nu1"]-ns1/Ntot[, "ns1"]+1)
# }
# .ll10 = function(par,Ntot,ns1,nu1,ns0,nu0){
#   bbll(par[c(1,2)],(Ntot[,c("nu1")]),nu1)+bbll(par[c(3,4)], (Ntot[,c("ns1")]), ns1)+
#     bbll(par[c(5,6)],(Ntot[,c("ns0")]),ns0)+bbll(par[c(7,8)], (Ntot[,c("nu0")]), nu0)+log(nu1/Ntot[, "nu1"]-ns1/Ntot[, "ns1"]+1)
# }
# .ll11 = function(par,Ntot,ns1,nu1,ns0,nu0){
#   bbll(par[c(1,2)],(Ntot[,c("nu1")]),nu1)+bbll(par[c(3,4)], (Ntot[,c("ns1")]), ns1)+
#     bbll(par[c(7,8)],(Ntot[,c("ns0")]+Ntot[,c("nu0")]),ns0+nu0)+log(nu1/Ntot[, "nu1"]-ns1/Ntot[, "ns1"]+1)
# }
.llm1 = function(par,Ntot,ns,nu){
  const(Ntot[,"ns"],ns)+const(Ntot[,"nu"],nu)+
    bbll(par[c(1,2)],Ntot[,"nu"],nu)+
    bbll(par[c(3,4)],Ntot[,"ns"],ns)
}
.llm0 = function(par,Ntot,ns,nu){
  const(Ntot[,"ns"],ns)+const(Ntot[,"nu"],nu)+
    bbll(par[c(1,2)],Ntot[,"ns"]+Ntot[,"nu"],ns+nu)
}
.cllm0m1 = function(par,Ntot,ns,nu){
  cbind(llm0(par,Ntot,ns,nu),llm1(par,Ntot,ns,nu))
}
.sumcllm0m1 = function(..., inds, pi_est) {
  sum(inds * t(log(pi_est) + t(cll(...))))
}

#' log-likelihood of each component
#' @details Calcualtes the  log-likelihood for each observation at each model component
#' @param par vector of parameters.
#' @param Ntot \code{matrix} integer vector of total trials. One row per subject. Should have four columns named "ns1" "ns0" "nu1" and "nu0"
#' @param ns1 \code{numeric} integer vector of successes in condition 1 treatment s.
#' @param nu1 \code{numeric} integer vector of successes in condition 1 treatment u.
#' @param ns0 \code{numeric} integer vector of successes in condition 0 treatment s.
#' @param nu0 \code{numeric} integer vector of successes in condition 0 treatment u.
#' @seealso \link{bbll} \link{MIMOSA2}
cll = function(par, Ntot, ns1, nu1, ns0, nu0) {
  cbind(
    .ll1(par, Ntot=Ntot, ns1=ns1, nu1=nu1, ns0=ns0, nu0=nu0),
   .ll2(par, Ntot=Ntot, ns1=ns1, nu1=nu1, ns0=ns0, nu0=nu0),
   .ll3(par, Ntot=Ntot, ns1=ns1, nu1=nu1, ns0=ns0, nu0=nu0),
   .ll4(par,Ntot=Ntot, ns1=ns1, nu1=nu1, ns0=ns0, nu0=nu0),
   .ll5(par,Ntot=Ntot, ns1=ns1, nu1=nu1, ns0=ns0, nu0=nu0),
   .ll6(par,Ntot=Ntot, ns1=ns1, nu1=nu1, ns0=ns0, nu0=nu0),
    .ll7(par,Ntot=Ntot, ns1=ns1, nu1=nu1, ns0=ns0, nu0=nu0),
   .ll8(par,Ntot=Ntot, ns1=ns1, nu1=nu1, ns0=ns0, nu0=nu0)
  )
}

#'Calculate the complete data log-likelihood across all observations.
#'
#' @param ... all model parameters and data (par, Ntot, ns0,ns1,nu0,nu1)
#' @param inds \code{matrix} of type \code{numeric} represents the max(z's), i.e. the class assignments of each observation to each component.
#' @param pi_est \code{numeric} the mixing proportions.
#' @seealso \link{bbll} \link{MIMOSA2}
sumcll = function(..., z, pi_est) {
  params = as.list(...)
  mus = unlist(params[c(1,3,5,7)])
  mus = invlogit(mus)
  -(sum(z * t(log1p(pi_est) + t(cll(...)))))+
    dgamma((params[[2]]),shape=11/4,rate=0.5,log=TRUE)+
    dgamma((params[[6]]),shape=11/4,rate=0.5,log=TRUE)+
    dgamma((params[[8]]),shape=11/4,rate=0.5,log=TRUE)+
    dgamma((params[[4]]),shape=11/4,rate=0.5,log=TRUE)-
  dnorm(mus[1]-mus[3]-mus[2]+mus[4],mean=1e-4,sd=5e-4,log=TRUE)-
    dnorm(mus[1]-mus[3],mean=1e-4,sd=5e-4,log=TRUE)-
    dnorm(mus[2]-mus[4],mean=1e-4,sd=5e-4,log=TRUE)-
    dnorm(mus[1]-mus[2],mean=1e-4,sd=5e-4,log=TRUE)-
    dnorm(mus[3]-mus[4],mean=1e-4,sd=5e-4,log=TRUE)
}


#' logit
#' @description Return the logit of p
#' @details  No bounds checking is done on the input
#' @param p \code{numeric} between (0,1)
#'
#' @return logit(p)
#' @export
#' @examples
#' logit(0.5)
logit = function (p) {
  log(p/(1 - p))
}

#' invlogit
#'
#' @description Return the inverse logit of x.
#' @details No bounds checking is done on the input
#' @param x \code{numeric} between (-Inf,Inf)
#'
#' @return \code{invlogit(x)}
#' @export
#' @examples
#' invlogit(0)
invlogit = function(x) {
  1/(1 + exp(-x))
}
