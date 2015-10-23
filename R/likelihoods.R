

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



#' Model component 1
#' 
#' @description Represents a vaccine-specific response where the baseline response may be positive, but 
#' lower than the post-vaccine response. Log-Likelihood for the model component where all p's are different.
#' @details The model hyperparameters are different for condition 1 treatment s, condition 0 treatment s, and condition treatment u.
#' @param par vector of parameters. 
#' @param Ntot \code{matrix} integer vector of total trials. One row per subject. Should have four columns named "ns1" "ns0" "nu1" and "nu0" 
#' @param ns1 \code{numeric} integer vector of successes in condition 1 treatment s.
#' @param nu1 \code{numeric} integer vector of successes in condition 1 treatment u.
#' @param ns0 \code{numeric} integer vector of successes in condition 0 treatment s.
#' @param nu0 \code{numeric} integer vector of successes in condition 0 treatment u.
#' @seealso \link{bbll} \link{MIMOSA2}
ll1 = function(par, Ntot, ns1, nu1, ns0, nu0) {
  bbll(par[1:2], Ntot[, "ns1"], ns1) + 
    bbll(par[c(5, 6)],  Ntot[, "nu1"] + 
           Ntot[, "nu0"], nu1 + nu0) + 
    bbll(par[c(3, 4)],  Ntot[, "ns0"], ns0)
}

#' Model component 2
#' @description Represents non-responders.
#' Log-Likelihood where all stims are equal to unstim.
#' @details The model hyperparameters are shared across treatments and conditions and are equal to the control.
#' @param par vector of parameters. 
#' @param Ntot \code{matrix} integer vector of total trials. One row per subject. Should have four columns named "ns1" "ns0" "nu1" and "nu0" 
#' @param ns1 \code{numeric} integer vector of successes in condition 1 treatment s.
#' @param nu1 \code{numeric} integer vector of successes in condition 1 treatment u.
#' @param ns0 \code{numeric} integer vector of successes in condition 0 treatment s.
#' @param nu0 \code{numeric} integer vector of successes in condition 0 treatment u.
#' @seealso \link{bbll} \link{MIMOSA2}
ll2 = function(par, Ntot, ns1, nu1, ns0, nu0) {
  bbll(par[c(5, 6)], rowSums(Ntot), ns1 + nu1 + ns0 + nu0)
}

#' Model component 3
#' @description Represents a non-responder where there may be a response in condition 0 but not in condition 1. 
#' @details Log-Likelihood for non-response where s0 is a response but s1 is not
#' @param par vector of parameters. 
#' @param Ntot \code{matrix} integer vector of total trials. One row per subject. Should have four columns named "ns1" "ns0" "nu1" and "nu0" 
#' @param ns1 \code{numeric} integer vector of successes in condition 1 treatment s.
#' @param nu1 \code{numeric} integer vector of successes in condition 1 treatment u.
#' @param ns0 \code{numeric} integer vector of successes in condition 0 treatment s.
#' @param nu0 \code{numeric} integer vector of successes in condition 0 treatment u.
#' @seealso \link{bbll} \link{MIMOSA2}
ll3 = function(par, Ntot, ns1, nu1, ns0, nu0) {
  bbll(par[c(3, 4)], Ntot[, "ns0"] , ns0) + 
    bbll(par[c(5, 6)],  Ntot[, "ns1"] + 
           Ntot[, "nu1"] + Ntot[, "nu0"], 
         ns1 + nu1 + nu0)
}

#' Model component 4
#' @description Represents non-specific response, where the two conditions are responses but are equal to each other.
#' @details Likelihood for non-specific response where stims are equal and generated from mu[s0]
#' @param par vector of parameters. 
#' @param Ntot \code{matrix} integer vector of total trials. One row per subject. Should have four columns named "ns1" "ns0" "nu1" and "nu0" 
#' @param ns1 \code{numeric} integer vector of successes in condition 1 treatment s.
#' @param nu1 \code{numeric} integer vector of successes in condition 1 treatment u.
#' @param ns0 \code{numeric} integer vector of successes in condition 0 treatment s.
#' @param nu0 \code{numeric} integer vector of successes in condition 0 treatment u.
#' @seealso \link{bbll} \link{MIMOSA2}
ll4 = function(par, Ntot, ns1, nu1, ns0, nu0) {
  bbll(par[c(3, 4)], Ntot[, "ns1"] + Ntot[, "ns0"], ns1 + ns0) + 
    bbll(par[c(5, 6)],  Ntot[, "nu1"] + Ntot[, "nu0"], nu1 + nu0)
}

#' Model component 5
#' @description Likelihood for response where stims at post-vaccine and no response at baseline mu[s0]
#' @details Represents a vaccine specific response where there is no response at condition 0, but a response at condition 1.
#' @param par vector of parameters. 
#' @param Ntot \code{matrix} integer vector of total trials. One row per subject. Should have four columns named "ns1" "ns0" "nu1" and "nu0" 
#' @param ns1 \code{numeric} integer vector of successes in condition 1 treatment s.
#' @param nu1 \code{numeric} integer vector of successes in condition 1 treatment u.
#' @param ns0 \code{numeric} integer vector of successes in condition 0 treatment s.
#' @param nu0 \code{numeric} integer vector of successes in condition 0 treatment u.
#' @seealso \link{bbll} \link{MIMOSA2}
ll5 = function(par, Ntot, ns1, nu1, ns0, nu0) {
  bbll(par[c(1, 2)], Ntot[, "ns1"] , ns1) + 
    bbll(par[c(5, 6)], Ntot[,"ns0"] + Ntot[, "nu1"] + 
           Ntot[, "nu0"], ns0 + nu1 + nu0)
}

#' Complete data log-likelihood
#' @details Calcualtes the complete data log-likelihood for each observation
#' @param par vector of parameters. 
#' @param Ntot \code{matrix} integer vector of total trials. One row per subject. Should have four columns named "ns1" "ns0" "nu1" and "nu0" 
#' @param ns1 \code{numeric} integer vector of successes in condition 1 treatment s.
#' @param nu1 \code{numeric} integer vector of successes in condition 1 treatment u.
#' @param ns0 \code{numeric} integer vector of successes in condition 0 treatment s.
#' @param nu0 \code{numeric} integer vector of successes in condition 0 treatment u.
#' @seealso \link{bbll} \link{MIMOSA2}
cll = function(par, Ntot, ns1, nu1, ns0, nu0) {
  cbind(
    ll1(par, Ntot, ns1, nu1, ns0, nu0),
    ll2(par, Ntot, ns1, nu1, ns0, nu0),
    ll3(par, Ntot, ns1, nu1, ns0, nu0),
    ll4(par, Ntot, ns1, nu1, ns0, nu0),
    ll5(par, Ntot, ns1, nu1, ns0, nu0)
  )
}

#'Sum of the complete data log-likelihood across all observations.
#' 
#' @param ... all model parameters and data (par, Ntot, ns0,ns1,nu0,nu1)
#' @param inds \code{matrix} of type \code{numeric} represents the max(z's), i.e. the class assignments of each observation to each component.
#' @param pi_est \code{numeric} the mixing proportions.
#' @seealso \link{bbll} \link{MIMOSA2}
sumcll = function(..., inds, pi_est) {
  sum(inds * t(log(pi_est) + t(cll(...))))
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