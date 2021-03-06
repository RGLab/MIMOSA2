% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/likelihoods.R
\name{cll}
\alias{cll}
\title{log-likelihood of each component}
\usage{
cll(par, Ntot, ns1, nu1, ns0, nu0)
}
\arguments{
\item{par}{vector of parameters.}

\item{Ntot}{\code{matrix} integer vector of total trials. One row per subject. Should have four columns named "ns1" "ns0" "nu1" and "nu0"}

\item{ns1}{\code{numeric} integer vector of successes in condition 1 treatment s.}

\item{nu1}{\code{numeric} integer vector of successes in condition 1 treatment u.}

\item{ns0}{\code{numeric} integer vector of successes in condition 0 treatment s.}

\item{nu0}{\code{numeric} integer vector of successes in condition 0 treatment u.}
}
\description{
log-likelihood of each component
}
\details{
Calcualtes the  log-likelihood for each observation at each model component
}
\seealso{
\link{bbll} \link{MIMOSA2}
}
