% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/likelihoods.R
\name{bbll}
\alias{bbll}
\title{Beta-binomial log likelihood}
\usage{
bbll(par, n, k)
}
\arguments{
\item{par}{vector of beta-binomial parameters \code{par[1] = logit(a / (a+b)), par[2] = log(a+b)}.}

\item{n}{\code{numeric} integer, the total trials. Can be a vector.}

\item{k}{\code{numeric} integer, the number of successes. Can be a vector.}
}
\value{
vector of beta-binomial log likelihoods for data \code{n,k} with parameters \code{par}
}
\description{
Beta-binomial log likelihood
}
\examples{
bbll(c(logit(0.5),log(100)),100,50)
}
