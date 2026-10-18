\name{residual_law}
\alias{residual_law}
\title{Bundle Residual Distribution Functions}
\description{Store a named parameter vector and residual distribution functions
without evaluating the functions or drawing random numbers.}
\usage{
residual_law(param = numeric(), density = NULL, cdf = NULL,
             quantile = NULL, rng = NULL, name = "custom")
}
\arguments{
\item{param}{Finite numeric vector with unique, nonempty parameter names.}
\item{density, cdf, quantile, rng}{Distribution functions, or \code{NULL}.
Use R distribution signatures: the variate, probability or sample count is the
first argument. Named parameters are matched to explicit formal arguments as in
\code{emhawkes::set_residual}. Log-density and tail options are supported there.}
\item{name}{A descriptive label, not a verified distribution-family identifier.}
}
\details{At least a density or random generator is required. Simulation needs a
generator; likelihood evaluation needs a density, and multivariate likelihoods
also need a CDF. Residual inference needs a CDF and quantile function.
The constructor does not establish consistency of user-supplied functions,
normalization of residual scale, existence of moments, or parameter identification.}
\value{An object of class \code{frpp_residual_law}.}
