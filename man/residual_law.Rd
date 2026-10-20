\name{residual_law}
\alias{residual_law}
\title{Bundle Residual Distribution Functions}
\description{Store a named parameter vector and residual distribution functions
without evaluating the functions or drawing random numbers.}
\usage{
residual_law(param = numeric(), dresidual = NULL, presidual = NULL,
             qresidual = NULL, rresidual = NULL, name = "custom",
             density = NULL, cdf = NULL, quantile = NULL, rng = NULL)
}
\arguments{
\item{param}{Finite numeric vector with unique, nonempty parameter names.}
\item{dresidual, presidual, qresidual, rresidual}{Residual density, cumulative
distribution, quantile and random-generation functions, or \code{NULL}.
Use R distribution signatures: the variate, probability or sample count is the
first argument. Named parameters are matched to explicit formal arguments as in
\code{emhawkes::set_residual}. Log-density and tail options are supported there.}
\item{density, cdf, quantile, rng}{Compatibility aliases for \code{dresidual},
\code{presidual}, \code{qresidual} and \code{rresidual}, respectively.
Supply at most one non-NULL function per pair; specifying both is an error,
even when they are the same function. No deprecation warning is emitted.}
\item{name}{A descriptive label, not a verified distribution-family identifier.}
}
\details{Supply only the functions needed for your operation. Construction
requires at least \code{dresidual} or \code{rresidual}. Simulation needs only
\code{rresidual}; univariate likelihood and fitting need only \code{dresidual}.
Multivariate likelihood and fitting additionally need \code{presidual}.
\code{infer = TRUE} needs \code{dresidual}, \code{presidual} and \code{qresidual},
and does not require \code{rresidual}. Missing functions are reported when an
operation requires them, before callbacks or random draws.
The original positional order of the first six arguments is unchanged.
Stored law objects keep their existing representation, and callback functions,
environments and defaults are preserved.
The constructor does not establish consistency of user-supplied functions,
normalization of residual scale, existence of moments, or parameter identification.}
\value{An object of class \code{frpp_residual_law}.}
\examples{
# Simulation needs only a residual generator.
law <- residual_law(c(rate = 1), rresidual = stats::rexp)
model <- fspec(methods::new("hspec", mu = 1, alpha = 0, beta = 1), law)
fsim(model, size = 4, seed = 73, lambda_component0 = 0)

# Density alone is enough for a univariate likelihood.
law <- residual_law(c(rate = 1), dresidual = stats::dexp)
model <- fspec(methods::new("hspec", mu = 1, alpha = 0, beta = 1), law)
stats::logLik(model, c(0, .2, .3), lambda_component0 = 0)
}
