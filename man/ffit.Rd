\name{ffit}
\alias{ffit}
\alias{fitted_fspec}
\alias{finitial_state}
\alias{fdiagnostics}
\alias{coef.ffit}
\alias{vcov.ffit}
\alias{logLik.ffit}
\alias{nobs.ffit}
\alias{fitted.ffit}
\alias{summary.ffit}
\alias{print.ffit}
\alias{print.summary.ffit}
\title{Fit a Flexible Residual Point Process}
\description{Estimate a composition model using an FRPP likelihood. Reuses the
public emhawkes parameter map and optimizer adapters and restores an fspec.}
\usage{
ffit(object, inter_arrival, type = NULL, mark = NULL,
     N = NULL, Nc = NULL, N0 = NULL, Nc0 = NULL,
     lambda_component0 = "auto", estimate_residual = FALSE,
     optimizer = c("maxLik", "stats"), method = "BFGS",
     backend = c("bridge", "state"), reduced = TRUE,
     fixed = NULL, constraints = NULL, constraint = NULL,
     valid = NULL, ...)
fitted_fspec(object)
finitial_state(object)
fdiagnostics(object)
}
\arguments{
\item{object}{An \code{fspec} for fitting, or an \code{ffit} result for the helpers.}
\item{inter_arrival, type, mark}{Observed intervals, event types and optional marks,
including an initial zero row.}
\item{N, Nc}{Optional cumulative count and mark-total paths.}
\item{N0, Nc0}{Initial counts and mark totals.}
\item{lambda_component0}{Initial kernel components. \code{"auto"} selects zero;
specify a scalar or conformable matrix when known. Hawkes mean initialization is not used.}
\item{estimate_residual}{Keep the residual distribution fixed by default.
Set \code{TRUE} to include its named parameters in the optimization.}
\item{optimizer}{Optional maxLik provider, or R's stats provider through the public
emhawkes adapter. The stats provider does not require maxLik.}
\item{method}{Optimization method, default BFGS. The stats provider also supports
Nelder-Mead, CG and L-BFGS-B.}
\item{backend}{FRPP likelihood backend. \code{state} retains the numeric-kernel
restrictions described in \code{\link{logLik.fspec}}.}
\item{reduced}{Share entries having equal starting values within each numeric
kernel matrix. Function parameters sharing a name are shared under either setting.}
\item{fixed}{Logical mask or parameter positions, passed to the optimizer.
With stats, names are also supported. The order follows emhawkes: mu, alpha,
beta, eta, impact, mark density, residual density/CDF/quantile; duplicated
function parameter names occur once. Fixed residual parameters are excluded
entirely when \code{estimate_residual = FALSE}.}
\item{constraints, constraint}{Alternative names for optimization constraints.
Supply at most one. Linear inequalities are \code{list(ineqA, ineqB)} with
\eqn{A\theta+B\ge0}; columns follow the full parameter order, including fixed
optimization parameters. Stats uses constrOptim and requires a strictly interior
start for movable inequalities; it does not support equalities.}
\item{valid}{Optional function of the named full parameter vector, returning
one nonmissing logical value. FALSE rejects a proposal with minus infinity.
Fixed raw residual parameters are also included. Validator and callback errors
propagate. Supply distribution-specific domain checks when needed.}
\item{...}{Named optimizer arguments accepted by the public hfit adapter,
including \code{control}, \code{grad}, \code{hess}, and stats box bounds.
For stats use \code{control$maxit}; maxLik uses \code{control$iterlim}.
FRPP objective, analytic Hawkes gradient and Hawkes final-Hessian selection
cannot be overridden through these arguments.}
}
\details{
The new fitting API delegates parameter sharing, optimization and kernel
reconstruction to public emhawkes functions; its objective always evaluates the
FRPP model. No internal emhawkes functions or Hawkes curvature stabilization are used.
Likelihood evaluation and default BFGS fitting do not generate random values.
Stochastic optimizers can consume randomness. User callbacks must be deterministic.

Residual parameters are fixed by default to avoid silently estimating two
confounded scales. Joint estimation requires an identifiable model: fixing or
normalizing the relevant scale remains the user's modeling choice. This API does
not automatically diagnose all forms of nonidentification or normalize arbitrary laws.
Matrix sharing is preserved on reconstruction; residual parameters shared among
density, CDF, quantile and generator functions are updated together.

The result supports coef, logLik, nobs, vcov, summary and fitted.
\code{fitted_fspec} and \code{fitted} return the estimated composition model;
\code{finitial_state} returns the initialization arguments for likelihood replay
or simulation. Save the result with saveRDS and restore it with readRDS.
Captured callback environments remain references; external state is not frozen.

Wald covariance is withheld on optimizer failure, missing or nonpositive
curvature, near detected constraint boundaries, and equality constraints.
The original provider result remains in \code{fit$fit}. A positive covariance is
reported as \code{available_unvalidated}: numerical step sensitivity, nonsmooth
densities, weak identification and actual confidence-interval coverage have not
been validated. A supplied analytic Hessian can avoid finite-difference error
for a smooth model. No eigenvalue clipping or absolute-variance correction is used.
The tzexp transition may give a nonsmooth optimum; consider profile likelihood
or bootstrap separately rather than interpreting a success code as valid inference.
}
\value{\code{ffit} returns a list of class \code{ffit}, including the provider fit,
estimated fspec, coefficients, covariance, likelihood, initialization and diagnostics.
\code{fitted_fspec} returns an fspec, \code{finitial_state} a list, and
\code{fdiagnostics} a list containing convergence and Wald-availability information.
Likelihood degrees of freedom count only active optimization parameters;
nobs excludes the initial row. Unavailable free-parameter covariance entries are NA.}
\examples{
model <- fspec(methods::new("hspec", mu = 1, alpha = 0, beta = 1),
  residual_law(c(rate = 1), stats::dexp, stats::pexp, stats::qexp, stats::rexp))
fit <- ffit(model, c(0, .2, .3, .5), optimizer = "stats",
  fixed = c(FALSE, TRUE, TRUE), lambda_component0 = 0)
coef(fit)
summary(fit)
fdiagnostics(fit)
estimated <- fitted(fit)
do.call(logLik, c(list(object = estimated, inter_arrival = c(0, .2, .3, .5)),
  finitial_state(fit)))
}
