\name{simulate.frpp_spec}
\alias{simulate.frpp_spec}
\alias{fsim}
\alias{logLik.frpp_spec}
\alias{simulate.fspec}
\alias{logLik.fspec}
\alias{print.frpp_path}
\title{Simulate and Evaluate the Composition Prototype}
\description{Use a compatibility bridge, or assemble a numeric FRPP likelihood
from public emhawkes kernel replay output and the residual law.}
\usage{
fsim(object, size = 100, lambda_component0 = "auto", N0 = NULL,
    Nc0 = NULL, seed = NULL, ...)
\method{simulate}{fspec}(object, nsim = 1, seed = NULL, size = 100,
    lambda_component0 = "auto", N0 = NULL, Nc0 = NULL, ...)
\method{logLik}{fspec}(object, inter_arrival, type = NULL, mark = NULL,
    N = NULL, Nc = NULL, N0 = NULL, Nc0 = NULL,
    lambda_component0 = "auto", infer = FALSE,
    backend = c("bridge", "state"), ...)
}
\arguments{
\item{object}{An \code{fspec}. Old \code{frpp_spec} objects are also supported.}
\item{nsim}{One path is supported. Use \code{size} for its length.}
\item{seed}{Optional seed. When supplied, the caller's RNG state is restored.}
\item{size}{Number of output rows including the initial row.}
\item{lambda_component0}{Initial kernel components. Automatic initialization
uses zero for the residual-distribution model. Supply an explicit state when known.
A Hawkes stationary mean is not automatically applied to FRPP.}
\item{N0, Nc0}{Initial counts and mark totals, passed to emhawkes.}
\item{inter_arrival, type, mark}{Observations starting with an initial zero row.}
\item{N, Nc}{Optional observed count paths.}
\item{infer}{Return driver states, likelihood and residual diagnostics.}
\item{backend}{\code{bridge} delegates to the existing FRPP likelihood.
\code{state} independently assembles the residual likelihood from the
kernel's public replay output. This requires the emhawkes development kernel
API; with earlier versions, use \code{bridge}.}
\item{...}{Additional public emhawkes arguments for the bridge.
The state prototype accepts none.}
}
\details{Simulation retains the existing R FRPP loop and requires equal decay
rates within each row. Every global event refreshes all candidate residuals.
The state backend currently supports numeric kernels, numeric mark effects and
kernel mappings, with equal decay rates within each row and no impact or
mark-density callbacks. Callback specifications remain available through the
bridge. Kernel states and their integrals are drivers; with a non-exponential
residual law they are not automatically the FRPP conditional intensities.

For event type j in a global inter-event interval, its contribution is
\deqn{\log v_j + \log f(A_j) + \sum_{k \ne j}\log S(A_k),}
where v is the pre-event kernel driver and A is its interval integral.
The public replay call returns kernel-only interval integrals and pre-event
sums. The FRPP layer adds baseline contributions. Full component histories are
retained only with \code{infer = TRUE}; interval integrals and sums use memory
proportional to path length times the number of event types. This backend does
not evaluate a Hawkes likelihood. \code{fsim()} is the primary simulation
function; the \code{simulate()} methods call it and return the same path.

The existing \code{emhawkes.engine} option retains its original meaning. In
particular, the bridge does not support a forced C++ FRPP likelihood.}
\value{Simulation returns a \code{frpp_path} list containing its composition model
and the original path fields other than the transitional \code{hspec} field.
Likelihood evaluation returns the existing named numeric value, or an inference
list. The state backend additionally returns \code{integrated_driver}.
No fitted degrees of freedom or confidence intervals are supplied.}
\examples{
kernel <- methods::new("hspec", mu = .3, alpha = .2, beta = 1)
model <- fspec(kernel, residual_law(c(shape = 2, rate = 2),
    stats::dgamma, stats::pgamma, stats::qgamma, stats::rgamma))
path <- fsim(model, size = 12, seed = 73, lambda_component0 = 0)
stats::logLik(model, path$inter_arrival, path$type, lambda_component0 = 0)
stats::logLik(model, path$inter_arrival, path$type, lambda_component0 = 0,
              backend = "state")
}
