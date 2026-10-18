\name{fkernel}
\alias{fkernel}
\alias{fresidual}
\alias{fkernel<-}
\alias{fresidual<-}
\title{Access Components of an FRPP Specification}
\description{Read or replace the kernel and residual law of an S4 \code{fspec}.
Legacy S3 specifications are converted before access. Replacement validates
the whole specification without executing its callbacks.}
\usage{
fkernel(object)
fkernel(object) <- value
fresidual(object)
fresidual(object) <- value
}
\arguments{
\item{object}{An \code{fspec} or a compatible legacy FRPP specification.}
\item{value}{For \code{fkernel}, an \code{hspec} without residual callbacks.
For \code{fresidual}, a law constructed by \code{residual_law}.}
}
\value{The requested component. Replacement returns an updated S4 specification.}
\seealso{\code{\link{fspec}}, \code{\link{residual_law}}}
\examples{
m <- fspec(methods::new("hspec", mu = .4, alpha = .1, beta = 1),
  residual_law(c(rate = 1), stats::dexp, stats::pexp, stats::qexp, stats::rexp))
fkernel(m)
fresidual(m)$name
fresidual(m) <- residual_law(c(shape = 2, rate = 2), stats::dgamma,
  stats::pgamma, stats::qgamma, stats::rgamma)
methods::validObject(m)
}
