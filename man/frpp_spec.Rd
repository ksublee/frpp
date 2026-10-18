\name{frpp_spec}
\alias{fspec}
\alias{as_fspec}
\alias{print.fspec}
\alias{frpp_spec}
\alias{as_frpp}
\alias{as_legacy_hspec}
\alias{print.frpp_spec}
\title{Compose a Kernel Specification and a Residual Law}
\description{Create an experimental FRPP object containing an emhawkes kernel and
a separate residual-law object. The result does not inherit from \code{hspec}.}
\usage{
fspec(kernel, residual)
frpp_spec(kernel, residual)
as_fspec(object, name = "legacy")
as_frpp(object, name = "legacy")
as_legacy_hspec(object)
}
\arguments{
\item{kernel}{An \code{hspec} without residual-distribution slots.}
\item{residual}{A law created by \code{residual_law}.}
\item{object}{For conversion, an existing FRPP \code{hspec} or a \code{frpp_spec}.
For the legacy adapter, a \code{frpp_spec}.}
\item{name}{A descriptive name for a converted legacy law.}
}
\details{The current model refreshes every dimension's candidate residual after
every event in the combined process. It does not retain separate residual clocks.
Conversion preserves callback functions and their default bindings without
executing them. Captured environments remain references, as in emhawkes.
The primary constructor is \code{fspec}, matching the \code{hspec} naming convention.
\code{frpp_spec} and \code{as_frpp} remain compatibility names. New objects have
classes \code{c("fspec", "frpp_spec")}; saved objects from the first prototype
can still use the old simulation and likelihood methods.
The adapter is transitional: it restores the old FRPP representation to use
existing public emhawkes operations. It is not a pure Hawkes specification.}
\value{\code{frpp_spec} and \code{as_frpp} return a composition object.
\code{as_legacy_hspec} returns an \code{hspec} carrying the residual callbacks.}
\seealso{\code{\link{residual_law}}, \code{\link{simulate.frpp_spec}}}
\examples{
kernel <- methods::new("hspec", mu = .3, alpha = .2, beta = 1)
law <- residual_law(c(shape = 2, rate = 2), stats::dgamma,
                    stats::pgamma, stats::qgamma, stats::rgamma,
                    name = "Gamma(2,2)")
model <- fspec(kernel, law)
model
stopifnot(!methods::is(model, "hspec"))
}
