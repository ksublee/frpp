library(frpp)
source(file.path("helpers","kernel-api.R"))
if (frpp_has_parameter_api) local({
 model <- fspec(methods::new("hspec",mu=1,alpha=0,beta=1),
  residual_law(c(rate=1),dexp,pexp,qexp,rexp))
 dt <- c(0,.2,.3,.5)
 # Supported numeric fits must use their own public registration/restoration.
 trace("fitted_hspec",where=asNamespace("emhawkes"),print=FALSE,
  tracer=quote(stop("legacy restoration invoked")))
 on.exit(if(inherits(get("fitted_hspec",asNamespace("emhawkes")),"functionWithTrace"))
  untrace("fitted_hspec",where=asNamespace("emhawkes")),add=TRUE)
 for (provider in c("stats",if(requireNamespace("maxLik",quietly=TRUE))"maxLik"))
  for (estimate_residual in c(FALSE,TRUE)) {
   fixed <- if(estimate_residual)c(TRUE,TRUE,TRUE,FALSE) else c(FALSE,TRUE,TRUE)
   fit <- ffit(model,dt,optimizer=provider,fixed=fixed,lambda_component0=0,
    estimate_residual=estimate_residual)
   replay <- do.call(logLik,c(list(object=fitted_fspec(fit),inter_arrival=dt),finitial_state(fit)))
   stopifnot(identical(as.numeric(logLik(fit)),as.numeric(replay)))
  }
 callback <- model;callback@kernel@rmark <- function()1
 error <- tryCatch(ffit(callback,dt,optimizer="stats",fixed=c(FALSE,TRUE,TRUE),
   lambda_component0=0),error=identity)
 stopifnot(inherits(error,"error"),grepl("legacy restoration invoked",conditionMessage(error),fixed=TRUE))
 untrace("fitted_hspec",where=asNamespace("emhawkes"))
 # Shared kernel/residual names retain exactly one legacy start warning.
 shared <- fspec(model@kernel,residual_law(c(mu1=.5),
   density=function(x,mu1,log=FALSE)dexp(x,rate=mu1,log=log)))
 warnings <- character()
 fit <- withCallingHandlers(ffit(shared,dt,optimizer="stats",estimate_residual=TRUE,
   fixed=c(TRUE,TRUE,TRUE),lambda_component0=0),warning=function(w) {
   warnings <<- c(warnings,conditionMessage(w));invokeRestart("muffleWarning")
 })
 stopifnot(sum(startsWith(warnings,"Shared parameters have different starting values"))==1L,
   identical(fitted_fspec(fit)@residual$param,c(mu1=1)))
 cat("Numeric fit restoration uses the public contract; callback fallback and sharing warnings preserved.\n")
})
