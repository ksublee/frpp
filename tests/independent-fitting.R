library(frpp)
source(file.path("helpers","kernel-api.R"))
model <- fspec(methods::new("hspec",mu=1,alpha=0,beta=1),
  residual_law(c(rate=1),dexp,pexp,qexp,rexp))
dt <- c(0,.2,.3,.5)
fails <- function(expr,text) {
  error <- tryCatch(force(expr),error=identity)
  stopifnot(inherits(error,"error"),grepl(text,conditionMessage(error),fixed=TRUE))
}
if (!frpp_has_parameter_api) {
  fails(ffit(model,dt,fit_backend="independent"),"parameter_plan API")
  a <- ffit(model,dt,optimizer="stats",lambda_component0=0,fixed=c(FALSE,TRUE,TRUE))
  b <- ffit(model,dt,optimizer="stats",lambda_component0=0,fixed=c(FALSE,TRUE,TRUE),fit_backend="auto")
  stopifnot(identical(coef(a),coef(b)),identical(vcov(a),vcov(b)),
    identical(b$fitting$used,"bridge"),!is.null(b$fitting$fallback_reason))
} else local({
  # Reuse the fixed observations and independent hfit reference, never a newly
  # generated expected estimate. Tests must select the new installed code.
  environment <- new.env(parent=globalenv())
  environment$ffit <- function(...) frpp::ffit(...,fit_backend="independent")
  sys.source("fitting-equivalence.R",envir=environment)
  environment <- new.env(parent=globalenv())
  environment$ffit <- function(...) frpp::ffit(...,fit_backend="auto")
  sys.source("fitting.R",envir=environment)
  environment <- new.env(parent=globalenv())
  environment$ffit <- function(...) frpp::ffit(...,fit_backend="independent")
  sys.source("fitting-session.R",envir=environment)
  # Prove supported fitting does not delegate optimization to hfit.
  trace("hfit",signature="hspec",where=asNamespace("emhawkes"),print=FALSE,
    tracer=quote(stop("hfit delegation forbidden")))
  traced <- TRUE
  on.exit(if(traced)untrace("hfit",signature="hspec",where=asNamespace("emhawkes")),add=TRUE)
  fit <- ffit(model,dt,optimizer="stats",lambda_component0=0,
    fixed=c(FALSE,TRUE,TRUE),fit_backend="independent")
  stopifnot(identical(fit$fitting$used,"independent"),
    inherits(fit$fit,"frpp_stats_fit"),is.null(attr(fit$fit,"emhawkes")))
  untrace("hfit",signature="hspec",where=asNamespace("emhawkes"))
  traced <- FALSE
  # Exact candidate validation order, derivative and curvature comparison.
  for (provider in c("stats",if(requireNamespace("maxLik",quietly=TRUE))"maxLik")) {
    for (method in if(provider=="stats")c("BFGS","CG","Nelder-Mead","L-BFGS-B") else "BFGS") {
      evaluations <- list()
      validator <- function(param) {evaluations[[length(evaluations)+1L]] <<- param;TRUE}
      args <- list(object=model,inter_arrival=dt,optimizer=provider,method=method,
        fixed=c(FALSE,TRUE,TRUE),lambda_component0=.1,N0=2,Nc0=3,valid=validator)
      if(method=="L-BFGS-B")args <- c(args,list(lower=c(.01,-Inf,-Inf),upper=rep(10,3)))
      set.seed(971);a <- suppressWarnings(do.call(ffit,args));rng <- .Random.seed;calls <- evaluations
      evaluations <- list();set.seed(971)
      b <- suppressWarnings(do.call(ffit,c(args,list(fit_backend="independent"))))
      H <- function(fit) if(provider=="stats")fit$fit$hessian else maxLik::hessian(fit$fit)
      stopifnot(identical(coef(a),coef(b)),identical(as.numeric(logLik(a)),as.numeric(logLik(b))),
        identical(vcov(a),vcov(b)),identical(H(a),H(b)),identical(finitial_state(a),finitial_state(b)),
        identical(fdiagnostics(a),fdiagnostics(b)),identical(calls,evaluations),identical(rng,.Random.seed))
    }
  }
  callback <- model;callback@kernel@rmark <- function()1
  fails(ffit(callback,dt,fit_backend="independent"),"event callbacks")
  a <- ffit(callback,dt,optimizer="stats",fixed=c(FALSE,TRUE,TRUE),lambda_component0=0)
  b <- ffit(callback,dt,optimizer="stats",fixed=c(FALSE,TRUE,TRUE),lambda_component0=0,fit_backend="auto")
  stopifnot(identical(coef(a),coef(b)),identical(b$fitting$used,"bridge"),
    identical(b$fitting$fallback_reason,"event callbacks require the bridge"),
    is.null(a$fitting$fallback_reason))
  # Global shared-name conflict is reported once, including all-fixed fits.
  shared <- fspec(model@kernel,residual_law(c(mu1=.5),
    density=function(x,mu1,log=FALSE)dexp(x,rate=mu1,log=log)))
  warnings <- character()
  fit <- withCallingHandlers(ffit(shared,dt,optimizer="stats",estimate_residual=TRUE,
    fixed=rep(TRUE,3),lambda_component0=0,fit_backend="independent"),warning=function(w) {
      warnings <<- c(warnings,conditionMessage(w));invokeRestart("muffleWarning")
    })
  stopifnot(sum(startsWith(warnings,"Shared parameters have different starting values"))==1L,
    identical(fitted(fit)@residual$param,c(mu1=1)))
  # Non-square mappings, eta/observed marks and named residual coordinates.
  mapped <- fspec(methods::new("hspec",mu=c(.4,.6),alpha=matrix(.03,2,4),
    beta=matrix(1.2,2,4),eta=matrix(.01,2,4),
    type_col_map=list(c(3L,1L),c(2L,3L,4L))),
    residual_law(c(shape=2,rate=2),dgamma,pgamma,qgamma,rgamma))
  times <- c(0,.2,.8,.4,.3,.7,.6,.2,.9,.5,.4,.8)
  types <- c(0,1,2,2,1,2,1,1,2,2,1,2)
  marks <- c(0,1,2,3,1,2,1,2,1,3,1,2)
  for (provider in c("stats",if(requireNamespace("maxLik",quietly=TRUE))"maxLik"))
    for (reduced in c(TRUE,FALSE)) {
      fixed <- rep(TRUE,if(reduced)7 else 28)
      fixed[1:2] <- FALSE
      args <- list(object=mapped,inter_arrival=times,type=types,mark=marks,
        optimizer=provider,reduced=reduced,fixed=fixed,estimate_residual=TRUE,
        N0=c(2,3),Nc0=c(4,5),lambda_component0=matrix(.1,2,4))
      a <- suppressWarnings(suppressMessages(do.call(ffit,args)))
      b <- suppressWarnings(suppressMessages(do.call(ffit,c(args,list(fit_backend="independent")))))
      stopifnot(identical(coef(a),coef(b)),identical(as.numeric(logLik(a)),as.numeric(logLik(b))),
        identical(vcov(a),vcov(b)),identical(finitial_state(a),finitial_state(b)),
        identical(fdiagnostics(a),fdiagnostics(b)),
        identical(fkernel(fitted(a))@eta,fkernel(fitted(b))@eta),
        identical(fkernel(fitted(a))@type_col_map,fkernel(fitted(b))@type_col_map))
    }
  fails(ffit(model,c(1,.2),fit_backend="independent"),"initial zero")
  fails(ffit(model,dt,lambda_component0=c(1,2),fit_backend="independent"),"kernel matrix dimensions")
  fails(ffit(model,dt,optimizer="stats",fit_backend="independent",valid=function(param)
    stop("validator exception")),"validator exception")
  cat("Independent fitting: exact references, provider paths, validation calls, saved sessions and fallbacks passed.\n")
})
