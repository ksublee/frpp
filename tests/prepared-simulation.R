library(frpp)
source(file.path("helpers","kernel-api.R"))
if (frpp_has_plan_api) local({
  model <- fspec(methods::new("hspec",mu=c(.3,.4),alpha=matrix(.02,2,4),
    beta=matrix(c(1,2),2,4),type_col_map=list(c(3L,1L),c(4L,2L))),
    residual_law(c(shape=2,rate=2),dgamma,pgamma,qgamma,rgamma))
  # The fixed specification must not be reconstructed on every event.
  constructions <- 0L
  trace("exp_kernel_spec",where=asNamespace("emhawkes"),print=FALSE,
    tracer=quote(constructions <<- constructions + 1L))
  on.exit(if (inherits(get("exp_kernel_spec",asNamespace("emhawkes")),"functionWithTrace"))
    untrace("exp_kernel_spec",where=asNamespace("emhawkes")),add=TRUE)
  # trace() resolves its counter in the global environment.
  assign("constructions",0L,envir=.GlobalEnv)
  on.exit(rm("constructions",envir=.GlobalEnv),add=TRUE)
  for (size in c(1L,2L,101L)) {
    assign("constructions",0L,envir=.GlobalEnv)
    invisible(fsim(model,size=size,seed=73,lambda_component0=.1,backend="state"))
    stopifnot(get("constructions",.GlobalEnv)==2L)
  }
  untrace("exp_kernel_spec",where=asNamespace("emhawkes"))
  # Default selection must stay on bridge, without preparing a kernel plan.
  trace("exp_kernel_plan",where=asNamespace("emhawkes"),print=FALSE,
    tracer=quote(stop("plan must not be invoked by the default bridge")))
  on.exit(if (inherits(get("exp_kernel_plan",asNamespace("emhawkes")),"functionWithTrace"))
    untrace("exp_kernel_plan",where=asNamespace("emhawkes")),add=TRUE)
  set.seed(531)
  default <- fsim(model,size=17,lambda_component0=.1)
  saved <- .Random.seed; following <- rnorm(3)
  set.seed(531)
  explicit <- fsim(model,size=17,lambda_component0=.1,backend="bridge")
  stopifnot(identical(default,explicit),identical(saved,.Random.seed),
    identical(following,rnorm(3)))
  untrace("exp_kernel_plan",where=asNamespace("emhawkes"))
  # Earlier development versions with step/wait but no plan keep working.
  ns <- asNamespace("frpp")
  original <- get("prepare_frpp_simulation",ns)
  old_api <- original
  scope <- new.env(parent=environment(original))
  scope$getNamespaceExports <- function(ns)
    setdiff(base::getNamespaceExports(ns),"exp_kernel_plan")
  environment(old_api) <- scope
  replace <- function(value) {
    unlockBinding("prepare_frpp_simulation",ns)
    assign("prepare_frpp_simulation",value,ns)
    lockBinding("prepare_frpp_simulation",ns)
  }
  observe <- function() {
    set.seed(532)
    path <- fsim(model,size=101,lambda_component0=.1,backend="state")
    list(path=path,rng=.Random.seed,following=rnorm(3))
  }
  reference <- observe()
  replace(old_api)
  on.exit(replace(original),add=TRUE)
  stopifnot(identical(reference,observe()))
  cat("Prepared-plan reuse, bridge default and earlier event-wise API contracts passed.\n")
})
