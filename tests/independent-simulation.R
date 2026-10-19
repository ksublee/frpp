library(frpp)
source(file.path("helpers", "kernel-api.R"))
fails <- function(expr, text) {
  error <- tryCatch(force(expr), error = identity)
  stopifnot(inherits(error, "error"), grepl(text, conditionMessage(error), fixed = TRUE))
}
laws <- list(
  residual_law(c(shape = 2, rate = 2), stats::dgamma, stats::pgamma, stats::qgamma, stats::rgamma),
  residual_law(c(rate = 1), stats::dexp, stats::pexp, stats::qexp, stats::rexp),
  residual_law(c(a = .7, ell = 1.6), emhawkes::dtzexp, emhawkes::ptzexp, emhawkes::qtzexp, emhawkes::rtzexp))
kernels <- list(
  methods::new("hspec", mu = .4, alpha = .12, beta = 1.2),
  methods::new("hspec", mu = c(.3, .4), alpha = matrix(.03, 2, 2), beta = matrix(c(1,2,1,2), 2)),
  methods::new("hspec", mu = c(0, .4), alpha = matrix(.04, 2, 3), beta = matrix(c(1,2,1,2,1,2), 2),
    eta = matrix(.02, 2, 3), type_col_map = list(c(3L, 1L), c(2L, 3L))),
  methods::new("hspec", mu = c(.2, .3), alpha = matrix(c(-.02,.03,.06,.01), 2),
    beta = matrix(c(1,2,1,2), 2), type_col_map = list(c(2L,1L), c(1L,2L))))
comparisons <- 0L
for (h in kernels) for (law in laws) for (seed in c(17, 203)) for (size in c(1L, 2L, 17L)) {
  model <- fspec(h, law)
  args <- list(object = model, size = size, lambda_component0 = .1,
    N0 = seq_len(h@dimens), Nc0 = seq_len(h@dimens) + .5)
  set.seed(seed)
  bridge <- do.call(fsim, c(args, list(backend = "bridge")))
  rng <- .Random.seed; following <- stats::rnorm(3)
  set.seed(seed)
  result <- do.call(fsim, c(args, list(backend = if (frpp_has_sim_api) "state" else "auto")))
  stopifnot(identical(result, bridge), identical(rng, .Random.seed),
    identical(following, stats::rnorm(3)))
  if (size == 1L) { set.seed(seed); stopifnot(identical(rng, .Random.seed)) }
  if (frpp_has_sim_api) {
    # Independent existing inference, including component ordering and mapping.
    for (engine in c("r", "cpp")) {
      options(emhawkes.engine = engine)
      inferred <- emhawkes::logLik(h, result$inter_arrival,
        result$type, result$mark, lambda_component0 = .1, N0 = args$N0,
        Nc0 = args$Nc0, infer = TRUE)
      for (name in c("lambda", "lambda_component", "rambda", "rambda_component"))
        stopifnot(identical(inferred[[name]], result[[name]]))
    }
  }
  comparisons <- comparisons + 1L
}
options(emhawkes.engine = "auto")
model <- fspec(kernels[[1]], laws[[1]])
set.seed(917); saved <- .Random.seed
if (!frpp_has_sim_api) {
  fails(fsim(model, backend = "state", seed = 18), "event-wise")
  stopifnot(identical(saved, .Random.seed))
} else {
  # A caller-supplied seed restores its kind and state, even after generator errors.
  a <- fsim(model, size = 10, seed = 18, backend = "state", lambda_component0 = 0)
  b <- simulate(model, size = 10, seed = 18, backend = "state", lambda_component0 = 0)
  stopifnot(identical(a, b), identical(saved, .Random.seed))
  bad <- model
  bad$residual$functions$rng <- function(n, shape, rate) {
    RNGkind("L'Ecuyer-CMRG"); runif(1); stop("generator failed")
  }
  kind <- RNGkind()
  fails(fsim(bad, seed = 18, backend = "state", lambda_component0 = 0), "generator failed")
  stopifnot(identical(saved, .Random.seed), identical(kind, RNGkind()))
  rm(".Random.seed", envir = .GlobalEnv)
  invisible(fsim(model, size = 1, backend = "state", lambda_component0 = 0))
  stopifnot(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
  fails(fsim(bad, seed = 18, backend = "state", lambda_component0 = 0), "generator failed")
  stopifnot(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE), identical(kind, RNGkind()))
  # No simulator delegation is allowed in the independent loop.
  local({
    trace("hsim", where = asNamespace("emhawkes"),
      tracer = quote(stop("legacy simulator invoked")), print = FALSE)
    on.exit(untrace("hsim", where = asNamespace("emhawkes")))
    invisible(fsim(model, size = 3, seed = 71, backend = "state", lambda_component0 = 0))
    invisible(fsim(model, size = 3, seed = 71, backend = "auto", lambda_component0 = 0))
    fails(fsim(model, size = 3, seed = 71, backend = "bridge", lambda_component0 = 0), "legacy simulator invoked")
  })
  # Initial/count expressions and additional arguments retain legacy forcing order.
  lazy <- function(backend) {
    events <- character()
    set.seed(415)
    value <- fsim(model, size = 8, seed = 73, backend = backend,
      lambda_component0 = { events <- c(events, "initial"); runif(1) },
      N0 = { events <- c(events, "N0"); runif(1) },
      Nc0 = { events <- c(events, "Nc0"); runif(1) })
    list(value = value, events = events, rng = .Random.seed)
  }
  stopifnot(identical(lazy("bridge"), lazy("state")),
    identical(lazy("bridge"), lazy("auto")))
  extra <- function(backend) {
    set.seed(416)
    path <- fsim(model, size = 8, seed = 73, backend = backend,
      lambda_component0 = 0, unused = { runif(1); FALSE })
    list(path = path, rng = .Random.seed)
  }
  stopifnot(identical(extra("bridge"), extra("auto")))
  # An invalid new post-event state fails before drawing the next residuals.
  calls <- 0L
  invalid <- fspec(methods::new("hspec", mu = 1, alpha = -.2, beta = 1),
    residual_law(rng = function(n) { calls <<- calls + 1L; rexp(n) }))
  failed <- function(backend) {
    calls <<- 0L; set.seed(417)
    error <- tryCatch(fsim(invalid, size = 3, backend = backend, lambda_component0 = 0), error = identity)
    stopifnot(inherits(error, "error"), calls == 1L)
    list(message = conditionMessage(error), calls = calls, rng = .Random.seed)
  }
  stopifnot(identical(failed("bridge"), failed("state")),
    identical(failed("bridge"), failed("auto")))
  # The residual generator receives the legacy dimension storage type too.
  typed_kernel <- kernels[[2]]; typed_kernel@dimens <- as.numeric(typed_kernel@dimens)
  seen <- character()
  typed <- fspec(typed_kernel, residual_law(rng = function(n) {
    seen <<- c(seen, typeof(n)); stats::rgamma(n, shape = 2, rate = 2)
  }))
  typed_run <- function(backend) {
    seen <<- character(); set.seed(419)
    path <- fsim(typed, size = 5, backend = backend, lambda_component0 = 0)
    list(path = path, seen = seen, rng = .Random.seed)
  }
  stopifnot(identical(typed_run("bridge"), typed_run("state")))
  # Existing deprecated input recycling retains its warning and output.
  recycled <- function(backend) {
    text <- character()
    set.seed(418)
    value <- withCallingHandlers(fsim(fspec(kernels[[3]], laws[[1]]), size = 3,
      backend = backend, lambda_component0 = c(.1, .2)),
      warning = function(w) { text <<- c(text, conditionMessage(w)); invokeRestart("muffleWarning") })
    list(value = value, warning = text, rng = .Random.seed)
  }
  stopifnot(identical(recycled("bridge"), recycled("state")))
  # Ties, zero intervals and finite total hazard retain the legacy event rule.
  deterministic <- residual_law(rng = function(n) rep(0, n))
  tied <- fspec(kernels[[2]], deterministic)
  a <- fsim(tied, size = 5, backend = "state", lambda_component0 = 0)
  b <- fsim(tied, size = 5, backend = "bridge", lambda_component0 = 0)
  stopifnot(identical(a, b), identical(a$type, c(0, rep(1, 4))))
  extinct <- fspec(methods::new("hspec", mu = 0, alpha = .1, beta = 1),
    residual_law(rng = function(n) rep(Inf, n)))
  warned <- function(backend) {
    text <- character()
    path <- withCallingHandlers(fsim(extinct, size = 5, backend = backend, lambda_component0 = 0),
      warning = function(w) { text <<- c(text, conditionMessage(w)); invokeRestart("muffleWarning") })
    list(path = path, warning = text)
  }
  stopifnot(identical(warned("state"), warned("bridge")))
  # Rejected support and invalid initial state fail before any RNG consumption.
  callback <- model; callback$kernel@mu <- function() matrix(.4)
  set.seed(712); saved <- .Random.seed
  fails(fsim(callback, backend = "state", seed = 18), "plain numeric")
  fails(fsim(model, backend = "state", lambda_component0 = -1), "Initial kernel")
  fails(fsim(model, backend = "state", unused = { stop("must not force") }), "additional")
  stopifnot(identical(saved, .Random.seed))
  set.seed(11); a <- fsim(callback, size = 5, lambda_component0 = 0); rng <- .Random.seed
  set.seed(11); b <- fsim(callback, size = 5, backend = "bridge", lambda_component0 = 0)
  stopifnot(identical(a, b), identical(rng, .Random.seed))
}
cat(comparisons, "paired simulation/RNG/following-draw comparisons; inference and backend contracts passed.\n")
