library(frpp)
stopifnot(!("package:emhawkes" %in% search()))
close <- function(x, y) stopifnot(isTRUE(all.equal(x, y, tolerance = 1e-10)))
fails <- function(expr, text) {
  error <- tryCatch(force(expr), error = identity)
  stopifnot(inherits(error, "error"), grepl(text, conditionMessage(error), fixed = TRUE))
}

kernels <- list(
  methods::new("hspec", mu = .4, alpha = .15, beta = 1.2),
  methods::new("hspec", mu = c(.3, .4), alpha = matrix(.06, 2, 2),
    beta = matrix(c(1.2, 1.7, 1.2, 1.7), 2)),
  methods::new("hspec", mu = c(.3, .4), alpha = matrix(.025, 2, 4),
    beta = matrix(rep(c(1.2, 1.7), 4), 2),
    type_col_map = list(c(1L, 3L), c(2L, 4L))))
laws <- list(
  residual_law(c(shape = 2, rate = 2), stats::dgamma, stats::pgamma,
    stats::qgamma, stats::rgamma, "Gamma"),
  residual_law(c(rate = 1), stats::dexp, stats::pexp,
    stats::qexp, stats::rexp, "Exponential"),
  residual_law(c(a = .7, ell = 1.6), emhawkes::dtzexp, emhawkes::ptzexp,
    emhawkes::qtzexp, emhawkes::rtzexp, "tzexp"))
comparisons <- 0L
maximum_likelihood_difference <- 0
for (kernel in kernels) for (law in laws) {
  model <- frpp_spec(kernel, law)
  stopifnot(!methods::is(model, "hspec"), methods::is(model, "fspec"),
    identical(model$kernel, kernel))
  legacy <- as_legacy_hspec(model)
  converted <- as_frpp(legacy)
  stopifnot(identical(as_legacy_hspec(converted), legacy),
    identical(converted$kernel, kernel))
  for (seed in c(1, 73)) for (size in c(1, 35)) {
    initial <- matrix(.08, kernel@dimens, ncol(kernel@beta))
    set.seed(seed)
    expected <- emhawkes::hsim(legacy, size = size,
      lambda_component0 = initial, N0 = 3, Nc0 = 7)
    expected_seed <- .Random.seed
    set.seed(seed)
    actual <- stats::simulate(model, size = size,
      lambda_component0 = initial, N0 = 3, Nc0 = 7)
    stopifnot(identical(.Random.seed, expected_seed))
    fields <- setdiff(names(expected), "hspec")
    stopifnot(identical(unclass(actual)[fields], unclass(expected)[fields]))
    args <- list(inter_arrival = actual$inter_arrival, type = actual$type,
      mark = actual$mark, lambda_component0 = initial, N0 = 3, Nc0 = 7)
    old <- suppressMessages(do.call(emhawkes::logLik, c(list(object = legacy), args)))
    bridge <- suppressMessages(do.call(stats::logLik, c(list(object = model), args)))
    state <- suppressMessages(do.call(stats::logLik,
      c(list(object = model, backend = "state"), args)))
    stopifnot(identical(old, bridge))
    close(state, old)
    maximum_likelihood_difference <- max(maximum_likelihood_difference, abs(state - old))
    comparisons <- comparisons + 1L
  }
}
cat(comparisons, "simulation/RNG and likelihood comparisons passed.\n")
cat("Maximum state-backend likelihood difference:",
    format(maximum_likelihood_difference, digits = 17), "\n")

# The default FRPP state must remain zero even when the Hawkes kernel is stable.
model <- frpp_spec(kernels[[1]], laws[[1]])
set.seed(97)
automatic <- suppressMessages(stats::simulate(model, size = 8))
automatic_seed <- .Random.seed
set.seed(97)
zero <- stats::simulate(model, size = 8, lambda_component0 = 0)
stopifnot(identical(unclass(automatic), unclass(zero)),
  identical(.Random.seed, automatic_seed))
for (backend in c("bridge", "state")) {
  auto <- suppressMessages(stats::logLik(model, c(0, .2, .3), backend = backend))
  explicit <- stats::logLik(model, c(0, .2, .3), lambda_component0 = 0, backend = backend)
  stopifnot(identical(auto, explicit))
  fails(stats::logLik(model, c(0, .2), lambda_component0 = "mean", backend = backend), "mean")
}

# Capability failures and zero-row work must not consume random draws.
no_rng <- frpp_spec(kernels[[1]], residual_law(c(rate = 1), stats::dexp))
before <- .Random.seed
fails(stats::simulate(no_rng, size = 3), "rng")
stopifnot(identical(before, .Random.seed))
fails(residual_law(c(rate = NA_real_), stats::dexp), "finite")
fails(residual_law(c(rtae = 1), stats::dexp), "rtae")
fails(frpp_spec(as_legacy_hspec(model), laws[[1]]), "as_frpp")

# The seed argument follows the standard simulate convention and restores even
# the absence of a seed. Constructors never initialize the RNG themselves.
saved <- .Random.seed
invisible(stats::simulate(model, seed = 18, size = 5, lambda_component0 = 0))
stopifnot(identical(saved, .Random.seed))
rm(".Random.seed", envir = .GlobalEnv)
invisible(frpp_spec(kernels[[1]], laws[[1]]))
stopifnot(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
invisible(stats::simulate(model, seed = 18, size = 1, lambda_component0 = 0))
stopifnot(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))

# Conversion preserves legacy generator-only laws and never evaluates callbacks.
generator_only <- kernels[[1]]
generator_only@rresidual <- function(n) stop("must not run during conversion")
converted <- as_frpp(generator_only)
stopifnot(identical(as_legacy_hspec(converted), generator_only))
fails(stats::logLik(converted, c(0, 1)), "density")

# Callback defaults, marks, counts, histories and external environments survive
# the bridge. They remain unsupported by the independent numeric state backend.
callback <- methods::new("hspec", mu = function(N, n, scale = .01) {
  if (n <= 1L) matrix(.4) else matrix(.4 + scale * N[n - 1L, 1])
}, alpha = .1, beta = 1,
  impact = function(lambda_component_n, extra = .02) lambda_component_n * extra,
  rmark = function(mean = .2) 1 + stats::rexp(1, 1 / mean))
callback_model <- frpp_spec(callback, laws[[1]])
set.seed(111)
reference <- emhawkes::hsim(as_legacy_hspec(callback_model), size = 15,
  lambda_component0 = 0, N0 = 2, Nc0 = 4)
saved <- .Random.seed
set.seed(111)
actual <- stats::simulate(callback_model, size = 15,
  lambda_component0 = 0, N0 = 2, Nc0 = 4)
fields <- setdiff(names(reference), "hspec")
stopifnot(identical(unclass(actual)[fields], unclass(reference)[fields]),
  identical(saved, .Random.seed))
args <- list(inter_arrival = actual$inter_arrival, type = actual$type,
  mark = actual$mark, lambda_component0 = 0, N0 = 2, Nc0 = 4, infer = TRUE)
stopifnot(identical(do.call(stats::logLik, c(list(object = callback_model), args)),
  do.call(emhawkes::logLik, c(list(object = as_legacy_hspec(callback_model)), args))))
fails(do.call(stats::logLik, c(list(object = callback_model, backend = "state"), args)), "numeric")

# Methods also dispatch when the user has emhawkes attached.
library(emhawkes)
close(logLik(model, c(0, .2, .3), lambda_component0 = 0),
      stats::logLik(model, c(0, .2, .3), lambda_component0 = 0))
