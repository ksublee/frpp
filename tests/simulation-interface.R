library(frpp)
model <- fspec(methods::new("hspec", mu = .4, alpha = .1, beta = 1),
  residual_law(c(shape = 2, rate = 2), stats::dgamma, stats::pgamma,
               stats::qgamma, stats::rgamma))
old <- structure(list(kernel = fkernel(model), residual = fresidual(model),
  refresh = model$refresh, version = 1L), class = "frpp_spec")
for (object in list(model, old)) for (size in c(1L, 20L)) {
  set.seed(18)
  seed <- .Random.seed
  a <- fsim(object, size = size, seed = 73, lambda_component0 = 0)
  b <- stats::simulate(object, size = size, seed = 73, lambda_component0 = 0)
  stopifnot(identical(a, b), identical(seed, .Random.seed))
  set.seed(73)
  a <- fsim(object, size = size, lambda_component0 = 0)
  seed <- .Random.seed
  next_draw <- stats::rnorm(3)
  set.seed(73)
  b <- stats::simulate(object, size = size, lambda_component0 = 0)
  stopifnot(identical(a, b), identical(seed, .Random.seed),
            identical(next_draw, stats::rnorm(3)))
}
fails <- function(expr, text) {
  err <- tryCatch(force(expr), error = identity)
  stopifnot(inherits(err, "error"), grepl(text, conditionMessage(err), fixed = TRUE))
}
set.seed(98)
seed <- .Random.seed
fails(fsim(NULL), "object must")
fails(fsim(model, size = 0, seed = 73), "size must")
fails(stats::simulate(model, nsim = 2, seed = 73), "nsim = 1")
stopifnot(identical(seed, .Random.seed))
bad <- model
bad$residual$functions$rng <- function(n, shape, rate) {
  RNGkind("L'Ecuyer-CMRG"); runif(1); stop("generator failed")
}
kind <- RNGkind()
fails(fsim(bad, size = 2, seed = 73, lambda_component0 = 0), "generator failed")
stopifnot(identical(seed, .Random.seed), identical(kind, RNGkind()))
rm(".Random.seed", envir = .GlobalEnv)
invisible(fsim(model, size = 1, seed = 73, lambda_component0 = 0))
stopifnot(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
fails(fsim(bad, size = 2, seed = 73, lambda_component0 = 0), "generator failed")
stopifnot(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE), identical(kind, RNGkind()))

# Count validation stays at the FRPP observation boundary after kernel replay.
dt <- c(0, .2, 0, .1)
types <- c(0, 1, 1, 1)
marks <- c(0, 2, 1, .5)
N <- matrix(3:6, ncol = 1)
Nc <- matrix(c(5, 7, 8, 8.5), ncol = 1)
args <- list(object = model, inter_arrival = dt, type = types, mark = marks,
  N = N, Nc = Nc, N0 = 3, Nc0 = 5, lambda_component0 = 0)
for (backend in c("bridge", "state")) {
  value <- do.call(stats::logLik, c(args, list(backend = backend)))
  full <- do.call(stats::logLik, c(args, list(backend = backend, infer = TRUE)))
  stopifnot(identical(unname(value), unname(full$loglikelihood)))
  invalid <- args; invalid$N[3, 1] <- 99
  fails(do.call(stats::logLik, c(invalid, list(backend = backend))), "N is inconsistent")
  invalid <- args; invalid$Nc0 <- 3
  fails(do.call(stats::logLik, c(invalid, list(backend = backend))), "Nc is inconsistent")
  invalid <- args; invalid$N0 <- numeric(0)
  fails(do.call(stats::logLik, c(invalid, list(backend = backend))), "N0 must")
}
