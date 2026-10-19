library(frpp)
source(file.path("helpers", "kernel-api.R"))
frpp_numeric_backend <- if (frpp_has_kernel_api) "state" else "bridge"
close <- function(x, y) stopifnot(isTRUE(all.equal(as.numeric(x), as.numeric(y),
                                                  tolerance = 1e-10)))

# Reconstruct the numeric driver directly from all past events, then integrate
# it by quadrature. This does not reuse the prototype's recursion or integral.
kernel <- methods::new("hspec", mu = c(.3, .4),
  alpha = matrix(c(.05, .07, .03, .08, .02, .04, .06, .03), 2, 4),
  beta = matrix(rep(c(1.2, 1.7), 4), 2), eta = matrix(.02, 2, 4),
  type_col_map = list(c(1L, 3L), c(2L, 4L)))
dt <- c(0, .2, .4, .1, .7, .3)
type <- c(0, 1, 2, 1, 1, 2)
mark <- c(0, 1, 2, 1.5, 1, 3)
initial <- matrix(seq_len(8) / 100, 2, 4)
arrival <- cumsum(dt)
driver <- function(time, past) {
  value <- initial * exp(-kernel@beta * time)
  for (j in past) {
    jump <- matrix(0, 2, 4)
    columns <- kernel@type_col_map[[type[j]]]
    jump[, columns] <- kernel@alpha[, columns] + kernel@eta[, columns] * (mark[j] - 1)
    value <- value + jump * exp(-kernel@beta * (time - arrival[j]))
  }
  as.vector(kernel@mu) + rowSums(value)
}
phi <- lambda <- matrix(0, length(dt), 2)
lambda[1, ] <- driver(0, integer())
expected <- 0
for (n in 2:length(dt)) {
  past <- if (n == 2L) integer() else seq.int(2L, n - 1L)
  lambda[n, ] <- driver(arrival[n], past)
  for (k in 1:2) {
    phi[n, k] <- stats::integrate(function(t) vapply(t,
      function(s) driver(s, past)[k], numeric(1)),
      arrival[n - 1L], arrival[n], rel.tol = 1e-12)$value
    expected <- expected + if (k == type[n]) {
      log(lambda[n, k]) + stats::dgamma(phi[n, k], 2, 2, log = TRUE)
    } else stats::pgamma(phi[n, k], 2, 2, lower.tail = FALSE, log.p = TRUE)
  }
}
model <- frpp_spec(kernel, residual_law(c(shape = 2, rate = 2),
  stats::dgamma, stats::pgamma, stats::qgamma, stats::rgamma))
maximum_error <- 0
for (engine in if (frpp_has_kernel_api) c("r", "cpp") else "r") {
  options(emhawkes.engine = engine)
  result <- stats::logLik(model, dt, type, mark, lambda_component0 = initial,
    infer = TRUE, backend = frpp_numeric_backend)
  # The released bridge does not expose per-interval integrated drivers.
  # Its intensities and likelihood still have independent all-past/quadrature references.
  if (frpp_has_kernel_api) {
    close(result$integrated_driver, phi)
    maximum_error <- max(maximum_error, abs(result$integrated_driver - phi))
  }
  close(result$lambda, lambda)
  close(result$loglikelihood, expected)
}
options(emhawkes.engine = "auto")
bridge <- stats::logLik(model, dt, type, mark, lambda_component0 = initial, infer = TRUE)
if (frpp_has_kernel_api) {
  state <- stats::logLik(model, dt, type, mark, lambda_component0 = initial,
    infer = TRUE, backend = "state")
  for (name in names(bridge)) close(state[[name]], bridge[[name]])
}

# Exp(1) recovers the Hawkes likelihood under the same explicit initial state.
exponential <- frpp_spec(kernel, residual_law(c(rate = 1), stats::dexp,
  stats::pexp, stats::qexp, stats::rexp))
close(stats::logLik(exponential, dt, type, mark, lambda_component0 = initial,
  backend = frpp_numeric_backend), emhawkes::logLik(kernel, dt, type, mark, lambda_component0 = initial))

# Refresh semantics: even a type that did not fire receives a fresh candidate
# at the next global event. With no excitation, targets / mu give exact waits.
draws <- list()
rng <- function(n) {
  value <- stats::rgamma(n, shape = 2, rate = 2)
  draws[[length(draws) + 1L]] <<- value
  value
}
constant <- methods::new("hspec", mu = c(.4, .6), alpha = matrix(0, 2, 2),
  beta = matrix(1, 2, 2))
refresh <- frpp_spec(constant, residual_law(rng = rng))
set.seed(211)
simulated <- stats::simulate(refresh, size = 9, lambda_component0 = 0)
stopifnot(length(draws) == 8L, all(lengths(draws) == 2L))
for (j in seq_along(draws)) {
  candidates <- draws[[j]] / c(.4, .6)
  close(simulated$inter_arrival[j + 1L], min(candidates))
  stopifnot(simulated$type[j + 1L] == which.min(candidates))
}

# Retaining each dimension's renewal clock produces a different likelihood for
# Gamma residuals. It must not be introduced as a seemingly harmless refactor.
dt2 <- c(0, .2, .3, .5, .1)
type2 <- c(0, 1, 2, 1, 2)
constant_model <- frpp_spec(constant, model$residual)
global <- as.numeric(stats::logLik(constant_model, dt2, type2,
  lambda_component0 = 0, backend = frpp_numeric_backend))
renewal <- 0
time <- cumsum(dt2)
for (k in 1:2) {
  times <- time[type2 == k]
  rate <- constant@mu[k, 1]
  renewal <- renewal + length(times) * log(rate) +
    sum(stats::dgamma(rate * diff(c(0, times)), 2, 2, log = TRUE)) +
    stats::pgamma(rate * (tail(time, 1) - tail(times, 1)), 2, 2,
      lower.tail = FALSE, log.p = TRUE)
}
stopifnot(abs(global - renewal) > 1e-3)

# A free residual scale can be confounded with kernel scale. No estimator is
# added in this prototype before a normalization/identification rule is agreed.
scaled <- vapply(c(.2, .4, .8), function(mu) {
  h <- methods::new("hspec", mu = mu, alpha = 0, beta = 1)
  m <- frpp_spec(h, residual_law(c(rate = .8 / mu), stats::dexp, stats::pexp))
  as.numeric(stats::logLik(m, c(0, .2, .3, .5), lambda_component0 = 0,
    backend = frpp_numeric_backend))
}, numeric(1))
close(scaled, rep(scaled[1], 3))

# Empty type histories, simultaneous events and log-scale tails remain defined.
for (law in list(exponential$residual, model$residual)) {
  for (observed_type in list(c(0, 1, 1), c(0, 2, 1))) {
    for (observed_dt in list(c(0, .1, 1000), c(0, 0, .2))) {
      m <- frpp_spec(constant, law)
      a <- stats::logLik(m, observed_dt, observed_type, lambda_component0 = 0, infer = TRUE)
      b <- stats::logLik(m, observed_dt, observed_type, lambda_component0 = 0,
        infer = TRUE, backend = frpp_numeric_backend)
      for (name in names(a)) close(a[[name]], b[[name]])
    }
  }
}
if (frpp_has_kernel_api) {
  cat("Independent quadrature maximum integral error:", format(maximum_error, digits = 17), "\n")
} else {
  cat("Per-interval integral output unavailable in bridge; independent intensity and likelihood references passed.\n")
}
cat("Global refresh vs renewal-clock log-likelihood difference:", global - renewal, "\n")

if (!frpp_has_kernel_api) {
  set.seed(712)
  saved_rng <- .Random.seed
  unavailable <- tryCatch(stats::logLik(model, dt, type, mark,
    lambda_component0 = initial, backend = "state"), error = identity)
  stopifnot(inherits(unavailable, "error"),
    grepl("development kernel API", conditionMessage(unavailable), fixed = TRUE),
    identical(saved_rng, .Random.seed))
}
