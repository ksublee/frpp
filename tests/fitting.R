library(frpp)
source(file.path("helpers", "kernel-api.R"))
stopifnot(!("package:emhawkes" %in% search()))
close <- function(x, y, tol = 1e-7) {
  if (!isTRUE(all.equal(x, y, tolerance = tol))) { print(x, digits = 17); print(y, digits = 17); stop("Analytic comparison failed") }
}
fails <- function(expr, text) {
  x <- tryCatch(force(expr), error = identity)
  stopifnot(inherits(x, "error"), grepl(text, conditionMessage(x), fixed = TRUE))
}
law <- residual_law(c(rate = 1), stats::dexp, stats::pexp, stats::qexp, stats::rexp)
model <- fspec(methods::new("hspec", mu = 1, alpha = 0, beta = 1), law)
stopifnot(inherits(model, "fspec"), !methods::is(model, "hspec"))
legacy_class <- structure(list(kernel = fkernel(model), residual = fresidual(model),
  refresh = model$refresh, version = 1L), class = "frpp_spec")
stopifnot(inherits(as_fspec(legacy_class), "fspec"))
stopifnot(identical(logLik(model, c(0, .2)), logLik(legacy_class, c(0, .2))))
dt <- c(0, .2, .3, .5)
providers <- c("stats", if (requireNamespace("maxLik", quietly = TRUE)) "maxLik")
precise_fit <- function(..., optimizer = "stats") ffit(..., optimizer = optimizer,
  control = if (optimizer == "stats") list(reltol = 1e-12, maxit = 1000, ndeps = 1e-4) else
    list(reltol = 1e-12, iterlim = 1000))
saved_model <- serialize(model, NULL)
for (provider in providers) {
  fit <- precise_fit(model, dt, optimizer = provider, lambda_component0 = 0,
    fixed = c(FALSE, TRUE, TRUE),
    hess = function(param) diag(c(-3 / param["mu1"]^2, 0, 0)))
  stopifnot(fdiagnostics(fit)$optimizer_success, inherits(fitted(fit), "fspec"),
    identical(fitted_fspec(fit)$residual$param, law$param), nobs(fit) == 3L,
    attr(logLik(fit), "df") == 1L, !"rate" %in% names(coef(fit)))
  close(unname(coef(fit)["mu1"]), 3)
  close(unname(vcov(fit)[1, 1]), 3, 1e-4)
  close(as.numeric(logLik(fit)), 3 * log(3) - 3)
  replay <- do.call(logLik, c(list(object = fitted_fspec(fit), inter_arrival = dt), finitial_state(fit)))
  close(as.numeric(logLik(fit)), as.numeric(replay), 1e-12)
  # Estimate a residual parameter with the kernel fixed: same analytic MLE.
  joint <- precise_fit(model, dt, optimizer = provider, estimate_residual = TRUE,
    fixed = c(TRUE, TRUE, TRUE, FALSE), lambda_component0 = 0)
  close(unname(coef(joint)["rate"]), 3)
  close(unname(fitted(joint)$residual$param["rate"]), 3)
  stopifnot(attr(logLik(joint), "df") == 1L)
  set.seed(91)
  actual <- simulate(fitted(joint), size = 8, lambda_component0 = 0)
  set.seed(91)
  expected <- simulate(fspec(model$kernel, residual_law(c(rate = unname(coef(joint)["rate"])),
    stats::dexp, stats::pexp, stats::qexp, stats::rexp)), size = 8, lambda_component0 = 0)
  stopifnot(identical(actual$inter_arrival, expected$inter_arrival))
}
stopifnot(identical(serialize(model, NULL), saved_model))
if (frpp_has_kernel_api) {
  state <- precise_fit(model, dt, optimizer = "stats", backend = "state",
    fixed = c(FALSE, TRUE, TRUE), lambda_component0 = 0)
  close(unname(coef(state)["mu1"]), 3)
}
converted <- as_fspec(as_legacy_hspec(model))
converted_fit <- precise_fit(converted, dt, optimizer = "stats", lambda_component0 = 0,
  fixed = c(FALSE, TRUE, TRUE))
close(unname(coef(converted_fit)["mu1"]), 3)
automatic <- suppressMessages(ffit(model, dt, optimizer = "stats", fixed = c(FALSE, TRUE, TRUE)))
stopifnot(all(finitial_state(automatic)$lambda_component0 == 0))
fails(ffit(model, dt, lambda_component0 = "mean"), "Hawkes mean")
fails(ffit(model, dt, constraint = list(), constraints = list()), "only one")
fails(ffit(model, dt, gradient = "analytic"), "unsupported argument")
fails(ffit(model, dt, optimizer = "stats", valid = function(param) 1), "logical")
fails(ffit(model, dt, optimizer = "stats", valid = function(param) stop("validator exploded")), "validator exploded")
validated <- precise_fit(model, dt, optimizer = "stats", fixed = c(FALSE, TRUE, TRUE),
  lambda_component0 = 0, valid = function(param) { stopifnot(param["rate"] == 1); TRUE })
close(unname(coef(validated)["mu1"]), 3)
# Flat beta cannot acquire a fabricated Wald covariance.
unidentified <- suppressWarnings(ffit(model, dt, optimizer = "stats", lambda_component0 = 0,
  fixed = c(FALSE, TRUE, FALSE)))
stopifnot(fdiagnostics(unidentified)$wald == "unavailable",
  all(is.na(vcov(unidentified)[c(1, 3), c(1, 3)])))
bound <- suppressWarnings(ffit(model, dt, optimizer = "stats", method = "L-BFGS-B",
  fixed = c(FALSE, TRUE, TRUE), lambda_component0 = 0,
  lower = c(0, -Inf, -Inf), upper = c(2, Inf, Inf)))
close(unname(coef(bound)["mu1"]), 2)
stopifnot(fdiagnostics(bound)$wald == "unavailable",
  fdiagnostics(bound)$reason == "constraint_boundary")
ineq <- list(ineqA = matrix(c(1, 0, 0), 1), ineqB = 0)
a <- ffit(model, dt, optimizer = "stats", fixed = c(FALSE, TRUE, TRUE),
  lambda_component0 = 0, constraints = ineq)
b <- ffit(model, dt, optimizer = "stats", fixed = c(FALSE, TRUE, TRUE),
  lambda_component0 = 0, constraint = ineq)
stopifnot(identical(coef(a), coef(b)))
# Fixed raw laws remain fixed even when their names match a kernel parameter.
collision <- fspec(model$kernel, residual_law(c(mu1 = 1),
  density = function(x, mu1, log = FALSE) dexp(x, rate = mu1, log = log),
  rng = function(n, mu1) rexp(n, rate = mu1)))
cf <- ffit(collision, dt, optimizer = "stats", fixed = c(FALSE, TRUE, TRUE), lambda_component0 = 0)
stopifnot(identical(fitted(cf)$residual$param, c(mu1 = 1)))
mapped <- fspec(methods::new("hspec", mu = c(.2, .2), alpha = matrix(.04, 2, 4),
  beta = matrix(1.2, 2, 4), type_col_map = list(c(1L, 3L), c(2L, 4L))), law)
path <- simulate(mapped, size = 30, seed = 17, lambda_component0 = .15, N0 = c(1, 2), Nc0 = c(3, 4))
for (reduced in c(TRUE, FALSE)) {
  mf <- ffit(mapped, path$inter_arrival, path$type, optimizer = "stats", reduced = reduced,
    fixed = rep(TRUE, if (reduced) 3 else 18), lambda_component0 = .15,
    N = path$N, Nc = path$Nc)
  stopifnot(attr(logLik(mf), "df") == 0, identical(fitted(mf)$kernel@type_col_map, mapped$kernel@type_col_map),
    identical(finitial_state(mf)$N0, c(1, 2)), identical(finitial_state(mf)$Nc0, c(3, 4)))
  close(as.numeric(logLik(mf)), as.numeric(logLik(mapped, path$inter_arrival, path$type,
    lambda_component0 = .15, N0 = c(1, 2), Nc0 = c(3, 4))), 1e-12)
}
impact <- function(param = c(k = .02), n, Nc) matrix(param["k"] * (Nc[n, 1] %% 2))
callbacks <- fspec(methods::new("hspec", mu = .4, alpha = .1, beta = 1.2,
  impact = impact, rmark = function() rpois(1, 1) + 1), law)
path <- simulate(callbacks, size = 12, seed = 27, lambda_component0 = .1, N0 = 2, Nc0 = 3)
cb <- ffit(callbacks, path$inter_arrival, path$type, path$mark, optimizer = "stats",
  fixed = rep(TRUE, 4), lambda_component0 = .1, N0 = 2, Nc0 = 3)
close(as.numeric(logLik(cb)), as.numeric(logLik(callbacks, path$inter_arrival, path$type,
  path$mark, lambda_component0 = .1, N0 = 2, Nc0 = 3)), 1e-12)
fails(ffit(callbacks, path$inter_arrival, path$type, path$mark, optimizer = "stats",
  backend = "state", lambda_component0 = .1), "numeric kernels")
cat("FRPP analytic fits, optimizer providers, restoration, aliases and unavailable Wald covariance: OK\n")
