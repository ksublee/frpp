library(frpp)
fails <- function(expr, text) {
  error <- tryCatch(force(expr), error = identity)
  stopifnot(inherits(error, "error"), grepl(text, conditionMessage(error), fixed = TRUE))
}
k <- methods::new("hspec", mu = .4, alpha = .1, beta = 1.2)
law <- residual_law(c(shape = 2, rate = 2), dgamma, pgamma, qgamma, rgamma)
m <- fspec(k, law)
stopifnot(isS4(m), methods::is(m, "fspec"), !methods::is(m, "hspec"),
  methods::validObject(m), identical(fkernel(m), k), identical(fresidual(m), law),
  identical(m$kernel, k), identical(m$residual, law), is.null(m$unknown),
  identical(as_fspec(m), m), identical(methods::new("fspec", kernel = k, residual = law), m))
fails(as_fspec(structure(1, class = "fspec")), "object must")
# Typed slots and semantic validity cover constructor, replacement and direct edits.
frozen <- serialize(m, NULL)
fails(fkernel(m) <- 1, "hspec")
wrong <- law; wrong$binding <- "unknown"
fails(fresidual(m) <- wrong, "binding")
fails(m$unknown <- 1, "Unknown fspec field")
stopifnot(identical(serialize(m, NULL), frozen))
for (kind in c("beta", "columns", "mapping", "residual", "refresh", "version")) {
 bad <- m
 if (kind == "beta") bad@kernel@beta[1, 1] <- -1
 if (kind == "columns") bad@kernel@alpha <- matrix(.1, 1, 2)
 if (kind == "mapping") bad@kernel@type_col_map <- list(2L)
 if (kind == "residual") bad@residual$functions$rng <- 1
 if (kind == "refresh") bad@refresh <- "dimension_only"
 if (kind == "version") bad@version <- 99L
 stopifnot(!isTRUE(methods::validObject(bad, test = TRUE)))
 set.seed(71); rng <- .Random.seed
 fails(fsim(bad, size = 2, seed = 19), "invalid class")
 stopifnot(identical(rng, .Random.seed))
}
replacement <- m
fkernel(replacement) <- methods::new("hspec", mu = .5, alpha = .2, beta = 1)
fresidual(replacement) <- residual_law(c(rate = 1), dexp, pexp, qexp, rexp)
stopifnot(fkernel(replacement)@mu == .5, identical(fkernel(m), k))
# Both serialized S3 generations upgrade without changing the scientific model.
for (classes in list("frpp_spec", c("fspec", "frpp_spec"))) {
 old <- structure(list(kernel = k, residual = law, refresh = m$refresh, version = 1L), class = classes)
 old <- unserialize(serialize(old, NULL))
 stopifnot(!isS4(old), identical(as_fspec(old), m),
   identical(fresidual(as_fspec(old)), fresidual(m)))
 set.seed(73); a <- fsim(old, size = 12, lambda_component0 = .1); rng <- .Random.seed
 set.seed(73); b <- fsim(m, size = 12, lambda_component0 = .1)
 stopifnot(identical(a, b), identical(rng, .Random.seed))
 old$version <- 99L
 fails(as_fspec(old), "Unsupported legacy")
}
# Introspection must not evaluate any callback or default expression.
env <- new.env(parent = baseenv()); env$calls <- 0L
rng <- eval(quote(function(n, param = stop("default evaluated")) {
 calls <- calls + 1L; stop("generator executed")
}), env)
functional <- methods::new("hspec", mu = function(n) stop("baseline executed"),
 alpha = .1, beta = 1.2, impact = function(...) stop("impact executed"))
bound <- functional; bound@rresidual <- rng
converted <- as_fspec(bound)
stopifnot(methods::validObject(converted), identical(as_legacy_hspec(converted), bound),
 identical(fresidual(converted)$functions$rng, rng), env$calls == 0L)
rm(".Random.seed", envir = .GlobalEnv)
invisible(fspec(functional, residual_law(rng = rng)))
invisible(as_fspec(bound))
stopifnot(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE), env$calls == 0L)
# Explicit stats S3 entry points and formal S4 generics share the same methods.
path <- fsim(m, size = 12, seed = 17, lambda_component0 = .1)
args <- list(object = m, inter_arrival = path$inter_arrival, type = path$type,
 lambda_component0 = .1)
expected <- do.call(stats::logLik, args)
stopifnot(identical(do.call(emhawkes::logLik, args), expected),
 identical(do.call(logLik, args), expected),
 identical(stats::simulate(m, size = 12, seed = 17, lambda_component0 = .1), path),
 identical(simulate(m, size = 12, seed = 17, lambda_component0 = .1), path))
invisible(capture.output(methods::show(m)))
invisible(capture.output(print(m)))
cat("S4 composition, validity, accessors, legacy conversion and dual dispatch: OK\n")
