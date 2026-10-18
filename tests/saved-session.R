library(frpp)
check_saved_session <- function() {
  libraries <- tempfile(fileext = ".rds")
  saved <- tempfile(fileext = ".rds")
  script <- tempfile(fileext = ".R")
  on.exit(unlink(c(libraries, saved, script)), add = TRUE)
  h <- methods::new("hspec", mu = .4, alpha = .15, beta = 1.2)
  m <- frpp_spec(h, residual_law(c(shape = 2, rate = 2), stats::dgamma,
    stats::pgamma, stats::qgamma, stats::rgamma))
  converted <- as_frpp(as_legacy_hspec(m))
  x <- stats::simulate(m, size = 10, seed = 217, lambda_component0 = .1)
  expected <- stats::logLik(m, x$inter_arrival, x$type, lambda_component0 = .1)
  old <- structure(list(kernel = fkernel(m), residual = fresidual(m),
    refresh = m$refresh, version = 1L), class = "frpp_spec")
  old_named <- old; class(old_named) <- c("fspec", "frpp_spec")
  saveRDS(list(raw = m, converted = converted, old = old, old_named = old_named,
    path = x, value = expected), saved)
  saveRDS(.libPaths(), libraries)
  writeLines(c(
    "args <- commandArgs(TRUE)",
    ".libPaths(readRDS(args[1]))",
    "library(frpp)",
    "stopifnot(!('package:methods' %in% search()), !('package:emhawkes' %in% search()))",
    "x <- readRDS(args[2])",
    "for (m in list(x$raw, x$converted, x$old, x$old_named)) {",
    "  stopifnot(isS4(as_fspec(m)), methods::validObject(as_fspec(m)))",
    "  for (backend in c('bridge', 'state')) {",
    "    value <- stats::logLik(m, x$path$inter_arrival, x$path$type,",
    "      lambda_component0 = .1, backend = backend)",
    "    stopifnot(isTRUE(all.equal(value, x$value, tolerance = 1e-10)))",
    "  }",
    "}",
    "h <- methods::new('hspec', mu = .2, alpha = 0, beta = 1)",
    "stopifnot(length(h@dimens) == 1L, h@dimens == 1)",
    "cat('Saved composition and legacy conversion in a minimal session: OK\\n')"
  ), script)
  executable <- file.path(R.home("bin"),
    if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  output <- suppressWarnings(system2(executable,
    c("--vanilla", "--default-packages=datasets,utils,grDevices,graphics,stats",
      shQuote(script), shQuote(libraries), shQuote(saved)),
    stdout = TRUE, stderr = TRUE))
  if (!is.null(attr(output, "status")))
    stop(paste(c("Fresh-session composition restoration failed:", output), collapse = "\n"), call. = FALSE)
}
check_saved_session()
