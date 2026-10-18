library(frpp)
check_fit_session <- function() {
  data <- tempfile(fileext = ".rds"); libs <- tempfile(fileext = ".rds"); script <- tempfile(fileext = ".R")
  on.exit(unlink(c(data, libs, script)), add = TRUE)
  model <- fspec(methods::new("hspec", mu = 1, alpha = 0, beta = 1),
    residual_law(c(rate = 1), stats::dexp, stats::pexp, stats::qexp, stats::rexp))
  dt <- c(0, .2, .3, .5)
  fit <- ffit(model, dt, optimizer = "stats", lambda_component0 = .1,
    fixed = c(FALSE, TRUE, TRUE), N0 = 2, Nc0 = 3)
  saveRDS(list(fit = fit, dt = dt), data); saveRDS(.libPaths(), libs)
  writeLines(c("args <- commandArgs(TRUE)", ".libPaths(readRDS(args[1]))", "library(frpp)",
    "stopifnot(!('package:methods' %in% search()), !('package:emhawkes' %in% search()))",
    "stopifnot(!('maxLik' %in% loadedNamespaces()))", "x <- readRDS(args[2])",
    "stopifnot(inherits(fitted(x$fit), 'fspec'))",
    "value <- do.call(logLik, c(list(object = fitted(x$fit), inter_arrival = x$dt), finitial_state(x$fit)))",
    "stopifnot(identical(as.numeric(value), as.numeric(logLik(x$fit))))",
    "stopifnot(identical(coef(x$fit), x$fit$coefficients), nobs(x$fit) == 3L)",
    "invisible(summary(x$fit)); invisible(vcov(x$fit))",
    "stopifnot(!('maxLik' %in% loadedNamespaces()))",
    "cat('Saved stats ffit restored in minimal frpp-only session: OK\\n')"), script)
  exe <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  output <- suppressWarnings(system2(exe, c("--vanilla",
    "--default-packages=datasets,utils,grDevices,graphics,stats", shQuote(script), shQuote(libs), shQuote(data)),
    stdout = TRUE, stderr = TRUE))
  if (!is.null(attr(output, "status"))) stop(paste(output, collapse = "\n"))
}
check_fit_session()
