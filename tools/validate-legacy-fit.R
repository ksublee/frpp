# Run from the frpp repository after installing frpp and emhawkes:
# Rscript tools/validate-legacy-fit.R <output-directory> [repetitions]
# Fits are sequential; alternate their order and exclude package loading,
# simulation and likelihood warm-up from the elapsed fitting times.
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || length(args) > 2L) stop("Supply an output directory and optional repetition count.")
out <- args[1L]
repetitions <- if (length(args) == 2L) as.integer(args[2L]) else 3L
stopifnot(length(repetitions) == 1L, !is.na(repetitions), repetitions >= 1L)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages(library(frpp))
stopifnot(requireNamespace("maxLik", quietly = TRUE))
cat("R:", R.version.string, "\nemhawkes:", as.character(packageVersion("emhawkes")),
  "\nfrpp:", as.character(packageVersion("frpp")), "\nmaxLik:", as.character(packageVersion("maxLik")), "\n")
cat("Paths:", find.package("emhawkes"), find.package("frpp"), "\n")
control <- list(iterlim = 1000L, reltol = 1e-8)
constraints <- list(ineqA = rbind(diag(6), c(0, 0, 0, 0, -1, 0)),
  ineqB = c(rep(0, 6), 3))
kernel <- methods::new("hspec", mu = c(0.3, 0.3),
  alpha = matrix(c(0.5, 0.8, 0.8, 0.5), 2, byrow = TRUE),
  beta = matrix(1.5, 2, 2))
law <- residual_law(c(a = 0.7, ell = 1.6), emhawkes::dtzexp,
  emhawkes::ptzexp, emhawkes::qtzexp, emhawkes::rtzexp, name = "tzexp")
model <- fspec(kernel, law)
legacy <- emhawkes::set_residual(kernel, param = c(a = 0.7, ell = 1.6),
  dresidual = emhawkes::dtzexp, presidual = emhawkes::ptzexp,
  qresidual = emhawkes::qtzexp, rresidual = emhawkes::rtzexp)
set.seed(1107)
path <- fsim(model, size = 10000, lambda_component0 = 0)
path_seed <- .Random.seed
set.seed(1107)
old_path <- emhawkes::hsim(legacy, size = 10000, lambda_component0 = 0)
fields <- setdiff(names(old_path), "hspec")
stopifnot(identical(unclass(path)[fields], unclass(old_path)[fields]),
  identical(path_seed, .Random.seed))
saveRDS(list(model = model, legacy = legacy, path = path, seed = path_seed,
  controls = control, constraints = constraints), file.path(out, "inputs.rds"))
validator_calls <- 0L
valid <- function(param) {
  validator_calls <<- validator_calls + 1L
  if (validator_calls %% 200L == 0L) {
    cat("Objective evaluations:", validator_calls, "at", format(Sys.time()), "\n"); flush.console()
  }
  tryCatch({emhawkes::dtzexp(1, a = param["a"], ell = param["ell"]); TRUE},
    error = function(e) FALSE)
}
legacy_objective <- function(object, inter_arrival, type, param) {
  if (!valid(param)) return(-Inf)
  suppressMessages(stats::logLik(object, inter_arrival = inter_arrival,
    type = type, lambda_component0 = 0))
}
# Warm both public likelihood routes; no fitting or initialization is timed.
for (i in 1:3) {
  a <- stats::logLik(legacy, path$inter_arrival, path$type, lambda_component0 = 0)
  b <- stats::logLik(model, path$inter_arrival, path$type, lambda_component0 = 0)
  stopifnot(identical(as.numeric(a), as.numeric(b)))
}
saved_rng <- .Random.seed
rows <- list(); estimates <- list(); fits <- list(); comparisons <- list()
for (repeat_id in seq_len(repetitions)) {
  order <- if (repeat_id %% 2L) c("hfit", "ffit") else c("ffit", "hfit")
  paired <- list()
  for (interface in order) {
    cat("Starting repetition", repeat_id, interface, "at", format(Sys.time()), "\n"); flush.console()
    validator_calls <- 0L
    warnings <- character()
    gc()
    elapsed <- system.time(fit <- withCallingHandlers(
      suppressMessages(if (interface == "hfit") emhawkes::hfit(legacy,
        path$inter_arrival, path$type, lambda_component0 = 0,
        mylogLik = legacy_objective, optimizer = "maxLik", method = "BFGS",
        reduced = TRUE, constraints = constraints, final_hessian = "optimizer",
        control = control) else ffit(model, path$inter_arrival, path$type,
        lambda_component0 = 0, estimate_residual = TRUE, valid = valid,
        optimizer = "maxLik", method = "BFGS", backend = "bridge",
        reduced = TRUE, constraints = constraints, control = control)),
      warning = function(w) {warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")}
    ))[["elapsed"]]
    stopifnot(identical(saved_rng, .Random.seed))
    raw <- if (interface == "ffit") fit$fit else fit
    restored <- if (interface == "ffit") fitted_fspec(fit) else emhawkes::fitted_hspec(fit)
    initial <- if (interface == "ffit") finitial_state(fit) else emhawkes::fit_initial_state(fit)
    inferred <- do.call(stats::logLik, c(list(object = restored,
      inter_arrival = path$inter_arrival, type = path$type, infer = TRUE), initial))
    ll <- as.numeric(stats::logLik(fit))
    stopifnot(isTRUE(all.equal(ll, inferred$loglikelihood, tolerance = 1e-10)))
    H <- maxLik::hessian(raw)
    information_eigen <- eigen(-(H+t(H))/2, symmetric = TRUE, only.values = TRUE)$values
    diagnosis <- if (interface == "ffit") fdiagnostics(fit) else NULL
    covariance <- suppressWarnings(stats::vcov(fit))
    se <- suppressWarnings(sqrt(diag(covariance)))
    coefficient <- stats::coef(fit)
    row <- data.frame(repetition = repeat_id, interface = interface, elapsed_seconds = elapsed,
      loglik = ll, restored_loglik = inferred$loglikelihood,
      restoration_difference = inferred$loglikelihood - ll,
      optimizer_code = maxLik::returnCode(raw),
      optimizer_message = maxLik::returnMessage(raw), validator_calls = validator_calls,
      iterations = raw$iterations, max_abs_score = max(abs(maxLik::gradient(raw))),
      minimum_information_eigenvalue = min(information_eigen),
      nonfinite_standard_errors = sum(!is.finite(se)),
      wald = if (is.null(diagnosis)) "provider_covariance" else diagnosis$wald,
      reason = if (is.null(diagnosis$reason)) "" else diagnosis$reason)
    rows[[length(rows)+1L]] <- row
    estimates[[length(estimates)+1L]] <- data.frame(repetition = repeat_id,
      interface = interface, parameter = names(coefficient), estimate = unname(coefficient),
      standard_error = unname(se))
    paired[[interface]] <- fit
    fits[[paste(repeat_id,interface,sep="_")]] <- fit
    saveRDS(list(fit = fit, warnings = warnings, time = elapsed, initial = initial,
      inference = inferred, information_eigenvalues = information_eigen),
      file.path(out, paste0(interface,"-",repeat_id,".rds")))
    write.csv(do.call(rbind,rows), file.path(out,"fits.csv"), row.names = FALSE)
    write.csv(do.call(rbind,estimates), file.path(out,"coefficients.csv"), row.names = FALSE)
    cat("Finished:"); print(row); if (length(warnings)) print(unique(warnings)); flush.console()
  }
  x <- paired$hfit; y <- paired$ffit
  coefficient_difference <- max(abs(stats::coef(x)-stats::coef(y)))
  likelihood_difference <- as.numeric(stats::logLik(y))-as.numeric(stats::logLik(x))
  comparisons[[repeat_id]] <- data.frame(repetition=repeat_id,
    coefficients_exact=identical(stats::coef(x),stats::coef(y)),
    maximum_coefficient_difference=coefficient_difference, likelihood_difference=likelihood_difference,
    raw_hessian_exact=identical(maxLik::hessian(x),maxLik::hessian(y$fit)),
    initial_exact=identical(emhawkes::fit_initial_state(x),finitial_state(y)))
  write.csv(do.call(rbind,comparisons),file.path(out,"comparisons.csv"),row.names=FALSE)
  print(comparisons[[repeat_id]])
  stopifnot(identical(names(stats::coef(x)),names(stats::coef(y))),
    isTRUE(all.equal(stats::coef(x),stats::coef(y),tolerance=1e-10)),
    isTRUE(all.equal(as.numeric(stats::logLik(x)),as.numeric(stats::logLik(y)),tolerance=1e-10)),
    identical(emhawkes::fit_initial_state(x),finitial_state(y)),
    isTRUE(fdiagnostics(y)$optimizer_success))
}
all_rows <- do.call(rbind,rows)
write.csv(aggregate(elapsed_seconds~interface,all_rows,median),
  file.path(out,"median-times.csv"),row.names=FALSE)
saveRDS(fits,file.path(out,"all-fits.rds"))
capture.output(sessionInfo(),file=file.path(out,"session-info.txt"))
cat("PASS: original and migrated 10000-row fits, restoration and RNG checks.\n")