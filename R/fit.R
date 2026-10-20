# Freeze distribution parameters without changing user functions or their
# original environments. Empty defaults keep them out of hfit's parameter map.
fixed_residual_spec <- function(model) {
  object <- as_legacy_hspec(model)
  density <- object@dresidual; cdf <- object@presidual
  quantile <- object@qresidual; rng <- object@rresidual
  if (!is.null(density)) object@dresidual <- function(x, param = numeric(), log = FALSE)
    density(x, log = log)
  if (!is.null(cdf)) object@presidual <- function(q, param = numeric(), lower.tail = TRUE, log.p = FALSE)
    cdf(q, lower.tail = lower.tail, log.p = log.p)
  if (!is.null(quantile)) object@qresidual <- function(p, param = numeric(), lower.tail = TRUE, log.p = FALSE)
    quantile(p, lower.tail = lower.tail, log.p = log.p)
  if (!is.null(rng)) object@rresidual <- function(size, param = numeric()) rng(size)
  object
}

# Explicit numeric registration. The bridge suppresses the sharing warning
# already emitted by hfit; independent fitting emits it during registration.
frpp_parameter_plan <- function(object, reduced, estimate_residual, warn_shared = FALSE) {
  if (!"parameter_plan" %in% getNamespaceExports("emhawkes")) return(NULL)
  h <- object@kernel; law <- object@residual
  values <- list(mu=h@mu, alpha=h@alpha, beta=h@beta, eta=h@eta)
  if (law$binding != "raw" || is.null(law$functions$density) ||
      !is.null(h@impact) || !is.null(h@rmark) || !is.null(h@dmark) ||
      any(vapply(values,function(x) !is.null(x) &&
        (!is.numeric(x) || is.complex(x) || is.object(x)),logical(1)))) return(NULL)
  values$residual <- if (estimate_residual) law$param else numeric()
  # hfit already emitted the legacy sharing warning before optimization.
  withCallingHandlers(getExportedValue("emhawkes","parameter_plan")(
    values, named="residual", reduced=reduced), warning=function(w) {
      if (!warn_shared && startsWith(conditionMessage(w),
          "Shared parameters have different starting values; using the first value for: "))
        invokeRestart("muffleWarning")
    })
}

restore_numeric_frpp <- function(object, plan, coefficients, estimate_residual) {
  values <- plan$restore(coefficients)
  kernel <- methods::new("hspec",mu=values$mu,alpha=values$alpha,
    beta=values$beta,eta=values$eta,type_col_map=object@kernel@type_col_map)
  law <- object@residual
  if (estimate_residual) {
    shared <- intersect(names(law$param),names(coefficients))
    law$param[shared] <- values$residual[shared]
  }
  fspec(kernel,law)
}

ffit <- function(object, inter_arrival, type = NULL, mark = NULL,
                 N = NULL, Nc = NULL, N0 = NULL, Nc0 = NULL,
                 lambda_component0 = "auto", estimate_residual = FALSE,
                 optimizer = c("maxLik", "stats"), method = "BFGS",
                 backend = c("bridge", "state"), reduced = TRUE,
                 fixed = NULL, constraints = NULL, constraint = NULL,
                 valid = NULL, ..., fit_backend = c("bridge", "independent", "auto")) {
  if (!is_fspec(object)) stop("object must be an fspec.", call. = FALSE)
  object <- as_fspec(object)
  require_law_functions(object, c("density", if (object@kernel@dimens > 1) "cdf"))
  if (!is.logical(estimate_residual) || length(estimate_residual) != 1L || is.na(estimate_residual))
    stop("estimate_residual must be TRUE or FALSE.", call. = FALSE)
  if (!is.null(valid) && !is.function(valid)) stop("valid must be a function or NULL.", call. = FALSE)
  optimizer <- match.arg(optimizer); backend <- match.arg(backend)
  fit_backend <- match.arg(fit_backend)
  if (!is.null(constraint)) {
    if (!is.null(constraints)) stop("Supply only one of constraint or constraints.", call. = FALSE)
    constraints <- constraint
  }
  dots <- list(...)
  if (length(dots) && (is.null(names(dots)) || anyNA(names(dots)) ||
      any(!nzchar(names(dots))) || anyDuplicated(names(dots))))
    stop("Additional optimizer arguments must have unique, nonempty names.", call. = FALSE)
  reserved <- c("mylogLik", "gradient", "jacobian", "final_hessian", "lambda0")
  if (any(names(dots) %in% reserved)) stop("ffit controls its FRPP objective and curvature; unsupported argument(s): ",
    paste(intersect(names(dots), reserved), collapse = ", "), ".", call. = FALSE)
  if (is.null(lambda_component0) || identical(lambda_component0, "auto")) {
    lambda_component0 <- 0
    message("Initial intensity components: auto -> zero for FRPP. Supplying lambda_component0 explicitly is recommended.")
  }
  if (is.character(lambda_component0))
    stop("FRPP fitting requires explicit initial components or 'auto'; a Hawkes mean is not used.", call. = FALSE)
  fitting_reason <- if (fit_backend == "bridge") NULL else frpp_fitting_reason(object)
  if (fit_backend == "independent" && !is.null(fitting_reason))
    stop("Independent FRPP fitting is unavailable: ",fitting_reason,". Use fit_backend = 'bridge'.",call.=FALSE)
  independent <- fit_backend != "bridge" && is.null(fitting_reason)
  if (independent) {
    execution <- run_frpp_independent(object,inter_arrival,type,mark,N,Nc,N0,Nc0,
      lambda_component0,estimate_residual,reduced,optimizer,method,constraints,fixed,
      valid,backend,dots)
    fit <- execution$fit
  } else {
    legacy <- if (estimate_residual) as_legacy_hspec(object) else fixed_residual_spec(object)
    # A user validator receives fitted names, plus fixed raw residual parameters.
    # FALSE rejects a numerical proposal; errors are not swallowed.
    residual_parameters <- if (object@residual$binding == "raw") object@residual$param else numeric()
    objective <- function(object, inter_arrival, type, mark, N, Nc, N0, Nc0,
                          lambda_component0, param) {
      if (!is.null(valid)) {
        full <- c(param, residual_parameters[setdiff(names(residual_parameters), names(param))])
        accepted <- valid(full)
        if (!is.logical(accepted) || length(accepted) != 1L || is.na(accepted))
          stop("valid must return one nonmissing logical value.", call. = FALSE)
        if (!accepted) return(-Inf)
      }
      stats::logLik(as_fspec(object), inter_arrival, type, mark,
        N = N, Nc = Nc, N0 = N0, Nc0 = Nc0, lambda_component0 = lambda_component0,
        backend = backend)
    }
    # Only public emhawkes APIs are reused: its sharing map, optimizers and
    # reconstruction. The objective is explicitly FRPP and never Hawkes curvature.
    args <- c(list(object = legacy, inter_arrival = inter_arrival, type = type, mark = mark,
      N = N, Nc = Nc, N0 = N0, Nc0 = Nc0, lambda_component0 = lambda_component0,
      mylogLik = objective, reduced = reduced, optimizer = optimizer, method = method,
      constraints = constraints, fixed = fixed, final_hessian = "optimizer"), dots)
    fit <- do.call(emhawkes::hfit, args)
  }
  if (is.null(fit)) stop("The optimizer did not return a fit.", call. = FALSE)
  estimate <- stats::coef(fit)
  registration <- if (independent) execution$registration else
    frpp_parameter_plan(object, reduced, estimate_residual)
  if (!is.null(registration)) {
    restored <- restore_numeric_frpp(object, registration, estimate, estimate_residual)
  } else {
    restored <- as_fspec(emhawkes::fitted_hspec(fit))
    if (object@residual$binding == "raw") {
      law <- object@residual
      if (estimate_residual) {
        shared <- intersect(names(law$param), names(estimate))
        law$param[shared] <- estimate[shared]
      }
      restored <- fspec(restored@kernel, law)
    } else if (!estimate_residual) restored <- fspec(restored@kernel, object@residual)
  }
  active <- if (optimizer == "stats") fit$active else maxLik::activePar(fit)
  code <- if (optimizer == "stats") fit$code else maxLik::returnCode(fit)
  message <- if (optimizer == "stats") fit$message else maxLik::returnMessage(fit)
  success <- if (optimizer == "stats") code == 0L else if (tolower(method) %in%
    c("nr", "newton-raphson", "bfgsr", "bfgs-r", "bhhh")) code %in% c(1L, 2L, 8L) else
    if (tolower(method) == "sann") NA else code == 0L
  covariance <- tryCatch(suppressWarnings(stats::vcov(fit)), error = function(e) NULL)
  reason <- NULL
  free <- which(active)
  if (!is.null(constraints$eqA) || !is.null(constraints$eqB)) reason <- "equality_constraints"
  if (is.null(reason) && !is.null(constraints$ineqA) && length(free)) {
    A <- constraints$ineqA; B <- constraints$ineqB
    distance <- drop(A %*% estimate + B)
    scale <- drop(abs(A) %*% pmax(abs(estimate), 1)) + abs(B)
    movable <- rowSums(abs(A[, active, drop = FALSE])) > 0
    if (any(distance[movable] <= sqrt(.Machine$double.eps) * pmax(scale[movable], 1)))
      reason <- "constraint_boundary"
  }
  if (is.null(reason) && optimizer == "stats" && !is.null(fit$boundary)) reason <- "constraint_boundary"
  if (is.null(reason) && !isTRUE(success)) reason <- "optimizer_not_converged"
  if (is.null(reason)) {
    if (!is.matrix(covariance) || !identical(dim(covariance), c(length(estimate), length(estimate))) ||
        any(!is.finite(covariance[active, active, drop = FALSE]))) reason <- "covariance_unavailable"
    else if (length(free)) {
      block <- covariance[active, active, drop = FALSE]
      if (!isTRUE(all.equal(block, t(block), tolerance = 1e-10)) ||
          inherits(tryCatch(chol(block), error = identity), "error")) reason <- "nonpositive_curvature"
    }
  }
  if (!is.null(reason)) {
    covariance <- matrix(0, length(estimate), length(estimate), dimnames = list(names(estimate), names(estimate)))
    covariance[active, active] <- NA_real_
    warning("FRPP Wald standard errors are unavailable (", reason,
      "). Inspect fdiagnostics(); consider profile likelihood or bootstrap.", call. = FALSE)
  }
  diagnostics <- list(optimizer = optimizer, method = method, backend = backend,
    optimizer_code = code, optimizer_message = message, optimizer_success = success,
    estimate_residual = estimate_residual, wald = if (is.null(reason)) "available_unvalidated" else "unavailable",
    reason = reason, active = active,
    note = "Optimizer curvature is not a validation of FRPP Wald coverage; no Hawkes Hessian stabilization is applied.")
  structure(list(fit = fit, model = restored, initial = if (independent) execution$initial else emhawkes::fit_initial_state(fit),
    coefficients = estimate, covariance = covariance, maximum = as.numeric(stats::logLik(fit)),
    nobs = length(inter_arrival) - 1L, diagnostics = diagnostics, call = match.call(),
    fitting = list(requested=fit_backend,used=if(independent)"independent" else "bridge",
      fallback_reason=fitting_reason)), class = "ffit")
}

fitted_fspec <- function(object) {
  if (!inherits(object, "ffit")) stop("object must be an ffit result.", call. = FALSE)
  as_fspec(object$model)
}
fdiagnostics <- function(object) {
  if (!inherits(object, "ffit")) stop("object must be an ffit result.", call. = FALSE)
  object$diagnostics
}
finitial_state <- function(object) {
  if (!inherits(object, "ffit")) stop("object must be an ffit result.", call. = FALSE)
  object$initial
}
coef.ffit <- function(object, ...) object$coefficients
vcov.ffit <- function(object, ...) object$covariance
logLik.ffit <- function(object, ...) structure(object$maximum, class = "logLik",
  df = sum(object$diagnostics$active), nobs = object$nobs)
nobs.ffit <- function(object, ...) object$nobs
fitted.ffit <- function(object, ...) fitted_fspec(object)
summary.ffit <- function(object, ...) {
  se <- sqrt(diag(stats::vcov(object)))
  z <- rep(NA_real_, length(se)); names(z) <- names(se)
  usable <- object$diagnostics$active & is.finite(se) & se > 0
  z[usable] <- stats::coef(object)[usable] / se[usable]
  structure(list(coefficients = cbind(Estimate = stats::coef(object), `Std. Error` = se,
    `z value` = z, `Pr(>|z|)` = 2 * stats::pnorm(-abs(z))), maximum = object$maximum,
    diagnostics = object$diagnostics), class = "summary.ffit")
}
print.summary.ffit <- function(x, ...) {
  d <- x$diagnostics
  cat("FRPP maximum likelihood fit (", d$optimizer, ", ", d$method, ", ", d$backend, ")\n", sep = "")
  cat("Log likelihood: ", format(x$maximum), "\nConvergence code: ", d$optimizer_code,
    " (", d$optimizer_message, ")\n", sep = "")
  stats::printCoefmat(x$coefficients, ...)
  if (any(!d$active)) cat("Fixed parameters: ", paste(rownames(x$coefficients)[!d$active], collapse = ", "), "\n", sep = "")
  cat("Wald covariance: ", d$wald, if (!is.null(d$reason)) paste0(" (", d$reason, ")"), "\n", sep = "")
  invisible(x)
}
print.ffit <- function(x, ...) { print(summary(x), ...); invisible(x) }
