# Observation and initialization arithmetic follows emhawkes 18986df.
# FRPP owns these fitting inputs; no internal emhawkes API is used.

frpp_type_to_N <- function(type, dimens, N0 = NULL){

  N <- matrix(rep(0, dimens*length(type)), ncol=dimens)
  colnames(N) <- paste0("N", 1:dimens)

  if (!is.null(N0)) N[1,] <- N0

  events <- seq_len(max(0L, length(type) - 1L)) + 1L
  N[cbind(events, type[events])] <- 1
  for (j in seq_len(dimens)) N[, j] <- cumsum(N[, j])

  N

}

frpp_type_mark_to_Nc <- function(type, mark, dimens, Nc0 = NULL){

  Nc <-  matrix(rep(0, dimens*length(type)), ncol=dimens)
  colnames(Nc) <- paste0("Nc", 1:dimens)

  if (!is.null(Nc0)) Nc[1,] <- Nc0

  events <- seq_len(max(0L, length(type) - 1L)) + 1L
  Nc[cbind(events, type[events])] <- mark[events]
  for (j in seq_len(dimens)) Nc[, j] <- cumsum(Nc[, j])

  Nc

}

frpp_count_initial <- function(value, dimens, name) {
  if (is.null(value)) return(rep(0, dimens))
  if (!is.numeric(value) || !(length(value) %in% c(1L, dimens)) || any(!is.finite(value))) {
    stop(name, " must be a finite scalar or contain one value per dimension.")
  }
  rep_len(as.numeric(value), dimens)
}

frpp_fit_observation <- function(dims, inter_arrival, type = NULL, mark = NULL,
                                  N = NULL, Nc = NULL, N0 = NULL, Nc0 = NULL,
                                  needs = c("N", "Nc")) {
  size <- length(inter_arrival)
  dimens <- dims
  if (!is.numeric(inter_arrival) || !is.null(dim(inter_arrival)) || !size ||
      any(!is.finite(inter_arrival)) || any(inter_arrival < 0) || inter_arrival[1] != 0) {
    stop("inter_arrival must be a finite, nonnegative vector starting with an initial zero.")
  }
  if (is.null(type)) {
    if (dimens != 1L) stop("type is required for a multivariate model; supply one event type per row of inter_arrival, including the initial row.", call. = FALSE)
    type <- c(0, rep(1, size - 1L))
  }
  if (is.null(mark)) mark <- c(0, rep(1, size - 1L))
  if (!is.numeric(type) || !is.null(dim(type)) || length(type) != size || any(!is.finite(type)) ||
      any(type[-1] < 1 | type[-1] > dimens | type[-1] != floor(type[-1]))) {
    stop("type must match inter_arrival (", size, " rows) and contain integer event types from 1 to ",
         dimens, " after the initial row.", call. = FALSE)
  }
  if (!is.numeric(mark) || !is.null(dim(mark)) || length(mark) != size || any(!is.finite(mark))) {
    stop("mark must be finite and match inter_arrival.")
  }
  normalize_path <- function(path, initial, name) {
    if (!is.null(path)) {
      if (!is.matrix(path) || !is.numeric(path) ||
          !identical(dim(path), c(size, as.integer(dimens))) || any(!is.finite(path))) {
        stop(name, " must be a finite matrix with one row per observation and one column per dimension.")
      }
      if (is.null(initial)) initial <- path[1, ]
    }
    initial <- frpp_count_initial(initial, dimens, paste0(name, "0"))
    if (!is.null(path) || name %in% needs) {
      expected <- if (name == "N") frpp_type_to_N(type, dimens, initial) else
        frpp_type_mark_to_Nc(type, mark, dimens, initial)
      if (is.null(path)) path <- expected else {
        if (!isTRUE(all.equal(as.numeric(path), as.numeric(expected), tolerance = sqrt(.Machine$double.eps)))) {
          stop(name, " is inconsistent with type, mark, or its initial values.")
        }
        if (is.null(colnames(path))) colnames(path) <- colnames(expected)
      }
    }
    list(path = path, initial = initial)
  }
  counts <- normalize_path(N, N0, "N")
  marked_counts <- normalize_path(Nc, Nc0, "Nc")
  list(inter_arrival = inter_arrival, type = type, mark = mark,
       N = counts$path, Nc = marked_counts$path, N0 = counts$initial, Nc0 = marked_counts$initial)
}

frpp_fit_initial <- function(value, dimens, kernel_cols, legacy = FALSE) {
  if (is.null(value)) return(NULL)
  size <- dimens * kernel_cols
  if (!is.numeric(value) || !length(value) || any(!is.finite(value))) {
    stop("lambda_component0 must contain finite numeric values.")
  }
  if (length(value) == 1L) return(matrix(value, dimens, kernel_cols))
  shape_ok <- is.null(dim(value)) || identical(dim(value), c(as.integer(dimens), as.integer(kernel_cols)))
  if (length(value) != size || !shape_ok) {
    expected <- paste0(dimens, " x ", kernel_cols)
    received <- if (is.null(dim(value))) paste0("length ", length(value)) else paste(dim(value), collapse = " x ")
    if (!legacy) stop("lambda_component0 must be a scalar or match the kernel matrix dimensions (",
                      expected, "; received ", received, ").", call. = FALSE)
    warning("lambda_component0 does not match the kernel dimensions (", expected, "; received ",
            received, "); recycling or truncating it is deprecated. Supply a scalar or a matrix with the correct dimensions.", call. = FALSE)
    value <- rep_len(as.numeric(value), size)
  }
  matrix(value, dimens, kernel_cols)
}

# Independent fitting is deliberately narrower than the compatibility bridge.
frpp_fitting_reason <- function(object) {
  if (!"parameter_plan" %in% getNamespaceExports("emhawkes"))
    return("the public emhawkes parameter_plan API is unavailable")
  h <- object@kernel; law <- object@residual
  values <- list(h@mu,h@alpha,h@beta,h@eta)
  if (law$binding != "raw") return("the residual law uses legacy bound functions")
  if (is.null(law$functions$density)) return("the residual density is missing")
  if (!is.null(h@impact) || !is.null(h@rmark) || !is.null(h@dmark))
    return("event callbacks require the bridge")
  if (any(vapply(values,function(x) !is.null(x) &&
    (!is.numeric(x) || is.complex(x) || is.object(x)),logical(1))))
    return("kernel parameters must be plain numeric values")
  NULL
}

run_frpp_independent <- function(object, inter_arrival, type, mark, N, Nc, N0, Nc0,
    initial, estimate_residual, reduced, optimizer, method, constraints, fixed,
    valid, backend, additional) {
  if (optimizer == "maxLik" && !requireNamespace("maxLik",quietly=TRUE))
    stop("optimizer = 'maxLik' requires the optional maxLik package. Install it or use optimizer = 'stats'.",call.=FALSE)
  registration <- frpp_parameter_plan(object,reduced,estimate_residual,warn_shared=TRUE)
  observation <- frpp_fit_observation(object@kernel@dimens,inter_arrival,type,mark,N,Nc,N0,Nc0)
  initial <- frpp_fit_initial(initial,object@kernel@dimens,ncol(object@kernel@beta))
  residual_parameters <- object@residual$param
  grad <- additional$grad; additional$grad <- NULL
  hess <- additional$hess; additional$hess <- NULL
  verbose <- additional$verbose; additional$verbose <- NULL
  objective <- function(param) {
    values <- registration$restore(param)
    if (any(!is.finite(values$beta)) || any(values$beta <= 0)) return(-Inf)
    for (value in values[c("mu","alpha","eta")])
      if (any(!is.finite(value))) return(-Inf)
    if (any(values$mu < 0)) return(-Inf)
    if (!is.null(valid)) {
      full <- c(param,residual_parameters[setdiff(names(residual_parameters),names(param))])
      accepted <- valid(full)
      if (!is.logical(accepted) || length(accepted)!=1L || is.na(accepted))
        stop("valid must return one nonmissing logical value.",call.=FALSE)
      if (!accepted) return(-Inf)
    }
    kernel <- methods::new("hspec",mu=values$mu,alpha=values$alpha,
      beta=values$beta,eta=values$eta,type_col_map=object@kernel@type_col_map)
    functions <- object@residual$functions
    names(functions) <- unname(residual_slots)
    # The public residual adapter preserves legacy log/tail/default conventions.
    bound <- do.call(emhawkes::set_residual,c(list(object=kernel,
      param=if(estimate_residual)values$residual else residual_parameters),functions))
    value <- do.call(stats::logLik,c(list(object=as_fspec(bound)),observation,
      list(lambda_component0=initial,backend=backend)))
    if (isTRUE(verbose)) {
      cat("Parameters : ",param,"\n");cat("Log likelihood : ",value,"\n")
    }
    value
  }
  # Match the bridge's existing matrix-sharing notice without exposing its map.
  shared <- vapply(list(mu=object@kernel@mu,alpha=object@kernel@alpha,
    beta=object@kernel@beta,eta=object@kernel@eta),function(x)
      !is.null(x) && isTRUE(reduced) && anyDuplicated(as.vector(x))>0L,logical(1))
  if (any(shared)) message("ffit: equal starting values within each matrix share a parameter ",
    "(reduced = TRUE; ",paste(names(shared)[shared],collapse=", "),").\n",
    "Use reduced = FALSE to estimate matrix entries separately.")
  if (optimizer == "maxLik") {
    fit <- do.call(maxLik::maxLik,c(list(logLik=objective,start=registration$start,
      grad=grad,hess=hess,method=method,constraints=constraints,fixed=fixed),additional))
  } else {
    fit <- run_frpp_stats_fit(objective,registration$start,grad,hess,constraints,
      method,c(list(fixed=fixed),additional),stable=FALSE)
    if (fit$code != 0L) warning("The stats optimizer did not report normal convergence (code ",
      fit$code,"). See fdiagnostics().",call.=FALSE)
  }
  list(fit=fit,registration=registration,
    initial=list(lambda_component0=initial,N0=observation$N0,Nc0=observation$Nc0))
}
