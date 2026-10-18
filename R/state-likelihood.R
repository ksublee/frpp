# Numeric FRPP likelihood using the public, baseline-independent kernel API.
# Keep the bridge for callback models and installations without that API.
state_loglik <- function(object, inter_arrival, type, mark, N, Nc, N0, Nc0,
                         initial, infer, ...) {
  kernel <- object@kernel
  values <- list(kernel@mu, kernel@alpha, kernel@beta, kernel@eta)
  if (any(vapply(values, function(x) !is.null(x) &&
    (!is.numeric(x) || is.object(x)), logical(1))) ||
      !is.null(kernel@impact) || !is.null(kernel@dmark))
    stop("The state backend currently requires numeric kernels without impact or mark-density callbacks.", call. = FALSE)
  if (any(vapply(seq_len(nrow(kernel@beta)),
    function(i) any(kernel@beta[i, ] != kernel@beta[i, 1]), logical(1))))
    stop("The state prototype requires equal beta values within each row.", call. = FALSE)
  if (length(list(...)))
    stop("The state prototype does not accept additional arguments.", call. = FALSE)
  if (is.null(initial) || identical(initial, "auto")) {
    initial <- 0
    message("Initial intensity components: auto -> zero for FRPP. ",
            "Supplying lambda_component0 explicitly is recommended.")
  }
  if (identical(initial, "mean"))
    stop("A Hawkes unconditional mean is not established for this FRPP law; supply an explicit state.", call. = FALSE)
  if (identical(initial, "zero")) initial <- 0
  if (!all(c("as_exp_kernel", "exp_kernel_replay") %in%
           getNamespaceExports("emhawkes")))
    stop("The state backend requires the emhawkes development kernel API; use backend = 'bridge' with earlier versions.", call. = FALSE)
  replay <- getExportedValue("emhawkes", "exp_kernel_replay")
  path <- replay(kernel, inter_arrival, type, mark,
    lambda_component0 = initial, history = infer)
  count <- nrow(path$kernel_before)
  dims <- as.integer(kernel@dimens)
  events <- if (is.null(type)) c(0, rep(1, count - 1L)) else type
  marks <- if (is.null(mark)) c(0, rep(1, count - 1L)) else mark
  validate_state_counts(events, marks, dims, N, Nc, N0, Nc0)
  # Delegate legacy function signatures and log/tail compatibility to the public
  # adapter; this does not execute callbacks or alter their default bindings.
  bound <- as_legacy_hspec(object)
  law <- emhawkes::set_residual(kernel, param = numeric(),
    dresidual = bound@dresidual, presidual = bound@presidual,
    qresidual = bound@qresidual)
  phi <- matrix(0, count, dims,
    dimnames = list(NULL, paste0("integrated_driver", seq_len(dims))))
  result <- 0
  mu <- as.vector(kernel@mu)
  lambda <- sweep(path$kernel_before, 2L, mu, "+")
  for (n in seq_len(count - 1L) + 1L) {
    phi[n, ] <- path$kernel_integral[n, ] + mu * inter_arrival[n]
    if (any(!is.finite(phi[n, ])) || any(phi[n, ] < 0))
      stop("The state prototype requires finite nonnegative integrated drivers.", call. = FALSE)
    contribution <- numeric(dims)
    for (k in seq_len(dims)) contribution[k] <- if (k == events[n]) {
      law@dresidual(phi[n, k], log = TRUE)
    } else {
      law@presidual(phi[n, k], lower.tail = FALSE, log.p = TRUE)
    }
    result <- result + log(unname(lambda[n, events[n]])) + sum(contribution)
  }
  if (!infer) return(c(loglikelihood = result))
  colnames(lambda) <- paste0("lambda", seq_len(dims))
  after <- matrix(0, count, dims)
  cols <- ncol(kernel@beta)
  for (k in seq_len(dims)) after[, k] <- mu[k] + rowSums(
    path$component_after[, (k - 1L) * cols + seq_len(cols), drop = FALSE])
  colnames(after) <- paste0("rambda", seq_len(dims))
  colnames(path$component_before) <- paste0("lambda", rep(seq_len(dims), each = cols), seq_len(cols))
  colnames(path$component_after) <- paste0("rambda", rep(seq_len(dims), each = cols), seq_len(cols))
  answer <- list(loglikelihood = result, lambda = lambda,
    lambda_component = path$component_before, rambda = after,
    rambda_component = path$component_after)
  for (k in seq_len(dims)) {
    indices <- which(events == k & seq_along(events) > 1L)
    residual <- exponential <- numeric(max(0L, length(indices) - 1L))
    for (j in seq_along(residual)) {
      interval <- seq.int(indices[j] + 1L, indices[j + 1L])
      log_survival <- sum(law@presidual(phi[interval, k],
        lower.tail = FALSE, log.p = TRUE))
      residual[j] <- law@qresidual(log_survival, lower.tail = FALSE, log.p = TRUE)
      exponential[j] <- -log_survival
    }
    answer[[paste0("res_process", k)]] <- residual
    answer[[paste0("eres_process", k)]] <- exponential
  }
  answer$integrated_driver <- phi
  answer
}

# Count paths do not enter numeric kernel recursion, but remain part of the
# model's observation contract. Match the public bridge's validation here.
validate_state_counts <- function(type, mark, dims, N, Nc, N0, Nc0) {
  count <- length(type)
  for (name in c("N", "Nc")) {
    value <- if (name == "N") N else Nc
    initial <- if (name == "N") N0 else Nc0
    if (!is.null(value)) {
      if (!is.matrix(value) || !is.numeric(value) ||
          !identical(dim(value), c(count, dims)) || any(!is.finite(value)))
        stop(name, " must be a finite matrix with one row per observation and one column per dimension.", call. = FALSE)
      if (is.null(initial)) initial <- value[1L, ]
    }
    if (is.null(initial)) initial <- rep(0, dims)
    if (!is.numeric(initial) || !length(initial) || any(!is.finite(initial)) ||
        !(length(initial) %in% c(1L, dims)))
      stop(name, "0 must be a finite scalar or one value per event type.", call. = FALSE)
    if (!is.null(value)) {
      expected <- matrix(0, count, dims)
      expected[1L, ] <- initial
      rows <- seq_len(count - 1L) + 1L
      expected[cbind(rows, type[rows])] <- if (name == "N") 1 else mark[rows]
      for (k in seq_len(dims)) expected[, k] <- cumsum(expected[, k])
      if (!isTRUE(all.equal(as.numeric(value), as.numeric(expected),
                           tolerance = sqrt(.Machine$double.eps))))
        stop(name, " is inconsistent with type, mark, or its initial values.", call. = FALSE)
    }
  }
  invisible(NULL)
}
