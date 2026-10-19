# Choose the entire event loop before any seed change or user RNG call.
frpp_simulation_reason <- function(object, extra) {
  api <- c("exp_kernel_step", "exp_kernel_wait")
  if (!all(api %in% getNamespaceExports("emhawkes"))) return("the event-wise emhawkes kernel API is unavailable")
  h <- object@kernel
  values <- list(h@mu, h@alpha, h@beta, h@eta)
  if (any(vapply(values, function(x) !is.null(x) &&
      (!is.numeric(x) || is.object(x)), logical(1)))) return("kernel parameters must be plain numeric matrices")
  if (!is.null(h@impact) || !is.null(h@rmark) || !is.null(h@dmark))
    return("event callbacks require the bridge")
  if (extra > 0L) return("additional simulation arguments require the bridge")
  NULL
}

frpp_simulation_components_valid <- function(component) {
  if (any(!is.finite(component))) return(FALSE)
  if (all(component >= 0)) return(TRUE)
  for (i in seq_len(nrow(component))) {
    weight <- component[i, ]
    if (all(weight >= 0)) next
    scale <- sum(abs(weight)); total <- sum(weight)
    if (!is.finite(scale) || !is.finite(total) ||
        total < -32 * .Machine$double.eps * max(1, scale)) return(FALSE)
  }
  TRUE
}

prepare_frpp_simulation <- function(object, size, initial, N0, Nc0) {
  if (!is.numeric(size) || length(size) != 1L || !is.finite(size) ||
      size < 1 || size != floor(size))
    stop("size must be a positive integer, including the initial row.", call. = FALSE)
  kernel <- getExportedValue("emhawkes", "as_exp_kernel")(object@kernel)
  d <- kernel$dimens
  counts <- function(x, name) {
    if (is.null(x)) return(rep(0, d))
    if (!is.numeric(x) || !(length(x) %in% c(1L, d)) || any(!is.finite(x)))
      stop(name, " must be a finite scalar or contain one value per dimension.", call. = FALSE)
    rep_len(as.numeric(x), d)
  }
  N0 <- counts(N0, "N0"); Nc0 <- counts(Nc0, "Nc0")
  selection <- NULL
  if (is.null(initial)) initial <- "auto"
  if (is.character(initial)) {
    if (length(initial) != 1L || is.na(initial) || !initial %in% c("auto", "mean", "zero"))
      stop("lambda_component0 must be 'auto', 'mean', 'zero', a finite scalar, or a kernel-shaped matrix.", call. = FALSE)
    if (initial == "mean")
      stop("Mean initialization is unavailable: no supported unconditional mean formula for this model. Supply lambda_component0 explicitly or use 'zero'.", call. = FALSE)
    selection <- paste0("Initial intensity components: ",
      if (initial == "auto") "auto -> zero (no supported unconditional mean formula for this model)" else "zero",
      ". Supplying lambda_component0 explicitly is recommended for this model.")
    initial <- 0
  }
  # Preserve legacy input shape and lazy argument evaluation inside the RNG scope.
  if (!is.numeric(initial) || !length(initial) || any(!is.finite(initial)))
    stop("lambda_component0 must contain finite numeric values.", call. = FALSE)
  count <- d * kernel$kernel_cols
  shape_ok <- is.null(dim(initial)) || identical(dim(initial), dim(kernel$beta))
  if (length(initial) != 1L && (length(initial) != count || !shape_ok)) {
    received <- if (is.null(dim(initial))) paste0("length ", length(initial)) else paste(dim(initial), collapse = " x ")
    warning("lambda_component0 does not match the kernel dimensions (", d, " x ", kernel$kernel_cols,
      "; received ", received, "); recycling or truncating it is deprecated. Supply a scalar or a matrix with the correct dimensions.", call. = FALSE)
    initial <- rep_len(as.numeric(initial), count)
  }
  initial <- matrix(initial, d, kernel$kernel_cols)
  mu <- object@kernel@mu
  # Fix the kernel/baseline once. Older event-wise APIs retain the checked path.
  # Invalid baselines keep their original error timing at the first wait call.
  has_plan <- "exp_kernel_plan" %in% getNamespaceExports("emhawkes")
  valid_mu <- is.numeric(mu) && !is.object(mu) && length(mu) %in% c(1L, d) &&
    all(is.finite(mu)) && all(mu >= 0)
  if (has_plan && valid_mu) {
    prepared <- getExportedValue("emhawkes", "exp_kernel_plan")(kernel, mu)
    step <- prepared$step; wait <- prepared$wait
  } else {
    step_api <- getExportedValue("emhawkes", "exp_kernel_step")
    wait_api <- getExportedValue("emhawkes", "exp_kernel_wait")
    step <- function(component, dt, type = NULL, mark = 1, integrate = TRUE)
      step_api(kernel, component, dt, type, mark, integrate)
    wait <- function(component, residual) wait_api(kernel, component, mu, residual)
  }
  component <- step(initial, 0, integrate = FALSE)$before
  if (!frpp_simulation_components_valid(component))
    stop("Initial kernel components must define nonnegative, decreasing excitation; signed components must satisfy the decay-rate cumulative-sum condition (see ?hsim).", call. = FALSE)
  # Reuse only the public residual-function binding contract, not its simulator.
  rng <- as_legacy_hspec(object)@rresidual
  list(kernel = kernel, mu = object@kernel@mu, component = component,
    N0 = N0, Nc0 = Nc0, size = size, selection = selection, rng = rng,
    draw_size = object@kernel@dimens,
    step = step, wait = wait)
}

simulate_frpp_state <- function(object, plan) {
  if (!is.null(plan$selection)) message(plan$selection)
  kernel <- plan$kernel; mu <- plan$mu
  d <- kernel$dimens; k <- kernel$kernel_cols; size <- plan$size
  component <- plan$component
  inter_arrival <- type <- numeric(size)
  mark <- c(0, rep(1, size - 1L))
  N <- Nc <- matrix(0, size, d)
  colnames(N) <- paste0("N", seq_len(d)); colnames(Nc) <- paste0("Nc", seq_len(d))
  N[1L, ] <- plan$N0; Nc[1L, ] <- plan$Nc0
  lambda <- rambda <- matrix(0, size, d)
  lambda_component <- rambda_component <- matrix(0, size, d * k)
  lambda[1L, ] <- rambda[1L, ] <- mu + rowSums(component)
  lambda_component[1L, ] <- rambda_component[1L, ] <- t(component)
  colnames(lambda) <- paste0("lambda", seq_len(d))
  colnames(rambda) <- paste0("rambda", seq_len(d))
  colnames(lambda_component) <- paste0("lambda", rep(seq_len(d), each = k), rep(seq_len(k), d))
  colnames(rambda_component) <- sub("^lambda", "rambda", colnames(lambda_component))
  generated <- 1L
  for (n in seq_len(size - 1L) + 1L) {
    # A bad post-event state must fail before consuming the next residual vector.
    if (!frpp_simulation_components_valid(component))
      stop("Simulation requires finite, nonnegative baselines and nonnegative, decreasing excitation; signed components must satisfy the decay-rate cumulative-sum condition (see ?hsim).", call. = FALSE)
    residual <- plan$rng(plan$draw_size)
    candidates <- plan$wait(component, residual)
    dt <- min(candidates, na.rm = TRUE)
    if (is.infinite(dt)) {
      warning(sprintf("Simulation stopped early: no further event occurs on this path. Returned %d of %.0f requested rows (including the initial row).", generated, size), call. = FALSE)
      break
    }
    event <- which(candidates == dt)[1L]
    inter_arrival[n] <- dt; type[n] <- event; generated <- n
    N[n, ] <- N[n - 1L, ]; N[n, event] <- N[n - 1L, event] + 1
    Nc[n, ] <- Nc[n - 1L, ]; Nc[n, event] <- Nc[n - 1L, event] + mark[n]
    advanced <- plan$step(component, dt, event, mark[n], integrate = FALSE)
    lambda_component[n, ] <- t(advanced$before)
    lambda[n, ] <- mu + rowSums(advanced$before)
    rambda_component[n, ] <- t(advanced$after)
    rambda[n, ] <- mu + rowSums(advanced$after)
    # Match history storage: no dimnames travel into the next event state.
    component <- matrix(advanced$after, nrow = d)
  }
  rows <- seq_len(generated)
  structure(list(model = object, inter_arrival = inter_arrival[rows],
    arrival = cumsum(inter_arrival[rows]), type = type[rows], mark = mark[rows],
    N = N[rows, , drop = FALSE], Nc = Nc[rows, , drop = FALSE],
    lambda = lambda[rows, , drop = FALSE], lambda_component = lambda_component[rows, , drop = FALSE],
    rambda = rambda[rows, , drop = FALSE], rambda_component = rambda_component[rows, , drop = FALSE]),
    class = "frpp_path")
}
