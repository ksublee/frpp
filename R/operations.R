fsim <- function(object, size = 100, lambda_component0 = "auto", N0 = NULL,
                  Nc0 = NULL, seed = NULL, ...) {
  if (!is_fspec(object))
    stop("object must be an fspec or frpp_spec.", call. = FALSE)
  object <- as_fspec(object)
  require_law_functions(object, "rng")
  # Validate support before setting a seed or entering the event loop.
  beta <- object@kernel@beta
  if (is.matrix(beta) && any(vapply(seq_len(nrow(beta)),
    function(i) any(beta[i, ] != beta[i, 1]), logical(1))))
    stop("FRPP simulation requires equal beta values within each row.", call. = FALSE)
  legacy <- as_legacy_hspec(object)
  if (!is.null(seed)) {
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
    old_kind <- RNGkind()
    on.exit({
      do.call(RNGkind, as.list(old_kind))
      if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv) else
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
          rm(list = ".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
    set.seed(seed)
  }
  path <- emhawkes::hsim(legacy, size = size,
    lambda_component0 = lambda_component0, N0 = N0, Nc0 = Nc0, ...)
  values <- unclass(path)
  values$hspec <- NULL
  structure(c(list(model = object), values), class = "frpp_path")
}

simulate.frpp_spec <- function(object, nsim = 1, seed = NULL, size = 100,
                               lambda_component0 = "auto", N0 = NULL,
                               Nc0 = NULL, ...) {
  if (!is.numeric(nsim) || length(nsim) != 1L || is.na(nsim) || nsim != 1)
    stop("This prototype returns one path; use nsim = 1 and specify size.", call. = FALSE)
  fsim(object, size = size, lambda_component0 = lambda_component0,
       N0 = N0, Nc0 = Nc0, seed = seed, ...)
}

print.frpp_path <- function(x, ...) {
  cat("FRPP path: ", length(x$inter_arrival) - 1L, " events, ",
      fkernel(x$model)@dimens, " event types.\n", sep = "")
  invisible(x)
}

logLik.frpp_spec <- function(object, inter_arrival, type = NULL, mark = NULL,
                            N = NULL, Nc = NULL, N0 = NULL, Nc0 = NULL,
                            lambda_component0 = "auto", infer = FALSE,
                            backend = c("bridge", "state"), ...) {
  object <- as_fspec(object)
  backend <- match.arg(backend)
  if (!is.logical(infer) || length(infer) != 1L || is.na(infer))
    stop("infer must be TRUE or FALSE.", call. = FALSE)
  require_law_functions(object, c("density",
    if (object@kernel@dimens > 1 || infer) "cdf", if (infer) "quantile"))
  if (backend == "state") return(state_loglik(object, inter_arrival, type,
    mark, N, Nc, N0, Nc0, lambda_component0, infer, ...))
  emhawkes::logLik(as_legacy_hspec(object), inter_arrival = inter_arrival,
    type = type, mark = mark, N = N, Nc = Nc, N0 = N0, Nc0 = Nc0,
    lambda_component0 = lambda_component0, infer = infer, ...)
}

simulate.fspec <- simulate.frpp_spec
logLik.fspec <- logLik.frpp_spec

# The S3 registrations serve stats:: calls and saved S3 objects. The S4
# methods serve the formal generics, including emhawkes::logLik().
methods::setMethod("simulate", "fspec",
  function(object, nsim = 1, seed = NULL, ...)
    simulate.fspec(object, nsim = nsim, seed = seed, ...))
methods::setMethod("logLik", "fspec",
  function(object, ...) logLik.fspec(object, ...))
