residual_slots <- c(density = "dresidual", cdf = "presidual",
                    quantile = "qresidual", rng = "rresidual")

# Creating a law never evaluates its callbacks or consumes random numbers.
residual_law <- function(param = numeric(), density = NULL, cdf = NULL,
                         quantile = NULL, rng = NULL, name = "custom") {
  if (!is.numeric(param) || is.complex(param) || !is.null(dim(param)) ||
      any(!is.finite(param)) || (length(param) &&
      (is.null(names(param)) || anyNA(names(param)) || any(names(param) == "") ||
       anyDuplicated(names(param)))))
    stop("param must be a finite numeric vector with unique, nonempty names.", call. = FALSE)
  callbacks <- list(density = density, cdf = cdf, quantile = quantile, rng = rng)
  if (!all(vapply(callbacks, function(x) is.null(x) || is.function(x), logical(1))))
    stop("Residual functions must be functions or NULL.", call. = FALSE)
  if (is.null(density) && is.null(rng))
    stop("Supply a density or a random generator for the residual law.", call. = FALSE)
  available <- unique(unlist(lapply(callbacks, function(f) {
    if (is.null(f)) character() else names(formals(f))[-1L]
  })))
  available <- setdiff(available, c("...", "log", "log.p", "lower.tail"))
  unused <- setdiff(names(param), available)
  if (length(unused)) stop("Residual parameter(s) do not match function arguments: ",
    paste(unused, collapse = ", "), ".", call. = FALSE)
  if (!is.character(name) || length(name) != 1L || is.na(name) || !nzchar(name))
    stop("name must be a nonempty string.", call. = FALSE)
  structure(list(name = name, param = param, functions = callbacks,
                 binding = "raw"), class = "frpp_residual_law")
}

fspec <- function(kernel, residual) {
  if (!methods::is(kernel, "hspec"))
    stop("kernel must be an emhawkes hspec object.", call. = FALSE)
  methods::validObject(kernel)
  if (any(!vapply(unname(residual_slots), function(s)
    is.null(methods::slot(kernel, s)), logical(1))))
    stop("kernel already contains a residual law; convert it with as_frpp().", call. = FALSE)
  if (!inherits(residual, "frpp_residual_law"))
    stop("residual must be created by residual_law().", call. = FALSE)
  structure(list(kernel = kernel, residual = residual,
    refresh = "all_dimensions_after_every_global_event", version = 1L),
    class = c("fspec", "frpp_spec"))
}

# Keep source compatibility with the first local prototype.
frpp_spec <- function(kernel, residual) fspec(kernel, residual)

# Preserve existing functions and their default bindings; do not inspect or
# reconstruct their enclosing environments or guess a distribution family.
as_frpp <- function(object, name = "legacy") {
  if (inherits(object, "frpp_spec")) {
    if (!inherits(object, "fspec")) class(object) <- c("fspec", class(object))
    return(object)
  }
  if (!methods::is(object, "hspec"))
    stop("object must be an hspec or frpp_spec object.", call. = FALSE)
  callbacks <- lapply(unname(residual_slots), function(s) methods::slot(object, s))
  names(callbacks) <- names(residual_slots)
  law <- residual_law(density = callbacks$density, cdf = callbacks$cdf,
    quantile = callbacks$quantile, rng = callbacks$rng, name = name)
  law$binding <- "bound"
  # A bound legacy law has function-specific defaults, not a claimed common
  # parameter vector. No information is inferred by running user code.
  law$param <- NULL
  kernel <- object
  for (s in unname(residual_slots)) methods::slot(kernel, s) <- NULL
  frpp_spec(kernel, law)
}

as_fspec <- function(object, name = "legacy") as_frpp(object, name)

# Transitional adapter. The resulting hspec is the old FRPP representation.
as_legacy_hspec <- function(object) {
  if (!inherits(object, "frpp_spec"))
    stop("object must be a frpp_spec.", call. = FALSE)
  law <- object$residual
  if (law$binding == "raw") {
    args <- law$functions
    names(args) <- unname(residual_slots)
    return(do.call(emhawkes::set_residual,
      c(list(object = object$kernel, param = law$param), args)))
  }
  kernel <- object$kernel
  for (key in names(residual_slots))
    methods::slot(kernel, residual_slots[[key]]) <- law$functions[[key]]
  kernel
}

print.frpp_spec <- function(x, ...) {
  cat("FRPP specification (fspec)\n",
      "Event types: ", x$kernel@dimens, "\n",
      "Residual law: ", x$residual$name, "\n",
      "Residual candidates refreshed after every global event.\n", sep = "")
  invisible(x)
}

print.fspec <- print.frpp_spec

require_law_functions <- function(object, needed) {
  absent <- needed[vapply(object$residual$functions[needed], is.null, logical(1))]
  if (length(absent)) stop("This operation requires residual function(s): ",
    paste(absent, collapse = ", "), ".", call. = FALSE)
}
