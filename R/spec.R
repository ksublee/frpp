residual_slots <- c(density = "dresidual", cdf = "presidual",
                    quantile = "qresidual", rng = "rresidual")

# Creating a law never evaluates its callbacks or consumes random numbers.
residual_law <- function(param = numeric(), dresidual = NULL, presidual = NULL,
                         qresidual = NULL, rresidual = NULL, name = "custom",
                         density = NULL, cdf = NULL, quantile = NULL, rng = NULL) {
  if (!is.numeric(param) || is.complex(param) || !is.null(dim(param)) ||
      any(!is.finite(param)) || (length(param) &&
      (is.null(names(param)) || anyNA(names(param)) || any(names(param) == "") ||
       anyDuplicated(names(param)))))
    stop("param must be a finite numeric vector with unique, nonempty names.", call. = FALSE)
  resolve <- function(primary, alias, primary_name, alias_name) {
    if (!is.null(primary) && !is.null(alias))
      stop("Supply only one of '", primary_name, "' or '", alias_name, "'.", call. = FALSE)
    if (is.null(primary)) alias else primary
  }
  density <- resolve(dresidual, density, "dresidual", "density")
  cdf <- resolve(presidual, cdf, "presidual", "cdf")
  quantile <- resolve(qresidual, quantile, "qresidual", "quantile")
  rng <- resolve(rresidual, rng, "rresidual", "rng")
  # Keep stored keys and callback identity for saved laws and existing accessors.
  callbacks <- list(density = density, cdf = cdf, quantile = quantile, rng = rng)
  if (!all(vapply(callbacks, function(x) is.null(x) || is.function(x), logical(1))))
    stop("Residual functions must be functions or NULL.", call. = FALSE)
  if (is.null(density) && is.null(rng))
    stop("Supply dresidual (density) or rresidual (rng) for the residual law.", call. = FALSE)
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

# Residual laws retain their existing S3 representation. Register their type
# for the containing S4 slot without changing callback environments.
methods::setOldClass("frpp_residual_law")

validate_residual_law <- function(law) {
  tryCatch({
    if (!inherits(law, "frpp_residual_law") || !is.list(law))
      stop("residual must be created by residual_law().")
    if (!is.character(law$binding) || length(law$binding) != 1L ||
        is.na(law$binding) || !law$binding %in% c("raw", "bound"))
      stop("Residual binding must be 'raw' or 'bound'.")
    if (!is.list(law$functions) ||
        !identical(names(law$functions), names(residual_slots)))
      stop("Residual functions must contain density, cdf, quantile and rng in that order.")
    if (law$binding == "bound" && !is.null(law$param))
      stop("A bound legacy residual law must have param = NULL.")
    # Reuse constructor validation; it inspects formals but runs no callbacks.
    invisible(do.call(residual_law, c(list(param = if (law$binding == "raw")
      law$param else numeric(), name = law$name), law$functions)))
    TRUE
  }, error = function(e) conditionMessage(e))
}

# Check stored numeric structure without evaluating functional parameters.
# hspec initialization checks values, but validObject(hspec) alone only checks
# declared slot types after direct slot edits.
validate_fspec_kernel <- function(kernel) {
  d <- kernel@dimens
  if (length(d) != 1L || !is.finite(d) || d < 1 || d != floor(d))
    return("kernel dimens must be a positive integer.")
  errors <- character(); columns <- integer()
  for (name in c("mu", "alpha", "beta", "eta")) {
    value <- methods::slot(kernel, name)
    if (is.function(value) || (name == "eta" && is.null(value))) next
    if (!is.matrix(value) || !is.numeric(value) || is.complex(value) ||
        any(!is.finite(value)) || nrow(value) != d || ncol(value) < 1L) {
      errors <- c(errors, paste0("kernel ", name, " must be a finite numeric matrix with dimens rows."))
      next
    }
    if (name == "mu") {
      if (ncol(value) != 1L || any(value < 0))
        errors <- c(errors, "kernel mu must have one nonnegative column.")
    } else columns <- c(columns, ncol(value))
    if (name == "beta" && any(value <= 0)) errors <- c(errors, "kernel beta must be positive.")
  }
  columns <- unique(columns)
  if (length(columns) > 1L) errors <- c(errors, "kernel alpha, beta and eta column counts must match.")
  mapping <- kernel@type_col_map
  if (is.null(mapping)) {
    if (length(columns) == 1L && columns != d)
      errors <- c(errors, "kernel type_col_map is required for a non-square kernel.")
  } else if (length(mapping) != d || !all(vapply(mapping, function(x)
    is.numeric(x) && !is.complex(x) && length(x) > 0L && all(is.finite(x)) &&
    all(x == floor(x) & x >= 1) && !anyDuplicated(x) &&
    (length(columns) != 1L || all(x <= columns)), logical(1)))) {
    errors <- c(errors, "kernel type_col_map must contain valid distinct column indices for each event type.")
  }
  if (length(errors)) errors else TRUE
}

methods::setClass("fspec",
  slots = c(kernel = "hspec", residual = "frpp_residual_law",
    refresh = "character", version = "integer"),
  prototype = list(refresh = "all_dimensions_after_every_global_event", version = 2L),
  validity = function(object) {
    errors <- character()
    kernel_valid <- methods::validObject(object@kernel, test = TRUE)
    if (!isTRUE(kernel_valid)) errors <- c(errors, kernel_valid)
    numeric_valid <- validate_fspec_kernel(object@kernel)
    if (!isTRUE(numeric_valid)) errors <- c(errors, numeric_valid)
    if (any(!vapply(unname(residual_slots), function(s)
        is.null(methods::slot(object@kernel, s)), logical(1))))
      errors <- c(errors, "kernel already contains a residual law; convert it with as_frpp().")
    law_valid <- validate_residual_law(object@residual)
    if (!isTRUE(law_valid)) errors <- c(errors, law_valid)
    if (!identical(object@refresh, "all_dimensions_after_every_global_event"))
      errors <- c(errors, "Only refresh after every global event is supported.")
    if (!identical(object@version, 2L)) errors <- c(errors, "Unsupported fspec version; expected 2L.")
    if (length(errors)) errors else TRUE
  })

fspec <- function(kernel, residual) {
  if (!methods::is(kernel, "hspec"))
    stop("kernel must be an emhawkes hspec object.", call. = FALSE)
  if (!inherits(residual, "frpp_residual_law"))
    stop("residual must be created by residual_law().", call. = FALSE)
  methods::new("fspec", kernel = kernel, residual = residual)
}

frpp_spec <- function(kernel, residual) fspec(kernel, residual)

is_fspec <- function(object) (isS4(object) && methods::is(object, "fspec")) ||
  (!isS4(object) && is.list(object) && inherits(object, c("fspec", "frpp_spec")))

# Preserve callback identity/defaults; never infer a family by executing code.
as_frpp <- function(object, name = "legacy") {
  if (isS4(object) && methods::is(object, "fspec")) {
    methods::validObject(object)
    return(object)
  }
  if (is_fspec(object)) {
    if (!is.null(object$version) && !identical(object$version, 1L))
      stop("Unsupported legacy fspec version; expected 1L.", call. = FALSE)
    if (!is.null(object$refresh) &&
        !identical(object$refresh, "all_dimensions_after_every_global_event"))
      stop("Only refresh after every global event is supported.", call. = FALSE)
    return(fspec(object$kernel, object$residual))
  }
  if (!methods::is(object, "hspec"))
    stop("object must be an fspec, legacy frpp_spec or hspec object.", call. = FALSE)
  callbacks <- lapply(unname(residual_slots), function(s) methods::slot(object, s))
  names(callbacks) <- names(residual_slots)
  law <- residual_law(density = callbacks$density, cdf = callbacks$cdf,
    quantile = callbacks$quantile, rng = callbacks$rng, name = name)
  law$binding <- "bound"
  law$param <- NULL
  kernel <- object
  for (s in unname(residual_slots)) methods::slot(kernel, s) <- NULL
  fspec(kernel, law)
}

as_fspec <- function(object, name = "legacy") as_frpp(object, name)

# Accessors are the preferred boundary for future kernel representation changes.
fkernel <- function(object) as_fspec(object)@kernel
fresidual <- function(object) as_fspec(object)@residual
`fkernel<-` <- function(object, value) {
  object <- as_fspec(object)
  object@kernel <- value
  methods::validObject(object)
  object
}
`fresidual<-` <- function(object, value) {
  object <- as_fspec(object)
  object@residual <- value
  methods::validObject(object)
  object
}

# Transitional list-style access remains available, including nested updates.
methods::setMethod("$", "fspec", function(x, name) {
  if (name %in% methods::slotNames(x)) methods::slot(x, name) else NULL
})
methods::setReplaceMethod("$", "fspec", function(x, name, value) {
  if (!name %in% methods::slotNames(x)) stop("Unknown fspec field: ", name, ".", call. = FALSE)
  methods::slot(x, name) <- value
  methods::validObject(x)
  x
})

as_legacy_hspec <- function(object) {
  if (!is_fspec(object)) stop("object must be an fspec or legacy frpp_spec.", call. = FALSE)
  object <- as_fspec(object)
  law <- object@residual
  if (law$binding == "raw") {
    args <- law$functions
    names(args) <- unname(residual_slots)
    return(do.call(emhawkes::set_residual,
      c(list(object = object@kernel, param = law$param), args)))
  }
  kernel <- object@kernel
  for (key in names(residual_slots))
    methods::slot(kernel, residual_slots[[key]]) <- law$functions[[key]]
  kernel
}

print.frpp_spec <- function(x, ...) {
  x <- as_fspec(x)
  cat("FRPP specification (S4 fspec)\n",
      "Event types: ", x@kernel@dimens, "\n",
      "Residual law: ", x@residual$name, "\n",
      "Residual candidates refreshed after every global event.\n", sep = "")
  invisible(x)
}
print.fspec <- print.frpp_spec
methods::setMethod("show", "fspec", function(object) { print.fspec(object); invisible(NULL) })

require_law_functions <- function(object, needed) {
  absent <- needed[vapply(object@residual$functions[needed], is.null, logical(1))]
  if (length(absent)) stop("This operation requires residual function(s): ",
    paste(paste0(unname(residual_slots[absent]), " (", absent, ")"), collapse = ", "),
    ".", call. = FALSE)
}
