# Run after installing the prototype and emhawkes 1.0.0 in the same library path.
# Usage: Rscript --vanilla tools/benchmark.R OUTPUT_DIRECTORY
library(frpp)
args <- commandArgs(TRUE)
if (length(args) != 1L) stop("Supply one output directory.")
dir.create(args[1], recursive = TRUE, showWarnings = FALSE)
options(emhawkes.engine = "auto")
kernel <- methods::new("hspec", mu = c(.3, .4), alpha = matrix(.025, 2, 4),
  beta = matrix(rep(c(1.2, 1.7), 4), 2),
  type_col_map = list(c(1L, 3L), c(2L, 4L)))
model <- frpp_spec(kernel, residual_law(c(shape = 2, rate = 2),
  stats::dgamma, stats::pgamma, stats::qgamma, stats::rgamma))
legacy <- as_legacy_hspec(model)
rows <- list()
for (size in c(100L, 2000L)) {
  path <- stats::simulate(model, size = size, seed = 73, lambda_component0 = .08)
  data <- list(inter_arrival = path$inter_arrival, type = path$type,
    mark = path$mark, lambda_component0 = .08)
  calls <- list(
    legacy = function() do.call(emhawkes::logLik, c(list(object = legacy), data)),
    bridge = function() do.call(stats::logLik, c(list(object = model), data)),
    state = function() do.call(stats::logLik,
      c(list(object = model, backend = "state"), data)))
  values <- vapply(calls, function(f) unname(f()), numeric(1))
  stopifnot(identical(values[["legacy"]], values[["bridge"]]),
    abs(values[["state"]] - values[["legacy"]]) < 1e-10)
  times <- matrix(NA_real_, 3, 3, dimnames = list(NULL, names(calls)))
  # Warm calls above are excluded. Alternate the execution order; no compilation
  # or other package checks should run concurrently with this measurement.
  for (repeat_id in 1:3) {
    order <- if (repeat_id %% 2L) names(calls) else rev(names(calls))
    for (name in order) times[repeat_id, name] <- system.time(calls[[name]]())[["elapsed"]]
  }
  for (name in names(calls)) {
    profile <- file.path(args[1], paste0(name, "-", size, "-allocations.out"))
    gc()
    Rprofmem(profile)
    tryCatch(calls[[name]](), finally = Rprofmem(NULL))
    allocations <- suppressWarnings(as.numeric(sub(" .*", "", readLines(profile))))
    rows[[length(rows) + 1L]] <- data.frame(size = size, backend = name,
      median_seconds = median(times[, name]),
      cumulative_allocation_bytes = sum(allocations, na.rm = TRUE),
      absolute_likelihood_difference = abs(values[[name]] - values[["legacy"]]))
  }
}
results <- do.call(rbind, rows)
write.csv(results, file.path(args[1], "benchmark.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(args[1], "session-info.txt"))
print(results, row.names = FALSE)
