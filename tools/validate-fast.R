# From the frpp repository, after installing both CURRENT source trees:
# Rscript tools/validate-fast.R results/fast [path/to/emhawkes] [--extended]
# Rscript tools/validate-fast.R results/release --release [--require-maxlik]
# Installed packages supply implementation; source trees supply tests.
args <- commandArgs(TRUE)
extended <- "--extended" %in% args
release <- "--release" %in% args
require_maxlik <- "--require-maxlik" %in% args
args <- args[!args %in% c("--extended", "--release", "--require-maxlik")]
if (length(args) < 1L || length(args) > 2L || (release && extended))
  stop("Supply an output directory, optionally an emhawkes source directory; --release and --extended cannot be combined.")
frpp_root <- normalizePath(".", mustWork = TRUE)
em_root <- if (release) NULL else normalizePath(if (length(args) == 2L)
  args[2L] else "../emhawkes", mustWork = TRUE)
roots <- c(frpp = frpp_root, if (!is.null(em_root)) c(emhawkes = em_root))
for (package in names(roots)) {
  description <- read.dcf(file.path(roots[[package]], "DESCRIPTION"))
  if (description[1L, "Package"] != package) stop("Unexpected source package: ", roots[[package]])
}
for (package in c("emhawkes", "frpp"))
  if (!requireNamespace(package, quietly = TRUE)) stop("Install the current ", package, " source first.")
has_api <- all(c("as_exp_kernel", "exp_kernel_replay") %in% getNamespaceExports("emhawkes"))
has_sim_api <- all(c("exp_kernel_step", "exp_kernel_wait") %in% getNamespaceExports("emhawkes"))
if (!release && !has_api)
  stop("This joint development suite requires the emhawkes shared-kernel API. Install the development source first; no tests were skipped.")
if (!release && !has_sim_api)
  stop("This joint development suite requires the event-wise emhawkes kernel API. No independent simulation checks were skipped.")
if (!release && !"exp_kernel_plan" %in% getNamespaceExports("emhawkes"))
  stop("This joint development suite requires the prepared emhawkes kernel API. No plan checks were skipped.")
if (require_maxlik && !requireNamespace("maxLik", quietly = TRUE))
  stop("Both optimizer providers were requested but maxLik is not installed.")
Sys.setenv(FRPP_TEST_PROFILE = if (release) "release" else "development")
if (require_maxlik) Sys.setenv(FRPP_TEST_REQUIRE_MAXLIK = "true")
out <- args[1L]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out, mustWork = TRUE)
# Discover every ordinary FRPP test: future regression files must not be missed.
frpp_tests <- sort(list.files(file.path(frpp_root, "tests"), "[.]R$", full.names = TRUE))
if (!length(frpp_tests)) stop("No FRPP tests found.")
em_tests <- if (release) character() else file.path(em_root, "tests", c("exp-kernel.R", "exp-kernel-step.R", "exp-kernel-plan.R",
  if (extended) c("simulate-interface.R", "fit-observations.R")))
tests <- c(frpp_tests, em_tests)
missing_tests <- tests[!file.exists(tests)]
if (length(missing_tests)) stop("Selected test files are missing: ", paste(missing_tests, collapse = ", "))
labels <- c(paste0("frpp/", basename(frpp_tests)),
  if (length(em_tests)) paste0("emhawkes/", basename(em_tests)))
libraries <- file.path(out, "libraries.rds")
saveRDS(.libPaths(), libraries)
wrapper <- file.path(out, "run-one.R")
writeLines(c("args <- commandArgs(TRUE)", ".libPaths(readRDS(args[1L]))",
  "sys.source(args[2L], envir = .GlobalEnv, chdir = TRUE)"), wrapper)
executable <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
has_maxlik <- requireNamespace("maxLik", quietly = TRUE)
coverage <- c(if (release)
  "Scope: all FRPP tests against the installed released emhawkes; no development source tests." else if (extended)
  "Scope: all FRPP tests, the shared kernel, and the additional simulate/observation-count development tests." else
  "Scope: all FRPP tests plus the public replay and event-wise emhawkes kernel tests.",
  paste0("State numerical comparisons: ", if (has_api) "required and enabled" else "unavailable; bridge and capability-error checks enabled"),
  paste0("Independent simulation: ", if (has_sim_api) "required and enabled" else "unavailable; bridge/auto fallback and capability-error checks enabled"),
  "Tests execute INSTALLED implementations. Reinstall both source trees after runtime edits; hashes record identity, not source/install equivalence.",
  paste0("maxLik: ", if (has_maxlik) "included" else "NOT INSTALLED; optional provider cases will not run"),
  "Not a replacement for full R CMD check, platform CI, large joint fits or coverage experiments.",
  "No fitting benchmark or performance claim: elapsed times include fresh R startup.")
writeLines(coverage, file.path(out, "coverage.txt"))
cat(paste(coverage, collapse = "\n"), "\n")
# Record both tested source and installed implementation, not just version labels.
source_hashes <- do.call(rbind, lapply(unname(roots), function(root) {
  paths <- c(file.path(root, c("DESCRIPTION", "NAMESPACE")),
    list.files(file.path(root, "R"), full.names = TRUE, recursive = TRUE),
    list.files(file.path(root, "src"), "[.](cpp|h)$", full.names = TRUE, recursive = TRUE),
    list.files(file.path(root, "tests"), "[.]R$", full.names = TRUE, recursive = TRUE))
  data.frame(path = paths, md5 = unname(tools::md5sum(paths)))
}))
write.csv(source_hashes, file.path(out, "source-md5.csv"), row.names = FALSE)
installed_hashes <- do.call(rbind, lapply(c("emhawkes", "frpp"), function(package) {
  root <- find.package(package)
  files <- list.files(root, full.names = TRUE, recursive = TRUE)
  files <- files[grepl("([.](rdb|rdx|dll|so|dylib)$|/DESCRIPTION$)", gsub("\\", "/", files, fixed = TRUE))]
  data.frame(package = package, version = as.character(packageVersion(package)),
    path = files, md5 = unname(tools::md5sum(files)))
}))
write.csv(installed_hashes, file.path(out, "installed-md5.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(out, "session-info.txt"))
rows <- list()
started <- proc.time()[["elapsed"]]
for (i in seq_along(tests)) {
  log <- file.path(out, paste0(gsub("/", "-", labels[i], fixed = TRUE), ".log"))
  cat("Running ", labels[i], " ... ", sep = ""); flush.console()
  before <- proc.time()[["elapsed"]]
  status <- system2(executable, c("--vanilla", shQuote(wrapper), shQuote(libraries),
    shQuote(tests[i])), stdout = log, stderr = log)
  elapsed <- proc.time()[["elapsed"]] - before
  rows[[i]] <- data.frame(test = labels[i], status = if (status == 0L) "PASS" else "FAIL",
    exit_code = status, elapsed_seconds = elapsed, log = log)
  write.csv(do.call(rbind, rows), file.path(out, "checks.csv"), row.names = FALSE)
  cat(rows[[i]]$status, " (", format(elapsed, digits = 4), " s)\n", sep = "")
  if (status != 0L) cat(paste(readLines(log, warn = FALSE), collapse = "\n"), "\n")
}
result <- do.call(rbind, rows)
total <- proc.time()[["elapsed"]] - started
writeLines(paste("Total test-loop seconds:", total), file.path(out, "elapsed.txt"))
if (any(result$status != "PASS")) stop("Validation failed; see checks.csv and individual logs.")
cat("PASS: ", nrow(result), " files in ", format(total, digits = 4), " seconds. Reports: ", out, "\n", sep = "")
