# From the frpp repository, after installing both CURRENT source trees:
# Rscript tools/validate-fast.R results/fast [path/to/emhawkes]
# Installed packages supply implementation; source trees supply tests.
args <- commandArgs(TRUE)
if (length(args) < 1L || length(args) > 2L)
  stop("Supply an output directory and optionally the emhawkes source directory.")
frpp_root <- normalizePath(".", mustWork = TRUE)
em_root <- normalizePath(if (length(args) == 2L) args[2L] else "../emhawkes", mustWork = TRUE)
for (entry in list(c(frpp_root, "frpp"), c(em_root, "emhawkes"))) {
  description <- read.dcf(file.path(entry[1L], "DESCRIPTION"))
  if (description[1L, "Package"] != entry[2L]) stop("Unexpected source package: ", entry[1L])
}
for (package in c("emhawkes", "frpp"))
  if (!requireNamespace(package, quietly = TRUE)) stop("Install the current ", package, " source first.")
if (!all(c("as_exp_kernel", "exp_kernel_replay") %in% getNamespaceExports("emhawkes")))
  stop("This joint development suite requires the emhawkes shared-kernel API. Install the development source first; no tests were skipped.")
out <- args[1L]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out, mustWork = TRUE)
# Discover every ordinary FRPP test: future regression files must not be missed.
frpp_tests <- sort(list.files(file.path(frpp_root, "tests"), "[.]R$", full.names = TRUE))
if (!length(frpp_tests)) stop("No FRPP tests found.")
em_tests <- file.path(em_root, "tests", c("exp-kernel.R", "simulate-interface.R", "fit-observations.R"))
tests <- c(frpp_tests, em_tests)
stopifnot(all(file.exists(tests)))
labels <- c(paste0("frpp/", basename(frpp_tests)), paste0("emhawkes/", basename(em_tests)))
libraries <- file.path(out, "libraries.rds")
saveRDS(.libPaths(), libraries)
wrapper <- file.path(out, "run-one.R")
writeLines(c("args <- commandArgs(TRUE)", ".libPaths(readRDS(args[1L]))",
  "sys.source(args[2L], envir = .GlobalEnv)"), wrapper)
executable <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
has_maxlik <- requireNamespace("maxLik", quietly = TRUE)
coverage <- c("Scope: all FRPP tests plus the three shared-interface emhawkes tests.",
  "Tests execute INSTALLED implementations. Reinstall both source trees after runtime edits; hashes record identity, not source/install equivalence.",
  paste0("maxLik: ", if (has_maxlik) "included" else "NOT INSTALLED; optional provider cases will not run"),
  "Not a replacement for full R CMD check, platform CI, large joint fits or coverage experiments.",
  "No fitting benchmark or performance claim: elapsed times include fresh R startup.")
writeLines(coverage, file.path(out, "coverage.txt"))
cat(paste(coverage, collapse = "\n"), "\n")
# Record both tested source and installed implementation, not just version labels.
source_hashes <- do.call(rbind, lapply(c(frpp_root, em_root), function(root) {
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
