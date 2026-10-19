# Shared test capability contract, not package runtime code.
frpp_has_kernel_api <- all(c("as_exp_kernel", "exp_kernel_replay") %in%
  getNamespaceExports("emhawkes"))
frpp_test_profile <- Sys.getenv("FRPP_TEST_PROFILE", "auto")
stopifnot(frpp_test_profile %in% c("auto", "release", "development"))
if (frpp_test_profile == "development" && !frpp_has_kernel_api)
  stop("Development validation requires the public emhawkes kernel API; no state tests may be omitted.")
if (identical(Sys.getenv("FRPP_TEST_REQUIRE_MAXLIK"), "true") &&
    !requireNamespace("maxLik", quietly = TRUE))
  stop("CI requires maxLik so both optimizer providers are tested.")
frpp_backends <- c("bridge", if (frpp_has_kernel_api) "state")
cat("emhawkes", as.character(utils::packageVersion("emhawkes")),
  "; profile:", frpp_test_profile, "; backends:", paste(frpp_backends, collapse = ", "), "\n")
if (!frpp_has_kernel_api)
  cat("State numerical comparisons unavailable in this emhawkes version; bridge checks and explicit capability-error checks remain active.\n")
