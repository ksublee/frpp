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

frpp_has_sim_api <- all(c("exp_kernel_step", "exp_kernel_wait") %in%
  getNamespaceExports("emhawkes"))
if (frpp_test_profile == "development" && !frpp_has_sim_api)
  stop("Development validation requires the event-wise emhawkes kernel API; independent simulation tests may not be omitted.")
frpp_has_plan_api <- "exp_kernel_plan" %in% getNamespaceExports("emhawkes")
if (frpp_test_profile == "development" && !frpp_has_plan_api)
  stop("Development validation requires the prepared emhawkes kernel API; no plan tests may be omitted.")
cat("Prepared kernel plan:", if (frpp_has_plan_api) "enabled" else "unavailable; existing API compatibility enabled", "\n")
cat("Independent simulation:", if (frpp_has_sim_api) "enabled" else "unavailable; bridge compatibility enabled", "\n")
