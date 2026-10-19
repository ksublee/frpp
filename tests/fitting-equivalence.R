library(frpp)
source(file.path("helpers", "kernel-api.R"))
# An identified small excitation fit; beta and residual shape/scale are fixed.
kernel <- methods::new("hspec", mu = .4, alpha = .2, beta = 1.2)
law <- residual_law(c(shape = 2, rate = 2), stats::dgamma, stats::pgamma,
  stats::qgamma, stats::rgamma)
model <- fspec(kernel, law)
legacy <- as_legacy_hspec(model)
# Fixed 65-row observations generated once by emhawkes 1.0.1 (0926075):
# seed 174, the kernel/law above, initial .1, N0=2, Nc0=3. Seventeen-digit
# decimals round-trip the original doubles. Do not regenerate from the code
# under test: a simulator change must not silently change this fitting input.
dt <- c(
  0, 3.7904061316193176, 1.6044939318316225,
  2.9827962340384748, 0.13515320933214961, 0.98914017778931795,
  3.9685676033834016, 5.1896487615632578, 2.8729921529914408,
  1.1659470959598766, 0.64289657871705475, 2.9229635525850957,
  2.025148313777049, 4.3816284923007585, 1.0570519824305626,
  1.2830307449358054, 0.24228866132642379, 0.74884656887378442,
  0.78919318262976013, 1.0247313918205558, 1.4990694774539342,
  0.44574425291112973, 0.84092078818385063, 6.1149354033847505,
  3.8322415804566416, 1.5377987200430558, 0.849199444325481,
  1.9927620440598304, 1.2966808150289164, 3.2310815884493724,
  3.0116538861584305, 1.6233520670084074, 0.8135295030671118,
  0.85739018914438114, 0.87158286611134106, 2.5796587837900282,
  2.1800243576095162, 1.8375943482758181, 1.4723009796824602,
  3.3755456845034488, 0.14236924627505909, 0.17203499053065774,
  1.9080675696396741, 1.3055051686064183, 0.73760511775011628,
  0.91301869617980913, 0.74935686034425397, 0.27476653328900436,
  1.7931518629571552, 0.91471299418933971, 3.3530176325674139,
  2.3906213377306242, 7.1429539194133884, 0.84011361127653061,
  3.4532651016015459, 1.5520323719547602, 2.3441379770977635,
  1.1978061040930725, 3.5986673211843323, 5.0177165380311468,
  1.2771579855980804, 4.1095874761758022, 1.1873101951905274,
  5.7781232599344872, 5.5978022213917953
)
path <- list(inter_arrival = dt, type = c(0, rep(1, length(dt) - 1L)))
set.seed(902)
rng <- .Random.seed
objective <- function(object, inter_arrival, type, mark, N, Nc, N0, Nc0,
                      lambda_component0, param) {
  stats::logLik(object, inter_arrival, type, mark, N = N, Nc = Nc,
    N0 = N0, Nc0 = Nc0, lambda_component0 = lambda_component0)
}
for (provider in c("stats", if (requireNamespace("maxLik", quietly = TRUE)) "maxLik")) {
  extra <- if (provider == "stats") list(method = "L-BFGS-B",
    lower = c(.02, 0, .1, .1, .1), upper = rep(10, 5),
    control = list(maxit = 200, factr = 1e7, pgtol = 1e-7)) else
    list(method = "BFGS", constraints = list(ineqA = rbind(c(1,0,0,0,0), c(0,1,0,0,0)),
      ineqB = c(-.02, 0)), control = list(iterlim = 200, reltol = 1e-10))
  args <- c(list(inter_arrival = path$inter_arrival, type = path$type,
    lambda_component0 = .1, N0 = 2, Nc0 = 3,
    fixed = c(FALSE, FALSE, TRUE, TRUE, TRUE), optimizer = provider), extra)
  a <- suppressMessages(do.call(emhawkes::hfit, c(list(object = legacy,
    mylogLik = objective, final_hessian = "optimizer"), args)))
  b <- suppressMessages(do.call(ffit, c(list(object = model, estimate_residual = TRUE), args)))
  H <- if (provider == "stats") a$hessian else maxLik::hessian(a)
  Hb <- if (provider == "stats") b$fit$hessian else maxLik::hessian(b$fit)
  stopifnot(identical(coef(a), coef(b)), identical(H, Hb),
    identical(as.numeric(logLik(a)), as.numeric(logLik(b))),
    identical(vcov(a), vcov(b)), identical(rng, .Random.seed),
    identical(emhawkes::fit_initial_state(a), finitial_state(b)),
    isTRUE(fdiagnostics(b)$optimizer_success),
    identical(fdiagnostics(b)$wald, "available_unvalidated"))
  # Restoration is a separate public call; the state backend uses the common
  # kernel API, rather than the compatibility likelihood. Keep its existing
  # arithmetic tolerance, including summation-order differences.
  args <- c(list(object = fitted_fspec(b), inter_arrival = dt,
    type = path$type, infer = TRUE), finitial_state(b))
  bridge <- do.call(stats::logLik, args)
  stopifnot(identical(as.numeric(logLik(b)), unname(bridge$loglikelihood)))
  if (frpp_has_kernel_api) {
    state <- do.call(stats::logLik, c(args, list(backend = "state")))
    stopifnot(isTRUE(all.equal(bridge, state[names(bridge)], tolerance = 1e-10)))
  }
  stopifnot(identical(rng, .Random.seed))
  cat(provider, ": exact coefficients, likelihood, raw Hessian, covariance, initial state and RNG; replay passed.\n")
}
