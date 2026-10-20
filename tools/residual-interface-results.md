# Residual-function names and optional inputs, 2026-10-20

## Interface

residual_law uses dresidual, presidual, qresidual and rresidual as its primary
function arguments. density, cdf, quantile and rng remain nondeprecated aliases.
The first six positional arguments retain their order. Different roles may use
different spellings. Two non-NULL functions for one role are rejected even if
identical; an explicit NULL does not conflict with the other spelling.

Stored function keys remain density/cdf/quantile/rng for compatibility with
saved laws and existing accessors. Functions, environments and lazy defaults are
retained unchanged. Neither construction nor validation invokes the callbacks
or creates a random seed. Numerical kernels, optimizers and backend defaults
are unchanged.

| Operation | Required functions |
|---|---|
| Construction | At least dresidual or rresidual |
| Simulation | rresidual |
| Univariate likelihood/fitting | dresidual |
| Multivariate likelihood/fitting | dresidual, presidual |
| infer = TRUE | dresidual, presidual, qresidual |

Missing-function errors name the primary argument and its alias. ffit now checks
its density/CDF requirements before optimization or user validator execution.
No omitted callback is inferred or synthesized. In particular, CDF/quantile-only
laws cannot be constructed; no current modeling operation can use those alone.

## Validation

- All 16 function-presence subsets: 12 constructible laws and four constructor
  rejections, with the expected success/error for each operation in 1D/2D and
  bridge/state backends when available.
- All 16 ways to mix primary and alias spellings, positional calls, conflict
  errors, explicit NULLs, invalid types, environment/default preservation and
  seed absence/preservation.
- Generator-only simulation exactly matches full-law path fields, final RNG
  state and subsequent draws. Density-only univariate fitting matches full-law
  coefficients, likelihood, covariance, initialization and diagnostics with
  both maxLik/stats and bridge/independent fitting when supported.
- An additional 32 partial-law fits estimating the residual rate (rather than
  freezing it) exactly match full-law fits across both providers and fit routes,
  including restored residual parameters and unchanged RNG.
- Missing density/CDF rejects fitting before a deliberately failing validator
  can run, with unchanged RNG. infer succeeds without a generator.
- Twelve actual laws and 48 records captured from the previous a453ea1 runtime
  exactly match new primary-name construction: saved law representation and
  applicable paths, likelihoods, inferred states/residuals, RNG and following
  draws. Existing positional and alias fixtures remain in the ordinary suite.
- Full R CMD check --no-manual --run-donttest with
  _R_CHECK_SUGGESTS_ONLY_=false passes against development and CRAN emhawkes,
  including all 13 FRPP test files and both vignette rebuilds.
- Quick development: 17 files. Quick release: 13 files. Full pkgdown site builds
  in a fresh default Rscript session without pre-attaching frpp.

README, function help, the examples vignette and website guides now use the
primary names and show generator-only simulation and density-only fitting.
No simulation/fitting speed or confidence-interval coverage claim is made.
