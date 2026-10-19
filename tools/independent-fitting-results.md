# Independent numeric fitting, 2026-10-19

## Boundary and defaults

`ffit(..., fit_backend = "independent")` prepares numeric observations and
initial states, registers kernel/residual coordinates through the public
emhawkes parameter plan, calls maxLik or stats directly and restores an fspec.
It does not call hfit or emhawkes internal functions. The local stats adapter
is adapted from emhawkes 18986df (GPL >= 2), retaining its fixed/free reduction,
constraint and box handling, numerical-difference arithmetic and covariance
rules. Future optimizer fixes need coordinated regression checks in both packages.

Both fitting execution and likelihood backend default to bridge. `fit_backend`
and likelihood `backend` are independent choices. Independent fitting supports
plain numeric mu/alpha/beta/eta, mappings and raw residual laws with a density;
the selected likelihood retains its own restrictions. Event callbacks, functional
kernels, legacy bound laws and older dependencies retain bridge fitting. Auto
selects support before optimization, records a reason and never retries a failed
fit on another path. Simulation defaults and RNG rules are unchanged.

Full parameter names/order, reduced/full sharing, conflicting-start warning and
first-value semantics, residual freeze/estimate behavior, fixed selectors,
constraints, optional derivatives, validator rejection/errors and explicit
initial states retain their existing meaning. New results add a fitting-route
field. Existing numerical diagnostics are unchanged. maxLik keeps its provider
class; independent stats fits use frpp_stats_fit with standard accessors. Saved
ffit accessors remain available; older saved fits need not contain route metadata.

## Verification

- Fixed 65-row observations: exact coefficients, likelihood, raw Hessian,
  covariance, initialization and RNG against actual hfit fits, both providers.
- Small paired fits: BFGS/CG/Nelder-Mead/L-BFGS-B, linear inequalities, box
  boundaries, fixed/estimated residuals, conflicting shared names, reduced/full
  mapped kernels, eta and observed marks. Numeric results and diagnostics match.
- Validator candidate sequences and RNG state match exactly. Unsupported
  callback requests error explicitly; auto retains the bridge and its reason.
- maxLik SANN also matches coefficients, curvature and following random draws;
  its existing global-seed restoration behavior is preserved in both paths.
- A traced hfit that errors on entry is never reached by independent fitting.
- Explicit gradient/Hessian checks pass for both likelihood backends and both
  providers. Disabled final curvature and named fixed parameters match.
- Saved independent stats ffit is restored in a fresh minimal session, without
  attaching emhawkes/methods or loading maxLik. Earlier saved results still work.
- Forty original simulation/likelihood/RNG records and two default fitting
  records exactly reproduce the pre-contract baseline.
- One actual independent 10,000-row refit matches eight previously computed
  hfit/ffit results exactly: coefficients, likelihood, raw Hessian/covariance,
  initial state, score, return code and iterations. It uses the same 338 validator
  calls and 55 iterations. Restored likelihood and final RNG match exactly.
  This is an accuracy check, not a timing benchmark or performance claim.
- Development and released emhawkes dependencies both pass R CMD check
  --no-manual --run-donttest with _R_CHECK_SUGGESTS_ONLY_=false, including all
  12 FRPP test files and both vignette rebuilds. Quick development: 16 files;
  quick release: 12 files. The computation guide renders successfully.

No tolerances were increased. State/bridge likelihood summation differences keep
the existing criterion; same-route fitting comparisons use exact equality.
The large tzexp fit still has unavailable Wald covariance for nonpositive
curvature. Reproduction does not establish identification or interval coverage.

## Corrections recorded during validation

The first released-dependency check caught missing fallback metadata: its
unsupported-route reason was overwritten by the separate Wald diagnostic
variable. Distinct variables now retain both reasons. Exact numerical results
were unaffected; final full checks use the corrected source. The large refit
preceded this metadata-only correction; optimization/objective code is unchanged.

An additional validation script initially guessed scalar fixed-parameter names
using matrix-style names. It now reads names from the actual registered result;
both providers' extra comparisons pass. Original failed logs are retained in
the external validation artifact directory.
