# Parameter registration and restoration, 2026-10-19

## Implemented boundary

emhawkes exports `parameter_plan(components, prefixes, named, reduced)`. It
registers model-independent numeric matrices/scalars, explicit named vectors,
functions and NULL components. A locked plan exposes the full starting vector
and component restoration. Matrix naming, reduced/full sharing, global name
sharing, conflicting-start warnings, function environments and default binding
follow the existing rules. Parameter defaults are evaluated in order; user
default expressions can have side effects. Nonfinite proposals remain for the
calling objective to handle. Optimizers, constraints and initial states are
outside the plan. Existing hspec map metadata keeps its original representation.

For numeric kernels without event callbacks and raw residual laws, frpp declares
kernel matrices and named residual coordinates, then assembles its fitted model
from the public plan. Fixed residual laws remain fixed. The explicit FRPP plan
is prepared after optimization, preserving existing argument-forcing and warning
order. Existing sharing warnings are emitted once by hfit.

Optimization and objective preparation still use public hfit. Callback/bound-law
models and older emhawkes versions retain fitted_hspec reconstruction. Raw
optimizer results and the ffit result fields remain available. Simulation's
bridge default and likelihood-engine options are unchanged.

## Exact comparisons

- Baseline installations: emhawkes ac4dc8c and frpp 8d0fdc1, Windows R 4.5.1.
- Twenty-six old/new parameter maps: starting values, names/order, component
  indices, keys, dimensions and restored values exactly equal. Includes reduced
  and full coordinates, non-square kernels, eta, and shared callback parameters.
- Forty simulation/likelihood/RNG reference records and two fixed-65-row
  optimizer records exactly equal. Coefficients, likelihood, raw Hessian,
  covariance, initial state and diagnostics match for stats and maxLik.
- A fresh 10,000-row hfit/ffit pair was re-estimated with matched maxLik controls.
  Coefficients, likelihood, raw Hessian and initial state match exactly between
  interfaces. Each result also matches all three corresponding earlier fits:
  six historical comparisons, including covariance, code/message, iterations
  and score. Both fresh fits used 338 validator calls and 55 iterations.
- All six historical large fits replay successfully, retaining the existing
  relative 1e-10 state/bridge arithmetic criterion and unchanged RNG. The known
  absolute likelihood sum-order difference is 1.12777343019843e-10; no tolerance
  was increased. Fresh paired fitting comparisons use exact equality.
- New tests prove numeric restoration does not call fitted_hspec, callback
  models retain that fallback, and colliding kernel/residual starts warn once.
- New stats and maxLik fits saved with this contract can be read under the
  previous emhawkes installation without parameter_plan. Coefficients,
  likelihood, covariance, initialization, diagnostics, subsequent simulation
  and RNG are exactly reproduced.

## Local checks

| Check | Result |
|---|---|
| Quick development suite | 15/15 files, 21.81 s including fresh R startup |
| Previous development dependency without parameter API | 11/11 FRPP files |
| CRAN emhawkes 1.0.1 | 11/11 FRPP files |
| emhawkes R CMD check | 50 test files and vignette rebuild; Status: OK |
| frpp with development emhawkes | 11 test files and both vignettes; Status: OK |
| frpp with CRAN emhawkes | 11 test files and both vignettes; Status: OK |
| Updated computation guides | Both rendered with code examples |

R CMD check uses --no-manual --run-donttest and
_R_CHECK_SUGGESTS_ONLY_=false. The PDF reference manual was not compiled.
Unrelated emhawkes working-tree edits were preserved and excluded from the
tested source snapshot. No commit, push or new platform CI is claimed here.

## Limits and validation-harness correction

This is the registration/restoration stage. It does not replace optimizer
adapters or legacy objective preparation. Callback migration remains separate.
The large example retains unavailable Wald covariance for nonpositive curvature;
exact reproduction is not a confidence-interval coverage claim.

The fresh large refit was an accuracy check, one pair rather than a timing
benchmark. Its elapsed values were deliberately omitted. The original optional
timing footer failed on all-NA values after both fits and their comparisons had
completed. The footer was corrected and checked against the actual saved
records; strict numerical comparisons and result finalization passed. No fitting
speed claim is made.
