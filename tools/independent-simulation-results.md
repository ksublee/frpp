# Independent numeric FRPP simulation, 2026-10-19

The initial implementation is recorded below. See the [prepared-plan follow-up](#prepared-plan-follow-up-2026-10-19) for current development behavior and paired timings.

## Implementation and support

`fsim(..., backend = "state")` now runs its own R event loop. It uses public
emhawkes `exp_kernel_wait()` and `exp_kernel_step()` for deterministic inversion,
decay, interval integration and mapped updates. emhawkes owns the common
mathematics; frpp owns residual draws, event selection, counts and output.
The package dependency remains frpp -> emhawkes; this is not a standalone engine.
Existing Hawkes and legacy FRPP loops and likelihood-engine options are unchanged.

Support: plain numeric mu/alpha/beta/eta, row-wise equal beta, multivariate and
non-square/mapped kernels, overlapping mappings, nonnegative baselines including
zero, and admissible signed components. Marks are one. Baseline, impact,
mark-generator and mark-density event callbacks remain in the bridge. Residual
generators retain their existing public parameter binding.

Bridge remains the default because state is slower on larger paths. Explicit
`auto` selects state when supported, otherwise bridge. Selection precedes any
seed change or draw; an error after a draw never retries another backend.
Inputs in `...` are not forced during selection. Initial/count expressions retain
their old evaluation order inside the seed scope, and the RNG generator receives
the original dimension argument type. Legacy initial-state recycling and its
warning are preserved. A supplied seed restores kind/state even on failure.

## Accuracy and local checks

- New daily fixture: four kernels x three residual laws x two seeds x sizes
  1, 2 and 17 = 72 combinations. Full paths, final RNG and following draws are
  identical to bridge. All four driver/component outputs independently restored
  by existing R and C++ likelihood engines are identical too.
- Existing baseline: 40 path/inference/RNG records and two optimizer records
  exactly match emhawkes cde77c9 / frpp 91c0572. The new capture explicitly uses
  state for numeric cases and bridge for event callbacks. Coefficients, raw
  Hessian, covariance, likelihood, initial state and RNG match for both providers.
- Actual vignette examples: 3,000-row univariate tzexp, 10,000-row bivariate
  tzexp and 5,000-row uniform. Development uses state explicitly; every common
  simulation field and final RNG matches the original examples exactly.
- A separate namespace trace rejects hsim and replay calls: 257 rows produce
  257 step calls (including initialization) and 256 wait calls, with exact output.
- Ties, simultaneous events, initial-only paths, finite total hazard/early stop,
  zero baseline, invalid post-event state before the next draw, lazy arguments,
  callbacks and fresh-session seed restoration are covered. A development run
  with replay but without event-wise APIs fails before testing.
- Existing numerical tolerances are unchanged. New deterministic step comparisons
  are exact; inversion quadrature uses 1e-11 with the existing 1e-12 root tolerance.
  State likelihood comparisons retain the existing relative 1e-10 criterion for
  arithmetic summation-order differences; this does not relax path comparisons.

Local R 4.5.1, Windows 11, maxLik 1.5-2.1:

| Check | Result |
|---|---|
| emhawkes R CMD check, 48 test files and vignette rebuild | Status: OK |
| frpp with development emhawkes, nine test files and both vignettes | Status: OK |
| frpp with CRAN emhawkes 1.0.1, nine files and both vignettes | Status: OK |
| Final quick development suite | 11/11, 20.13 s full command |
| Final quick release compatibility suite | 9/9, 12.67 s full command |
| Shared-kernel website guide, including the new API example | knit passed |

Full checks used `--no-manual --run-donttest` and `_R_CHECK_SUGGESTS_ONLY_=false`.
PDF manuals were not compiled. Tests use installed packages; the built source
was compared with current runtime/test files. Source/installed hashes and logs
are retained by the quick runner. GitHub workflows supply platform results for
these commits; the table above records local checks.

## Timing and memory

Whole fsim calls, excluding installation/compilation; one full warm-up per
backend, three repetitions in alternating order, same model/seed/initial state.
The baseline is the already optimized legacy R FRPP path from cde77c9, not an
old pre-optimization release. Every timed output matches exactly. No other
heavy local validation was run concurrently. Values are medians in seconds.

| Model | Rows, including initial | Bridge | State |
|---|---:|---:|---:|
| univariate | 32 | <0.01 | 0.01 |
| univariate | 1,000 | 0.10 | 0.15 |
| univariate | 10,000 | 1.25 | 1.90 |
| mapped | 32 | <0.01 | 0.02 |
| mapped | 1,000 | 0.25 | 0.33 |
| mapped | 10,000 | 2.61 | 4.11 |

Mapped is three event types and six kernels. The smallest cases reach the
roughly 0.01 s timer resolution; zero readings are not evidence of zero cost.
At 10,000 rows, state is about 1.5-1.6 times slower. Public wrappers add
per-event specification/input validation and R calls; no precise attribution
percentage is claimed. This work establishes the model boundary, not a speedup.

Recorded cumulative Rprofmem allocation bytes, one separate full call after
warm-up for each model/backend:

| Model | Bridge bytes | State bytes |
|---|---:|---:|
| univariate, 1,000 rows | 193,088 | 186,288 |
| mapped, 1,000 rows | 2,434,432 | 2,394,672 |

These records are not peak resident memory, which was not measured. Full
returned histories still scale with rows x event types x kernel columns.
No meaningful memory reduction is claimed.

A separate high-frequency 1 ms Rprof trial produced converted summary times
inconsistent with observed wall time, so it was excluded from timing and
attribution evidence. All timing rows above are unprofiled calls. A scratch
benchmark launcher path was corrected before measurement; it did not affect
package code or the retained benchmark results.

## Remaining work

Keep the independent loop opt-in until prepared common-kernel execution or
another measured change removes its overhead. Event callbacks and general
marks still use bridge. Fitting still uses public emhawkes parameter mapping,
optimization and reconstruction. No residual slots or saved-object formats
were removed, and no third engine package was introduced.


## Prepared-plan follow-up, 2026-10-19

The development simulator now reuses exported `exp_kernel_plan(kernel, mu)`.
The locked plan captures numeric kernel parameters, mappings, row rates and
baseline once. Its step/wait methods still validate changing components and
arguments. Current state remains an explicit argument; the plan stores no
event state, histories or RNG. Existing public step/wait functions retain their
checked contracts. Older event-wise installations use the original checked
calls, and releases without those APIs retain bridge support.

For 1,000 output rows, kernel construction falls from 2,000 calls to two.
Residual refresh, inversion and root tolerance, row-sum order, decay and jump
arithmetic, event selection, counts, output attributes and seed restoration
remain unchanged. Bridge remains the default. This production implementation
retains dynamic component checks in every public method.

### Exact-result checks

- Existing 72 paired settings pass with complete paths, final RNG, following
  draws and independent R/C++ driver reconstruction exactly equal to bridge.
- An additional 108 settings span 1/3/8 types, 3/9/24 kernel columns, three laws,
  three seeds and 1/2/101/2,000 rows. All outputs and RNG match exactly.
  Independent R/C++ reconstruction in 18 comparisons is exact for all four
  driver/component fields.
- Twelve invalid-residual/baseline records match the previous installed state
  implementation, including errors, warnings, RNG and following draws.
- Prepared/public API comparisons cover 864 marked/unmarked steps, signed
  states, overlapping mappings, tiny and zero intervals, optional integrals,
  zero baseline and infinite targets. Captured inputs, independent repeated
  calls, serialization, input failures and no-RNG behavior pass.
- Tests explicitly enforce the bridge default, constant preparation count and
  compatibility with the preceding event-wise API. All ten FRPP test files also
  pass against that actual earlier installation, including the 72 state cases.

### Same-machine timings

Baseline: emhawkes `5e96cde` and frpp `32e5713`. After: the checked source
snapshot for this change, recorded in source.json and validation hashes.
Windows, R 4.5.1; whole fsim calls including preparation and seed handling;
installation, compilation and R process startup excluded. Separate persistent
R workers load the before/after installations; only one executes at a time.
Each model/size/backend gets a full warm-up, then five alternating/reordered
repetitions. The 32-row cases use 40 calls per repetition. No other heavy local
job ran during timing. All timed paths match each worker's reference exactly;
the cross-backend equality checks above run independently.

| Model | Rows | Bridge seconds | State before | State after | Time reduction |
|---|---:|---:|---:|---:|---:|
| univariate | 32 | 0.00450 | 0.00600 | 0.00475 | 20.8% |
| univariate | 1,000 | 0.11000 | 0.17000 | 0.12000 | 29.4% |
| univariate | 10,000 | 1.04000 | 1.56000 | 1.18000 | 24.4% |
| mapped | 32 | 0.00750 | 0.01025 | 0.00775 | 24.4% |
| mapped | 1,000 | 0.20000 | 0.30000 | 0.22000 | 26.7% |
| mapped | 10,000 | 2.02000 | 2.83000 | 2.19000 | 22.6% |

The prepared loop remains opt-in. Smaller-case per-call differences are at
millisecond scale; timer resolution and normal run variation limit precision.

Recorded cumulative Rprofmem bytes in separate warmed 1,000-row calls:

| Model | Bridge bytes | State before | State after |
|---|---:|---:|---:|
| univariate | 193,088 | 186,288 | 187,288 |
| mapped | 2,434,432 | 2,394,672 | 2,395,672 |

These are recorded cumulative allocations, not peak resident memory. Peak
memory was not measured. Full output histories retain their existing scaling.

### Local package checks

- emhawkes: 49 test files and vignette rebuild, R CMD check Status: OK.
- frpp with development emhawkes: ten test files and both vignettes, Status: OK.
- frpp with CRAN emhawkes 1.0.1: ten test files and both vignettes, Status: OK.
- Quick development suite: 13/13 files; CRAN and earlier-API suites: 10/10 each.
- Updated shared-kernel and FRPP computation guides render successfully.

Checks use --no-manual --run-donttest and _R_CHECK_SUGGESTS_ONLY_=false.
The PDF reference manual was not compiled. Unrelated emhawkes working-tree
changes were preserved and excluded from the tested source snapshot.
