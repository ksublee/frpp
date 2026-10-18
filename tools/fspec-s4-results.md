# S4 fspec validation

2026-10-18, Windows / R 4.5.1, maxLik 1.5-2.1, frpp 0.0.0.9000.
The pre-edit package was installed separately before the source changed.
Package libraries, source hashes, installed-code hashes and logs were retained
with the validation artifacts.

## Representation and compatibility

`fspec` is now an independent S4 composition class with four typed slots:
`kernel` (hspec), `residual` (the existing S3 residual law), `refresh` and
representation `version = 2L`. It does not inherit hspec. The dependency remains
frpp -> emhawkes. Baseline still belongs to the embedded hspec; extracting a
separate baseline/common-kernel representation is deferred.

`fkernel()` and `fresidual()` read and replace model components. Replacement
validates the complete specification. Existing `$kernel`/`$residual` access,
`fspec()`, `frpp_spec()`, `fsim()`, `simulate()`, `logLik()` and `ffit()` remain
available. Both historical S3 specification classes (version 1L) convert via
`as_fspec()`; public operations and fitted-model extraction also upgrade them.
Use `methods::is(model, "fspec")` for new class checks. Arbitrary list operations
and changing `class(model)` are not an S4 migration interface.

Validity checks cover stored numeric dimensions, finite values, positive beta,
nonnegative baseline, mappings, residual structure and global refresh semantics.
They do not execute callbacks, evaluate callback defaults, establish
stationarity or validate confidence-interval coverage. Functional outputs and
residual-specific parameter domains still need model-specific validation.
The fitting result remains S3; its extracted model is S4.

## Results

| Check | Result |
|---|---|
| Before/after path and likelihood snapshots | 40/40 exactly equal, including all common path fields, inference histories, final RNG and following draws; model representation excluded |
| Before/after fits | Both stats and maxLik: coefficients, likelihood attributes, raw Hessian, covariance, initial state, diagnostics and RNG exactly equal |
| Actual pre-edit saved specifications | 10 S3 models converted with kernel and residual objects unchanged |
| Actual pre-edit saved fitted result | Extracted as valid S4; restored likelihood exactly equal |
| Minimal-session tests | New and both legacy S3 model classes restored without attaching emhawkes or methods; stats fits also restore without loading maxLik |
| Fast suite | 11/11 files; loop 14.86 s, full command 15.58 s including parent R startup/report preparation |
| Full package check | R CMD check --no-manual, default `_R_CHECK_SUGGESTS_ONLY_=false`: Status OK; eight FRPP test files and both vignettes/rebuild included |
| Release-source emhawkes 1.0.1 | Simulation/RNG, bridge likelihood, both fitting providers and explicit missing-state-API error passed |
| Saved 10,000-row fits | All six restored and independently replayed in 28.16 s without reoptimization; original reference artifacts preserved |

The quick suite includes 72 exact R/C++ state comparisons and independent
all-past-event/quadrature oracles. The additional S4 test covers slot types,
semantic validation, replacement, invalid inputs before RNG use, old-class
conversion, preserved callback identity/defaults and both S3/S4 dispatch paths.
Before the change, the ten-file quick loop took 13.91 s on this machine. These
single-run elapsed observations have different test counts and are not a
performance benchmark. S4 conversion is a representation/validation change.

The six large saved likelihoods remain -12411.1075525012. Independent replay
still differs by 1.12777343019843e-10 in absolute log likelihood. Repeating both
summation orders with the same component values reproduces both totals exactly.
The existing relative 1e-10 replay comparison was retained; no tolerance was
increased. Curvature limitations of those fits are unchanged.

## Issues found during validation

- An old S3 object named fspec can satisfy `methods::is(x, "fspec")`. Checking
  `isS4(x)` as well prevents treating its list fields as S4 slots.
- A minimal session needs an explicit import of `methods::show` for the exported
  show method. Added the import and kept the clean-session regression checks.
- Comparing separately generated compatibility callback closures with
  `identical()` also compares their newly allocated environments. Conversion
  instead checks the original residual functions exactly and checks generated
  adapters through output, likelihood and RNG equivalence. Bound legacy
  callback objects still round-trip exactly.
- In the release-source Poisson smoke test, default numerical optimization
  returned mu = 3.0000004074078754 rather than the analytic 3. This was exactly
  the same before and after S4 conversion. Both providers returned exactly 3
  with the known analytic derivatives. The smoke check uses those derivatives;
  the paired numerical-optimizer tests remain in the quick suite. No tolerance
  was widened and no numerical optimizer behavior was changed.

The checks above were completed locally before publication. Platform CI is
separate from these local results. The PDF manual was not compiled
(`--no-manual`). The existing hspec/engine/event-loop
implementation and user callback binding rules were not changed.
