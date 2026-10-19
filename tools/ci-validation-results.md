# Cross-package CI validation, 2026-10-19

No package calculation code, public API or dependency requirement changed.
The FRPP working tree and index were clean before this task. Existing emhawkes
work was left untouched.

## Local results

Windows 11, R 4.5.1, maxLik 1.5-2.1. Separate libraries contain the actual CRAN
emhawkes 1.0.1 source and the published development API at emhawkes cde77c9
(1.0.1.9000). CRAN source MD5: `7c8d439adceab2a4742de3086c05585d`.

| Check | CRAN dependency | Development dependency |
|---|---|---|
| Quick suite | 8/8 files | 9/9 files |
| Full command elapsed | 10.92 s | 18.15 s |
| R CMD check, both vignettes rebuilt | Status: OK | Status: OK |
| Optimizers | stats and maxLik | stats and maxLik |

Exact simulation fields, final RNG and paired fitting coefficients, likelihood,
raw Hessian, covariance and initial state checks passed. Independent state
comparisons keep 1e-10 and public-kernel arithmetic checks keep 1e-12; no
tolerance was widened. These elapsed times are single validation runs,
not a performance benchmark.

A development run against the CRAN library failed before testing, as required.
In isolated memory sessions, an extra simulation RNG draw and a 1e-5 interval
integral perturbation each failed the corresponding assertion. No installed
package files were modified by these fault checks.

## Compatibility corrections

Existing tests assumed the optional development API was always installed.
The release profile now exercises the bridge and verifies the explicit state
capability error without RNG consumption. The bridge exposes no interval
integrals and cannot be forced through the C++ Hawkes likelihood engine for a
residual-law model. Those development-only checks are required in development
jobs; release intensities and likelihood still use independent all-past-event
and quadrature references. An intermediate R test-output syntax error was
caught and corrected before the final passing runs.

YAML structure and all workflow R blocks parsed successfully. The new workflow
has six OS/dependency combinations, requires both optimizers and retains
quick reports and dependency provenance even when a later check fails.
See [validation instructions](validation.md) for trigger and coverage limits.

Cross-platform results are reported by the Package compatibility workflow on
GitHub; the local results above do not imply that remote jobs have passed.

The first macOS CI logs reported the Node.js 20 deprecation warning from the
copied v4 checkout/upload actions. The new compatibility workflow now uses
verified Node.js 24 releases: checkout v7.0.1 and upload-artifact v7.0.2. Their
input contracts were checked before the update. Package and test code did not
change in this follow-up; the final CI run validates the updated actions.
