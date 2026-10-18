# Flexible Residual Point Processes

**frpp** combines exponential kernels from
[emhawkes](https://ksublee.github.io/emhawkes/) with a user-specified residual
distribution. It provides model specification, simulation, likelihood evaluation
and maximum likelihood estimation.

This is a development package. Its model object contains a kernel and a residual
law; it does not inherit from the Hawkes specification.

## Start here

- [Get started](articles/get-started.html): installation, a first model and simulation.
- [Examples](articles/frpp-examples.html): one- and two-dimensional models migrated from emhawkes.
- [Legacy emhawkes examples](articles/emhawkes-compatibility.html): the complete original FRPP examples and fitting code.
- [Fitting and inference](articles/fitting-and-inference.html): optimizers, identification, diagnostics and restoration.
- [10,000-row fit validation](articles/legacy-fit-validation.html): matched original and migrated estimates, curvature and elapsed times.
- [Model and computation](articles/model-and-computation.html): residual refresh rules, kernel mapping and the emhawkes boundary.
- [Function reference](reference/index.html): arguments and return values.

## Main functions

| Task | Function |
|---|---|
| Specify a residual law | `residual_law()` |
| Compose a model | `fspec()` |
| Simulate a path | `fsim()` or `stats::simulate()` |
| Evaluate the likelihood | [`logLik()`](reference/operations.html) |
| Fit parameters | `ffit()` |
| Inspect a fit | `summary()`, `coef()`, `vcov()`, `fdiagnostics()` |
| Restore a fitted model | `fitted_fspec()` and `finitial_state()` |

Simulation redraws every dimension's candidate residual after every event
in the combined process. This is a model rule, not an independent renewal clock
for each dimension.

## Current implementation

The default likelihood backend and simulator use public emhawkes interfaces.
The experimental `state` backend assembles the FRPP likelihood from a shared
kernel replay API available in GitHub development emhawkes 1.0.1.9000.
Examples detect whether it is installed and label unavailable comparisons.

Residual parameters remain fixed during fitting unless joint estimation is
requested. A successful optimizer or a positive covariance matrix does not
establish valid confidence-interval coverage. The guides explain these limits.
