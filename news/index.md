# Changelog

## frpp 0.0.0.9000

- Add a quick joint regression runner with a fixed small fitting
  fixture, isolated test sessions and source/installed-code provenance.
  Keep saved-fit replay and expensive full re-estimation as separate
  validation steps.

- Convert `fspec` to an independent S4 composition class with typed
  slots, semantic validation, component accessors and S3/S4 method
  dispatch. Upgrade earlier S3 specifications and fitted models while
  retaining list-style component access, callback bindings, calculation
  and RNG behavior.

- Verify the 10,000-row migrated fit against the original emhawkes
  interface in three paired repetitions; document exact estimates and
  optimizer curvature, restored likelihoods, unavailable Wald covariance
  and elapsed fitting times.

- Preserve the full original emhawkes FRPP examples in a compatibility
  vignette, alongside the new-interface examples. Document the planned
  migration while retaining simple residual-law uses and Hawkes residual
  diagnostics in emhawkes.

- Add [`fsim()`](https://ksublee.github.io/frpp/reference/operations.md)
  as the primary simulation function. Existing
  [`simulate()`](https://rdrr.io/r/stats/simulate.html) methods call it,
  preserving path and RNG behavior and old serialized classes.

- Assemble numeric state-backend likelihoods using the experimental
  public `emhawkes::exp_kernel_replay()` interface, without Hawkes
  likelihood evaluation or full component histories when inference is
  not requested. This requires the emhawkes development interface; the
  bridge remains available with 1.0.0.

- Use [`fspec()`](https://ksublee.github.io/frpp/reference/frpp_spec.md)
  as the primary composition constructor. Keep
  [`frpp_spec()`](https://ksublee.github.io/frpp/reference/frpp_spec.md)
  and old serialized classes compatible; add
  [`as_fspec()`](https://ksublee.github.io/frpp/reference/frpp_spec.md).

- Add [`ffit()`](https://ksublee.github.io/frpp/reference/ffit.md) with
  maxLik and stats optimizers, optional joint residual estimation,
  shared/fixed parameters, constraints and a domain validator.

- Add fitted-model restoration, initialization replay, coefficient,
  likelihood, covariance, summary and diagnostic methods. Residual
  parameters stay fixed by default; invalid Wald covariance is withheld
  rather than repaired.
